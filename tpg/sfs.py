#!/usr/bin/env python3

#=============================
#        Marek Gorski
#=============================

"""Skyflat selector (sfs) - picks the sky flat field for the evening or morning session.

The session runs while the Sun is between `skyflat.hsun` limits (e.g. +1 -> -15 deg in the
evening, -15 -> +1 in the morning), on a 1-minute grid like tpg's night. Each field from the
fields file gets a per-minute `visibility` array narrowed step by step (altitude, Sun azimuth,
Moon, wind), merged under "all". A field must stay visible for the whole session; among those
the one closest on average to the twilight sky's gradient null point wins - the point ~13 deg
from zenith on the anti-solar vertical (Chromey & Hasselbacher 1996, PASP 108, 944), where the
sky brightness gradient is smallest. If no field stays visible the whole session, the one
visible for the longest part of it is taken and the result is marked as `fallback`.

Usage:
    poetry run sfs zb08 2026/10/26 [--evening | --morning] [--wind DEG] [--loud]
"""

import argparse
import copy
import sys

import ephem
import numpy

from .config import load_config, resolve_path, deep_update
from .dome_eq import dome_eq_azimuth

# zenith distance of the twilight sky gradient null point, on the anti-solar vertical
NULL_POINT_ZD = 13.0

SESSIONS = ("evening", "morning")


class NoSkyflatField(ValueError):
    """No field is visible at any minute of the session - the CLI reports it without a traceback."""


class SkyflatSelector():
    def __init__(self, tel, date, session, wind=None, loud=False, fields_file=None, cfg=None):
        if session not in SESSIONS:
            raise ValueError(f'session must be one of {SESSIONS}, not {session!r}')
        self.tel = tel
        self.date = str(date).replace("-", "/")
        self.session = session
        self.wind = wind
        self.loud = loud
        self.fields_file = fields_file
        self.cfg = cfg              # None = load config.yaml + local config in Initiate

        self.msg = ""
        self.fields = []
        self.sessionTime = []       # minute grid of the session, ephem.Date
        self.result = None

    def go(self):
        self.Initiate()
        self.LoadFields()
        self.MakeTime()
        self.CalcFields()
        self.MaskAltitude()
        self.MaskSun()
        self.MaskMoon()
        self.MaskWind()
        self.Rank()
        return self.result

    def log_msg(self, txt, always_loud=True):
        self.msg = self.msg + txt + "\n"
        if self.loud or always_loud:
            print(txt)

    def Initiate(self):
        if self.cfg is None:
            self.cfg = load_config(log=lambda txt: self.log_msg(txt, always_loud=False))

        # skyflat block from the top level, overridden field by field by the telescope's own
        self.scfg = copy.deepcopy(self.cfg.get("skyflat", {}))
        deep_update(self.scfg, self.cfg[self.tel].get("skyflat", {}))

        self.hmin = float(self.scfg["hmin"])
        self.hmax = float(self.scfg["hmax"])
        self.sun_az_dist = float(self.scfg["sun_az_dist"])
        self.moon_dist = float(self.scfg.get("moon_dist", self.cfg[self.tel]["moon_dist"]))
        self.wind_dist = float(self.scfg.get("wind_dist", self.cfg[self.tel]["wind_dist"]))
        self.hsun_low, self.hsun_high = sorted(float(h) for h in self.scfg["hsun"])

    def LoadFields(self):
        fname = self.fields_file or resolve_path(self.scfg["fields_file"])
        self.fields = load_fields(fname, log=self.log_msg)
        self.log_msg(f'{len(self.fields)} fields loaded from {fname}', always_loud=False)

    def MakeTime(self):
        self.oca = ephem.Observer()
        self.oca.lon = str(self.cfg["obs_lon"])
        self.oca.lat = str(self.cfg["obs_lat"])
        self.oca.elev = float(self.cfg["obs_elev"])

        # same convention as tpg: the date is the night starting that evening
        sun = ephem.Sun()
        self.oca.date = ephem.Date(self.date + " 12:00:00")
        self.oca.horizon = str(self.hsun_high)
        evening_high = self.oca.next_setting(sun, use_center=True)
        self.oca.horizon = str(self.hsun_low)
        evening_low = self.oca.next_setting(sun, use_center=True)

        if self.session == "evening":
            t0, t1 = evening_high, evening_low
        else:
            self.oca.date = evening_low
            self.oca.horizon = str(self.hsun_low)
            t0 = self.oca.next_rising(sun, use_center=True)
            self.oca.horizon = str(self.hsun_high)
            t1 = self.oca.next_rising(sun, use_center=True)

        n = int((t1 - t0) / ephem.minute) + 1
        self.sessionTime = [ephem.Date(t0 + i * ephem.minute) for i in range(n)]
        self.start_time, self.end_time = t0, t1

        txt = f'{self.session} session {ephem.Date(t0)} - {ephem.Date(t1)} UT ({n} min), ' \
              f'Sun {self.hsun_high:+.0f} / {self.hsun_low:+.0f} deg'
        self.log_msg(txt, always_loud=False)

    def CalcFields(self):
        sun = ephem.Sun()
        moon = ephem.Moon()

        sun_alt, sun_az = [], []
        for x in self.sessionTime:
            self.oca.date = x
            sun.compute(self.oca)
            sun_alt.append(numpy.degrees(float(sun.alt)))
            sun_az.append(numpy.degrees(float(sun.az)))
        self.sun = {"alt": numpy.array(sun_alt), "az": numpy.array(sun_az)}

        dome = self.cfg[self.tel].get("dome")

        for field in self.fields:
            s = ephem.FixedBody()
            s._ra = field["ra"]
            s._dec = field["dec"]

            alt, az, ha, moon_sep = [], [], [], []
            for x in self.sessionTime:
                self.oca.date = x
                s.compute(self.oca)
                moon.compute(self.oca)
                alt.append(numpy.degrees(float(s.alt)))
                az.append(numpy.degrees(float(s.az)))
                ha.append(numpy.degrees(float(s.ha)))
                moon_sep.append(numpy.degrees(float(ephem.separation(s, moon))))

            field["calc"] = {"alt": numpy.array(alt), "az": numpy.array(az),
                             "moon_sep": numpy.array(moon_sep)}

            # GEM: as in tpg, wind is checked against the dome azimuth for both sides of pier
            if dome and dome.get("mount") == "gem":
                field["calc"]["dome_az"] = [
                    dome_eq_azimuth(az, alt, ha, r_dome=dome["r_dome"], spx=dome["spx"],
                                    spy=dome["spy"], gem=dome["gem"], side_of_pier=side,
                                    latitude=self.cfg["obs_lat"])
                    for side in (0, 1)
                ]

            field["visibility"] = {"all": numpy.ones(len(self.sessionTime), dtype=bool)}

    def _apply_mask(self, field, key, mask):
        field["visibility"][key] = mask
        field["visibility"]["all"] = field["visibility"]["all"] & mask

    def MaskAltitude(self):
        for field in self.fields:
            alt = field["calc"]["alt"]
            self._apply_mask(field, "h_min", alt >= self.hmin)
            self._apply_mask(field, "h_max", alt <= self.hmax)

    def MaskSun(self):
        for field in self.fields:
            mask = az_dist(field["calc"]["az"], self.sun["az"]) >= self.sun_az_dist
            self._apply_mask(field, "sun_az", mask)

    def MaskMoon(self):
        for field in self.fields:
            self._apply_mask(field, "moon_dist", field["calc"]["moon_sep"] >= self.moon_dist)

    def MaskWind(self):
        if self.wind is None:
            return

        self.log_msg(f'WIND condition: {self.wind}', always_loud=False)

        for field in self.fields:
            az_list = field["calc"].get("dome_az", [field["calc"]["az"]])
            mask = numpy.ones(len(self.sessionTime), dtype=bool)
            for az in az_list:
                mask &= az_dist(az, self.wind) >= self.wind_dist
            self._apply_mask(field, "wind", mask)

    def Rank(self):
        null_alt = 90.0 - NULL_POINT_ZD
        null_az = (self.sun["az"] + 180.0) % 360.0

        for field in self.fields:
            dist = ang_dist(field["calc"]["alt"], field["calc"]["az"], null_alt, null_az)
            field["score"] = float(dist.mean())
            field["fraction"] = float(field["visibility"]["all"].mean())

        # visible for the whole session first, then by distance to the null point
        self.fields.sort(key=lambda f: (-f["fraction"], f["score"]))

        if self.loud:
            self._printFieldLog()

        best = self.fields[0] if self.fields else None
        if best is None or best["fraction"] == 0:
            raise NoSkyflatField(f'no sky flat field visible during the {self.session} session '
                                 f'of {self.date} on {self.tel}')

        self.result = {
            "name": best["name"],
            "ra": best["ra"],
            "dec": best["dec"],
            "line": best["line"],
            "session": self.session,
            "start": self.start_time,
            "end": self.end_time,
            "score": best["score"],
            "fraction": best["fraction"],
            "fallback": best["fraction"] < 1.0,
        }
        if self.result["fallback"]:
            self.log_msg(f'WARNING: no field visible for the whole session, '
                         f'{best["name"]} visible for {100 * best["fraction"]:.0f}%')

    def _printFieldLog(self):
        print(f'{"field":10} {"alt":>11} {"az":>11} {"sun_az":>6} {"moon":>5} {"vis":>5} {"score":>6}  rejected by')
        sun_az = self.sun["az"]
        for f in self.fields:
            c = f["calc"]
            failed = [k for k, v in f["visibility"].items() if k != "all" and not v.all()]
            print(f'{f["name"]:10} {c["alt"].min():5.1f}-{c["alt"].max():5.1f} '
                  f'{c["az"].min():5.1f}-{c["az"].max():5.1f} '
                  f'{az_dist(c["az"], sun_az).min():6.1f} {c["moon_sep"].min():5.1f} '
                  f'{100 * f["fraction"]:4.0f}% {f["score"]:6.1f}  {",".join(failed)}')


def load_fields(fname, log=print):
    """Read `SKYFLAT name ra dec [key=value ...]` lines; `#` comments and empty lines skipped."""
    fields = []
    with open(fname, 'r') as file:
        for line in file:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            tokens = line.split()
            if len(tokens) < 4 or tokens[0].upper() != "SKYFLAT":
                log(f'Line ERROR: {line}')
                continue
            try:
                ephem.hours(tokens[2])
                ephem.degrees(tokens[3])
            except ValueError:
                log(f'Line ERROR (coordinates): {line}')
                continue
            kwargs = dict(t.split("=", 1) for t in tokens[4:] if "=" in t)
            fields.append({"name": tokens[1], "ra": tokens[2], "dec": tokens[3],
                           "kwargs": kwargs, "line": line})
    return fields


def az_dist(az1, az2):
    """Azimuth difference folded to 0-180 deg."""
    d = numpy.abs(numpy.asarray(az1, dtype=float) - numpy.asarray(az2, dtype=float)) % 360.0
    return numpy.minimum(d, 360.0 - d)


def ang_dist(alt1, az1, alt2, az2):
    """Angular distance between two (alt, az) directions [deg]."""
    alt1, az1, alt2, az2 = (numpy.radians(numpy.asarray(a, dtype=float)) for a in (alt1, az1, alt2, az2))
    cos_d = numpy.sin(alt1) * numpy.sin(alt2) + numpy.cos(alt1) * numpy.cos(alt2) * numpy.cos(az1 - az2)
    return numpy.degrees(numpy.arccos(numpy.clip(cos_d, -1.0, 1.0)))


def main():
    argparser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    argparser.add_argument("tel", type=str, choices=['wk06', 'zb08', 'jk15', 'iris', 'test'], help='Telescope id')
    argparser.add_argument("date", type=str, help='yyyy/mm/dd date of the night (evening)')
    session = argparser.add_mutually_exclusive_group()
    session.add_argument("--evening", action="store_true", help='evening session only')
    session.add_argument("--morning", action="store_true", help='morning session only')
    argparser.add_argument("--wind", type=float, default=None, help='wind direction (to avoid)')
    argparser.add_argument("--fields_file", type=str, default=None, help='specify sky flat fields file')
    argparser.add_argument("--loud", action="store_true", help="print all fields with masks and scores")
    args = argparser.parse_args()

    if args.evening:
        sessions = ["evening"]
    elif args.morning:
        sessions = ["morning"]
    else:
        sessions = list(SESSIONS)

    for session in sessions:
        s = SkyflatSelector(args.tel, args.date, session, wind=args.wind, loud=args.loud,
                            fields_file=args.fields_file)
        try:
            r = s.go()
        except NoSkyflatField as e:
            sys.stderr.write(f"sfs: {e}\n")
            sys.exit(1)
        note = f', visible {100 * r["fraction"]:.0f}% of session' if r["fallback"] else ""
        print(f'# {session} {ephem.Date(r["start"])} - {ephem.Date(r["end"])} UT, '
              f'null point distance {r["score"]:.1f} deg{note}')
        print(r["line"])


if __name__ == "__main__":
    main()
