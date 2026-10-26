"""Testy skyflat selectora (sfs).

Te same zasady co w conftest.py: config podajemy w calosci (bez config.yaml i lokalnego
override'u), a oczekiwania zalezne od efemeryd wyprowadzamy z `field["calc"]` / `s.sun`.
"""

import ephem
import numpy
import pytest

from tpg.sfs import SkyflatSelector, NoSkyflatField, load_fields, az_dist, ang_dist, NULL_POINT_ZD

from conftest import SITE, TEL, TEL_CFG, DATE

SKYFLAT = {
    "hsun": [-15, 1],
    "hmin": 40,
    "hmax": 70,
    "sun_az_dist": 100,
}

# fragment master_files/sky_flat_fields.txt - pola rozlozone co ~2 h w RA
FIELDS = [
    "SKYFLAT BF00b   00:43:50.0  -29:24:00.0",
    "SKYFLAT BF02b   02:41:10.0  -32:07:00.0",
    "SKYFLAT BF04b   04:23:50.0  -30:11:00.0",
    "SKYFLAT BF06b   06:34:12.0  -28:14:00.0",
    "SKYFLAT BF08b   08:39:14.0  -27:58:00.0",
    "SKYFLAT BF10b   10:21:16.0  -27:04:00.0",
    "SKYFLAT BF12b   12:29:30.0  -27:55:00.0",
    "SKYFLAT BF14b   14:27:10.0  -31:20:00.0",
    "SKYFLAT BF16b   16:38:34.0  -32:34:00.0",
    "SKYFLAT BF18b   18:37:22.0  -36:47:30.0",
    "SKYFLAT BF20b   20:36:38.0  -29:17:00.0",
    "SKYFLAT BF22b   22:23:50.0  -30:58:00.0",
]

# nigdy nie wschodzi wyzej niz ~5 deg na OCM
NEVER_UP = "SKYFLAT north 12:00:00.0 +60:00:00.0"


@pytest.fixture
def make_sfs(tmp_path):
    """Zbuduj selector na podanych liniach pliku pol, z configiem odcietym od maszyny."""

    def _make(lines=FIELDS, *, session="evening", date=DATE, skyflat=None, tel_cfg=None, **kwargs):
        fields = tmp_path / "fields.txt"
        fields.write_text("\n".join(lines) + "\n")

        cfg = dict(SITE)
        cfg["skyflat"] = dict(SKYFLAT, fields_file=str(fields))
        cfg[TEL] = dict(TEL_CFG)
        if skyflat:
            cfg[TEL]["skyflat"] = skyflat
        if tel_cfg:
            cfg[TEL].update(tel_cfg)
        return SkyflatSelector(TEL, date, session, cfg=cfg, **kwargs)

    return _make


def sun_alt(s, t):
    sun = ephem.Sun()
    s.oca.date = t
    sun.compute(s.oca)
    return numpy.degrees(float(sun.alt))


class TestFields:
    def test_comments_and_bad_lines_are_skipped(self, tmp_path):
        f = tmp_path / "fields.txt"
        f.write_text("# komentarz\n\n"
                     "SKYFLAT ok 01:00:00 -30:00:00 h_min=50\n"
                     "OBJECT not_a_flat 01:00:00 -30:00:00\n"
                     "SKYFLAT short 01:00:00\n"
                     "SKYFLAT bad xx:00:00 -30:00:00\n")
        errors = []

        fields = load_fields(str(f), log=errors.append)

        assert [x["name"] for x in fields] == ["ok"]
        assert fields[0]["kwargs"] == {"h_min": "50"}
        assert fields[0]["line"] == "SKYFLAT ok 01:00:00 -30:00:00 h_min=50"
        assert len(errors) == 3


class TestSession:
    @pytest.mark.parametrize("session, first, last", [("evening", 1, -15), ("morning", -15, 1)])
    def test_session_spans_sun_altitude_range(self, make_sfs, session, first, last):
        s = make_sfs(session=session)
        s.Initiate()
        s.MakeTime()

        assert sun_alt(s, s.start_time) == pytest.approx(first, abs=0.1)
        assert sun_alt(s, s.end_time) == pytest.approx(last, abs=0.1)
        assert s.sessionTime[0] == s.start_time
        assert s.sessionTime[-1] <= s.end_time
        # 16 deg Slonca to na OCM okolo godziny, zawsze ponizej dwoch
        assert 50 < len(s.sessionTime) < 120

    def test_morning_follows_evening(self, make_sfs):
        ev = make_sfs(session="evening")
        mo = make_sfs(session="morning")
        for s in (ev, mo):
            s.Initiate()
            s.MakeTime()

        assert ev.end_time < mo.start_time < ev.end_time + 1

    def test_telescope_block_overrides_skyflat_fields(self, make_sfs):
        s = make_sfs(skyflat={"hmin": 30})
        s.Initiate()

        assert s.hmin == 30
        assert s.hmax == 70


class TestMasks:
    def run(self, s):
        s.Initiate()
        s.LoadFields()
        s.MakeTime()
        s.CalcFields()
        s.MaskAltitude()
        s.MaskSun()
        s.MaskMoon()
        s.MaskWind()
        return s

    def test_altitude_and_sun_masks_follow_calc(self, make_sfs):
        s = self.run(make_sfs())

        for f in s.fields:
            c = f["calc"]
            numpy.testing.assert_array_equal(f["visibility"]["h_min"], c["alt"] >= 40)
            numpy.testing.assert_array_equal(f["visibility"]["h_max"], c["alt"] <= 70)
            numpy.testing.assert_array_equal(f["visibility"]["sun_az"], az_dist(c["az"], s.sun["az"]) >= 100)

    def test_all_is_and_of_masks(self, make_sfs):
        s = self.run(make_sfs(wind=0.0))

        for f in s.fields:
            v = f["visibility"]
            expected = v["h_min"] & v["h_max"] & v["sun_az"] & v["moon_dist"] & v["wind"]
            numpy.testing.assert_array_equal(v["all"], expected)

    def test_moon_mask_uses_telescope_moon_dist(self, make_sfs):
        s = self.run(make_sfs(tel_cfg={"moon_dist": 60}))

        for f in s.fields:
            numpy.testing.assert_array_equal(f["visibility"]["moon_dist"], f["calc"]["moon_sep"] >= 60)

    def test_no_wind_no_wind_mask(self, make_sfs):
        s = self.run(make_sfs())

        assert all("wind" not in f["visibility"] for f in s.fields)

    def test_wind_masks_field_azimuth(self, make_sfs):
        s = self.run(make_sfs(wind=0.0))
        f = s.fields[0]
        az = f["calc"]["az"]
        wind = float(az[len(az) // 2])

        s = self.run(make_sfs(wind=wind))
        f = s.fields[0]

        numpy.testing.assert_array_equal(f["visibility"]["wind"], az_dist(f["calc"]["az"], wind) >= 95)
        assert not f["visibility"]["wind"][len(az) // 2]

    def test_gem_dome_azimuth_for_both_sides(self, make_sfs):
        dome = {"mount": "gem", "r_dome": 2050, "spx": -110, "spy": -110, "gem": 670}
        s = self.run(make_sfs(tel_cfg={"dome": dome}, wind=0.0))

        for f in s.fields:
            assert len(f["calc"]["dome_az"]) == 2
            expected = numpy.ones(len(s.sessionTime), dtype=bool)
            for az in f["calc"]["dome_az"]:
                expected &= az_dist(az, 0.0) >= 95
            numpy.testing.assert_array_equal(f["visibility"]["wind"], expected)


class TestRank:
    @pytest.mark.parametrize("session", ["evening", "morning"])
    def test_best_fully_visible_field_closest_to_null_point(self, make_sfs, session):
        s = make_sfs(session=session)
        r = s.go()

        full = [f for f in s.fields if f["visibility"]["all"].all()]
        assert full, "test assumes at least one field is visible the whole session"
        best = min(full, key=lambda f: f["score"])

        assert r["name"] == best["name"]
        assert r["fallback"] is False
        assert r["fraction"] == 1.0
        assert r["line"].split()[1] == r["name"]

        c = best["calc"]
        expected = ang_dist(c["alt"], c["az"], 90 - NULL_POINT_ZD, (s.sun["az"] + 180) % 360).mean()
        assert r["score"] == pytest.approx(expected)

    def test_field_is_opposite_the_sun(self, make_sfs):
        for session in ("evening", "morning"):
            s = make_sfs(session=session)
            s.go()
            f = s.fields[0]

            assert az_dist(f["calc"]["az"], s.sun["az"]).min() >= 100

    def test_fallback_to_longest_visible(self, make_sfs):
        # 2 deg pasma wysokosci - zadne pole nie utrzyma sie w nim przez cala sesje
        s = make_sfs(skyflat={"hmin": 55, "hmax": 57})
        r = s.go()

        assert r["fallback"] is True
        assert 0 < r["fraction"] < 1
        assert r["fraction"] == max(f["fraction"] for f in s.fields)

    def test_nothing_visible_raises(self, make_sfs):
        s = make_sfs([NEVER_UP])

        with pytest.raises(NoSkyflatField):
            s.go()


class TestGeometry:
    def test_az_dist_wraps(self):
        assert az_dist(350, 10) == pytest.approx(20)
        assert az_dist(10, 350) == pytest.approx(20)
        assert az_dist(0, 180) == pytest.approx(180)

    def test_ang_dist(self):
        assert ang_dist(90, 0, 90, 123) == pytest.approx(0, abs=1e-6)
        assert ang_dist(0, 0, 0, 90) == pytest.approx(90)
        assert ang_dist(77, 10, 70, 10) == pytest.approx(7)
