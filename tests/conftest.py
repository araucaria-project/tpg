"""Wspolne fixture'y i pomocniki dla testow tpg.

Dwie zasady, ktore warto zachowac dopisujac testy:

1. Config jest odcinany od maszyny. `~/.config/tpg/local_config.yaml` istnieje na
   maszynach deweloperskich i nadpisuje cale bloki teleskopow, wiec test, ktory go
   nie odetnie, przechodzi u jednej osoby i pada u drugiej.

2. Tam gdzie oczekiwanie zalezy od efemeryd, wyprowadzamy je z `ob["calc"]` zamiast
   wpisywac magiczna liczbe. Testujemy logike masek, nie poprawnosc `ephem`.
"""

import ephem
import numpy
import pytest

from tpg.telescope_plan_generator import TelescopePlanGenerator

TEL = "test"
DATE = "2026/08/12"

# obiekt widoczny przez ~8 h tej nocy na tym stanowisku - domyslny "kon roboczy"
RA = "21:24:00.2"
DEC = "-30:16:43.7"

SITE = {
    "obs_lat": -24.598,
    "obs_lon": -70.1963,
    "obs_elev": 2817,
    # celowo nieistniejace - zaden test nie moze czytac ani pisac poza tmp_path
    "local_config": "/nonexistent/local_config.yaml",
    "plan_catalog": "/nonexistent/plans/",
}

# jak blok "test" w config.yaml, z roznicami: moon_dist=0, wind_dist=95 (config.yaml ma 40 od 2026-11).
# Maski musza byc izolowane - przy produkcyjnym 40 deg Ksiezyc z akurat tej daty
# zaklocalby testy, ktore z Ksiezycem nie maja nic wspolnego.
TEL_CFG = {
    "def_read_mode": 2,
    "overhead": 15,
    "hmin": 35,
    "hmax": 80,
    "hsun": -18,
    "wind_dist": 95,
    "moon_dist": 0,
}

# kolejnosc z go() - zmiana tutaj musi isc w parze ze zmiana w go()
CASCADE = (
    "ObjectMask",
    "CalcObject",
    "MaskVisibility",
    "MaskMoon",
    "MaskWind",
    "MaskTwilight",
    "MaskCycle",
    "MaskStartEnd",
    "MaskPhaseStartEnd",
    "MaskPhase",
)


def obj(name="target", ra=RA, dec=DEC, seq="1/Ic/30", priority=5, **kw):
    """Zloz linie master file. Dodatkowe atrybuty przez kwargs: obj(cycle=3)."""
    parts = [f"OBJECT {name} {ra} {dec} seq={seq} priority={priority}"]
    parts += [f"{k}={v}" for k, v in kw.items()]
    return " ".join(parts)


@pytest.fixture
def make_gen(tmp_path):
    """Zbuduj generator na podanych liniach master file, z configiem odcietym od maszyny."""

    def _make(lines, *, date=DATE, ut=None, tel_cfg=None, **kwargs):
        master = tmp_path / "master.txt"
        master.write_text("\n".join(lines) + "\n")

        arg_date = [date] if ut is None else [date, ut]
        kwargs.setdefault("save_plan", False)
        kwargs.setdefault("use_ctc", False)

        p = TelescopePlanGenerator(TEL, arg_date, master_file=str(master), **kwargs)
        p.Initiate()

        # Initiate() wczytal prawdziwy config.yaml i lokalny override - zastepujemy
        # go w calosci. Nic miedzy _get_config() a koncem Initiate() z cfg nie korzysta.
        cfg = dict(SITE)
        cfg[TEL] = dict(TEL_CFG, data_file=str(tmp_path / "data") + "/")
        if tel_cfg:
            cfg[TEL].update(tel_cfg)
        p.cfg = cfg
        return p

    return _make


@pytest.fixture
def write_lc(tmp_path):
    """Zapisz krzywa blasku tam, gdzie MaskCycle/MaskPhase jej szukaja."""

    def _write(name, band, jds, quality=None):
        d = tmp_path / "data" / name.lower() / band / "light-curve"
        d.mkdir(parents=True, exist_ok=True)
        f = d / f"{name.lower()}_{band}_diff_light_curve.txt"
        q = [0] * len(jds) if quality is None else quality
        rows = "\n".join(f"{jd:.6f} {qq}" for jd, qq in zip(jds, q))
        f.write_text(f"jd_obs quality\n{rows}\n")
        return f

    return _write


def run_cascade(p, upto=None):
    """LoadObjects + MakeTime + kaskada masek, opcjonalnie do wskazanej wlacznie."""
    p.LoadObjects()
    p.MakeTime()
    for step in CASCADE:
        getattr(p, step)()
        if step == upto:
            break
    return p


def vis(p, i=0, key="all"):
    """Tablica widocznosci obiektu i pod danym kluczem, jako bool."""
    return numpy.asarray(p.ob[i]["visibility"][key], dtype=bool)


def minutes(p, i=0, key="all"):
    """Ile minut nocy zostaje widocznych."""
    return int(vis(p, i, key).sum())


def night_jd(p):
    """Siatka minut nocy w JD - do budowania oczekiwan fazowych."""
    return numpy.array([ephem.julian_date(ephem.Date(t)) for t in p.nightTime])
