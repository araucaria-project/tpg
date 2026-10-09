"""Azymut kopuly dla montazu GEM (wk06) i maska wiatru liczona na nim."""

import numpy
import pytest
import yaml

import pyaraucaria.dome_eq as pa_dome_eq

from conftest import SITE, obj, run_cascade, vis
from tpg.dome_eq import dome_eq_azimuth
from tpg.telescope_plan_generator import TelescopePlanGenerator, deep_update

# geometria wk06, jak w tpg/config.yaml i TOI
DOME = {"mount": "gem", "r_dome": 2050, "spx": -110, "spy": -110, "gem": 670}


class TestGeometry:
    @pytest.mark.parametrize("side", [0, 1])
    def test_same_as_pyaraucaria(self, monkeypatch, side):
        """Lokalna kopia ma liczyc dokladnie to, co pyaraucaria - ma ja kiedys zastapic.

        pyaraucaria bierze czas "teraz", wiec podstawiamy jej az/alt i czas gwiazdowy;
        przy ra=0 kat godzinny to po prostu czas gwiazdowy.
        """
        rng = numpy.random.default_rng(1)
        for az, alt, ha in rng.uniform([0, 20, -180], [360, 90, 180], size=(200, 3)):
            monkeypatch.setattr(pa_dome_eq, "ra_dec_2_az_alt", lambda **kw: (az, alt))
            monkeypatch.setattr(pa_dome_eq, "site_sidereal_time", lambda **kw: ha)

            expected, _ = pa_dome_eq.dome_eq_azimuth(
                ra=0, dec=0, r_dome=2050, spx=-110, spy=-110, gem=670, side_of_pier=side,
                latitude=SITE["obs_lat"], longitude=SITE["obs_lon"], elevation=SITE["obs_elev"],
            )
            got = dome_eq_azimuth(az, alt, ha, r_dome=2050, spx=-110, spy=-110, gem=670,
                                  side_of_pier=side, latitude=SITE["obs_lat"])

            assert got == pytest.approx(expected % 360, abs=1e-9)

    def test_without_offsets_dome_follows_telescope(self):
        az = numpy.linspace(0, 359, 50)
        got = dome_eq_azimuth(az, 45.0, 30.0, r_dome=2050, spx=0, spy=0, gem=0,
                              side_of_pier=0, latitude=SITE["obs_lat"])

        assert got == pytest.approx(az)


class TestWindOnDome:
    def test_alt_az_telescope_has_no_dome_azimuth(self, make_gen):
        p = run_cascade(make_gen([obj()]), upto="CalcObject")

        assert "dome_az" not in p.ob[0]["calc"]

    def test_both_sides_of_pier_are_masked(self, make_gen):
        bazowy = run_cascade(make_gen([obj()], tel_cfg={"dome": DOME}), upto="MaskVisibility")
        dome_az = numpy.array(bazowy.ob[0]["calc"]["dome_az"])
        kierunek = float(dome_az[0][len(dome_az[0]) // 2])

        p = run_cascade(make_gen([obj()], wind=kierunek, tel_cfg={"dome": DOME}), upto="MaskWind")

        odleglosc = numpy.abs(kierunek - dome_az)
        odleglosc = numpy.where(odleglosc > 180, 360 - odleglosc, odleglosc)
        zly = (odleglosc < 95).any(axis=0)
        assert numpy.array_equal(vis(p, key="wind"), vis(bazowy) & ~zly)

    def test_dome_changes_the_mask(self, make_gen):
        """Kontrola, ze test powyzej nie przechodzi trywialnie - kopula != teleskop."""
        bazowy = run_cascade(make_gen([obj()], tel_cfg={"dome": DOME}), upto="CalcObject")
        az = numpy.array(bazowy.ob[0]["calc"]["az"])
        dome_az = numpy.array(bazowy.ob[0]["calc"]["dome_az"])

        roznica = numpy.abs((dome_az - az + 180) % 360 - 180)
        assert roznica.max() > 5


class TestConfigMerge:
    def test_deep_update_keeps_fields_missing_locally(self):
        base = {"wk06": {"wind_dist": 30, "hmin": 35, "dome": {"r_dome": 2050, "gem": 670}}}
        deep_update(base, {"wk06": {"hmin": 40, "dome": {"gem": 700}}, "plan_catalog": "/x"})

        assert base == {
            "wk06": {"wind_dist": 30, "hmin": 40, "dome": {"r_dome": 2050, "gem": 700}},
            "plan_catalog": "/x",
        }

    def test_local_config_overrides_single_fields(self, tmp_path, monkeypatch):
        local = tmp_path / "local_config.yaml"
        local.write_text(yaml.safe_dump({"wk06": {"master_file": "/somewhere/wk06.txt"}}))
        monkeypatch.setattr("os.path.expanduser", lambda path: str(local))

        p = TelescopePlanGenerator("wk06", ["2026/08/12"])
        p._get_config()

        assert p.cfg["wk06"]["master_file"] == "/somewhere/wk06.txt"
        # pola, ktorych lokalny plik nie podal, zostaja z tpg/config.yaml
        assert p.cfg["wk06"]["dome"]["mount"] == "gem"
        assert "wind_dist" in p.cfg["wk06"]
