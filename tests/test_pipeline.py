"""Wczytywanie master file, siatka czasu, liczenie dlugosci slotu."""

import ephem
import pytest

from conftest import TEL, obj, run_cascade


class TestTimeGrid:
    def test_resolution_is_one_minute(self, make_gen):
        p = make_gen([obj()])
        p.LoadObjects()
        p.MakeTime()

        steps = [b - a for a, b in zip(p.nightTime, p.nightTime[1:])]
        assert all(abs(s - ephem.minute) < 1e-9 for s in steps)

    def test_grid_spans_the_night(self, make_gen):
        p = make_gen([obj()])
        p.LoadObjects()
        p.MakeTime()

        # noc w Las Campanas w sierpniu - rzedu 10 h
        assert 9 * 60 < len(p.nightTime) < 12 * 60
        assert p.nightTime[0] > p.sunset
        assert p.nightTime[-1] > p.start_time

    def test_grid_edges_are_off_by_one_minute(self, make_gen):
        """Charakteryzacja, nie zyczenie.

        MakeTime() inkrementuje przed dopisaniem, wiec siatka zaczyna sie minute po
        starcie planu (pierwsza minuta przepada) i konczy do minuty za wschodem.
        Operacyjnie bez znaczenia, ale test ma zlapac moment, w ktorym ktos to zmieni.
        """
        p = make_gen([obj()])
        p.LoadObjects()
        p.MakeTime()

        assert p.nightTime[0] == pytest.approx(p.start_time + ephem.minute)
        assert p.sunrise < p.nightTime[-1] <= p.sunrise + ephem.minute

    def test_no_ut_starts_at_sunset(self, make_gen):
        p = make_gen([obj()])
        p.LoadObjects()
        p.MakeTime()

        assert "start_at_sunset" in p.flags
        assert p.start_time == p.sunset

    def test_explicit_ut_starts_there(self, make_gen):
        p = make_gen([obj()], ut="03:00:00")
        p.LoadObjects()
        p.MakeTime()

        assert "start_at_sunset" not in p.flags
        assert p.start_time == ephem.Date("2026/08/12 03:00:00")

    def test_grids_are_parallel(self, make_gen):
        """nightObject/nightFree/nightObChange musza miec dlugosc siatki minut."""
        p = make_gen([obj()])
        p.LoadObjects()
        p.MakeTime()

        n = len(p.nightTime)
        assert len(p.nightObject) == n
        assert len(p.nightFree) == n
        assert len(p.nightObChange) == n
        assert all(p.nightFree)


class TestFixtureIsolation:
    def test_local_config_does_not_leak(self, make_gen):
        """Jesli ten test padnie, wszystkie pozostale sa warte tyle co maszyna, na
        ktorej biegna. Initiate() czyta ~/.config/tpg/local_config.yaml - fixture
        musi go nadpisac."""
        p = make_gen([obj()])

        assert p.cfg[TEL]["moon_dist"] == 0
        assert p.cfg[TEL]["hmin"] == 35
        assert p.cfg["local_config"] == "/nonexistent/local_config.yaml"

    def test_data_file_points_into_tmp(self, make_gen, tmp_path):
        p = make_gen([obj()])

        assert p.cfg[TEL]["data_file"].startswith(str(tmp_path))


class TestLoading:
    def test_valid_line_is_loaded(self, make_gen):
        p = make_gen([obj(name="alfa")])
        p.LoadObjects()

        assert len(p.ob) == 1
        assert p.ob[0]["name"] == "alfa"

    def test_comments_and_blanks_are_skipped(self, make_gen):
        p = make_gen(["# komentarz", "", "   ", obj(name="alfa")])
        p.LoadObjects()

        assert [o["name"] for o in p.ob] == ["alfa"]

    def test_unparseable_line_does_not_kill_the_run(self, make_gen):
        p = make_gen(["OBJECT polamane bez wspolrzednych", obj(name="alfa")])
        p.LoadObjects()

        assert [o["name"] for o in p.ob] == ["alfa"]

    def test_hjd0_defaults_when_period_given(self, make_gen):
        p = make_gen([obj(P=1.5)])
        p.LoadObjects()

        assert float(p.ob[0]["hjd0"]) == 2460000

    def test_hjd0_is_not_invented_without_period(self, make_gen):
        p = make_gen([obj()])
        p.LoadObjects()

        assert "hjd0" not in p.ob[0]


class TestSlotTime:
    def test_computed_from_seq(self, make_gen):
        p = make_gen([obj(seq="10/Ic/60")])
        p.LoadObjects()

        # 10 x 60 s ekspozycji + narzuty; w minutach, wiec wyraznie ponad 10
        assert p.ob[0]["slotTime"] > 10

    def test_longer_sequence_gives_longer_slot(self, make_gen):
        p = make_gen([obj(name="krotki", seq="1/Ic/30"), obj(name="dlugi", seq="20/Ic/30")])
        p.LoadObjects()

        krotki, dlugi = p.ob
        assert dlugi["slotTime"] > krotki["slotTime"]

    def test_ob_time_overrides_seq(self, make_gen):
        p = make_gen([obj(seq="20/Ic/300", ob_time=7)])
        p.LoadObjects()

        assert float(p.ob[0]["slotTime"]) == 7

    def test_overhead_from_config_is_used(self, make_gen):
        """Ten sam seq przy wiekszym overhead musi dac dluzszy slot."""
        maly = make_gen([obj(seq="5/Ic/30")], tel_cfg={"overhead": 15})
        maly.LoadObjects()
        duzy = make_gen([obj(seq="5/Ic/30")], tel_cfg={"overhead": 60})
        duzy.LoadObjects()

        assert duzy.ob[0]["slotTime"] > maly.ob[0]["slotTime"]


class TestVisibilityBookkeeping:
    def test_every_object_gets_an_all_key(self, make_gen):
        p = run_cascade(make_gen([obj(name="a"), obj(name="b")]))

        for o in p.ob:
            assert "all" in o["visibility"]
            assert len(o["visibility"]["all"]) == len(p.nightTime)

    def test_cascade_only_narrows(self, make_gen):
        """Zadny krok kaskady nie moze dolozyc widocznej minuty."""
        import numpy

        from conftest import CASCADE

        p = make_gen([obj(t_start="03:00:00", sunset_dt=1)])
        p.LoadObjects()
        p.MakeTime()

        poprzednie = None
        for step in CASCADE:
            getattr(p, step)()
            biezace = numpy.asarray(p.ob[0]["visibility"]["all"], dtype=bool)
            if poprzednie is not None:
                dodane = biezace & ~poprzednie
                assert not dodane.any(), f"{step} poszerzyl widocznosc"
            poprzednie = biezace
