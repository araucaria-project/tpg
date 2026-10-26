"""Kaskada masek poza fazowaniem: selekcja, geometria, czas."""

import ephem
import numpy
import pytest

from conftest import obj, minutes, night_jd, run_cascade, vis


class TestSelection:
    def test_tag_required_drops_untagged(self, make_gen):
        p = run_cascade(make_gen([obj()], tag="stromgren"), upto="ObjectMask")

        assert minutes(p) == 0

    def test_tag_required_keeps_matching(self, make_gen):
        p = run_cascade(make_gen([obj(tag="stromgren")], tag="stromgren"), upto="ObjectMask")

        assert minutes(p) == len(p.nightTime)

    def test_no_tag_option_ignores_tags(self, make_gen):
        p = run_cascade(make_gen([obj(tag="stromgren")]), upto="ObjectMask")

        assert minutes(p) == len(p.nightTime)

    @pytest.mark.parametrize("status", ["finished", "inactive"])
    def test_inactive_status_drops_object(self, make_gen, status):
        p = run_cascade(make_gen([obj(status=status)]), upto="ObjectMask")

        assert minutes(p) == 0

    def test_unknown_status_does_not_drop(self, make_gen):
        p = run_cascade(make_gen([obj(status="ongoing")]), upto="ObjectMask")

        assert minutes(p) == len(p.nightTime)

    def test_fwhm_max_drops_when_seeing_worse(self, make_gen):
        p = run_cascade(make_gen([obj(fwhm_max=2.0)], fwhm=3.0), upto="ObjectMask")

        assert minutes(p) == 0

    def test_fwhm_max_keeps_when_seeing_better(self, make_gen):
        p = run_cascade(make_gen([obj(fwhm_max=2.0)], fwhm=1.0), upto="ObjectMask")

        assert minutes(p) == len(p.nightTime)

    def test_fwhm_min_drops_when_seeing_better(self, make_gen):
        p = run_cascade(make_gen([obj(fwhm_min=2.0)], fwhm=1.0), upto="ObjectMask")

        assert minutes(p) == 0

    def test_fwhm_limits_inert_without_option(self, make_gen):
        p = run_cascade(make_gen([obj(fwhm_min=2.0, fwhm_max=3.0)]), upto="ObjectMask")

        assert minutes(p) == len(p.nightTime)

    def test_done_uobi_drops_object(self, make_gen):
        p = run_cascade(make_gen([obj(uobi="abc123")], done_uobi=["abc123"]), upto="ObjectMask")

        assert minutes(p) == 0

    def test_other_done_uobi_keeps_object(self, make_gen):
        p = run_cascade(make_gen([obj(uobi="abc123")], done_uobi=["xyz789"]), upto="ObjectMask")

        assert minutes(p) == len(p.nightTime)

    def test_night_monitor_ignores_done_uobi(self, make_gen):
        p = run_cascade(make_gen([obj(uobi="abc123", tag="night_monitor")], done_uobi=["abc123"], night_monitor=True),
                        upto="ObjectMask")

        assert minutes(p) == len(p.nightTime)


class TestAltitude:
    def test_config_limits_are_applied(self, make_gen):
        p = run_cascade(make_gen([obj()]), upto="MaskVisibility")

        alt = numpy.array(p.ob[0]["calc"]["alt"])
        oczekiwane = (alt >= 35) & (alt <= 80)
        assert numpy.array_equal(vis(p), oczekiwane)

    def test_object_limits_override_config(self, make_gen):
        p = run_cascade(make_gen([obj(h_min=60, h_max=70)]), upto="MaskVisibility")

        alt = numpy.array(p.ob[0]["calc"]["alt"])
        assert numpy.array_equal(vis(p), (alt >= 60) & (alt <= 70))

    def test_tighter_limits_never_give_more_time(self, make_gen):
        luzny = run_cascade(make_gen([obj(h_min=20)]), upto="MaskVisibility")
        ciasny = run_cascade(make_gen([obj(h_min=60)]), upto="MaskVisibility")

        assert minutes(ciasny) < minutes(luzny)


class TestMoon:
    def test_distance_zero_masks_nothing(self, make_gen):
        p = run_cascade(make_gen([obj(min_moon_dist=0)]), upto="MaskMoon")

        assert vis(p, key="min_moon_dist").all()

    def test_distance_180_masks_everything(self, make_gen):
        p = run_cascade(make_gen([obj(min_moon_dist=180)]), upto="MaskMoon")

        assert minutes(p) == 0

    def test_distance_follows_the_moon_minute_by_minute(self, make_gen):
        """Regresja opisu w README: separacja liczy sie per minuta, nie raz na noc."""
        p = run_cascade(make_gen([obj(min_moon_dist=50)]), upto="MaskMoon")

        sep = numpy.array(p.ob[0]["calc"]["moon_sep"])
        assert numpy.array_equal(vis(p, key="min_moon_dist"), sep >= 50)

    def test_max_moon_phase_ignores_altitude(self, make_gen):
        """Maska patrzy wylacznie na oswietlenie tarczy - wysokosc Ksiezyca nie gra roli."""
        p = run_cascade(make_gen([obj(max_moon_phase=101)]), upto="MaskMoon")
        assert vis(p, key="max_moon_phase").all()

        p = run_cascade(make_gen([obj(max_moon_phase=-1)]), upto="MaskMoon")
        assert not vis(p, key="max_moon_phase").any()


class TestWind:
    def test_inert_without_option(self, make_gen):
        p = run_cascade(make_gen([obj()]), upto="MaskWind")

        assert "wind" not in p.ob[0]["visibility"]

    # 0 i 359 sa tu po to, zeby zlapac zawijanie przez polnoc
    @pytest.mark.parametrize("kierunek", [0.0, 90.0, 180.0, 270.0, 359.0])
    def test_forbidden_sector_is_masked(self, make_gen, kierunek):
        bazowy = run_cascade(make_gen([obj()]), upto="MaskVisibility")
        az = numpy.array(bazowy.ob[0]["calc"]["az"])

        p = run_cascade(make_gen([obj()], wind=kierunek), upto="MaskWind")

        odleglosc = numpy.abs(kierunek - az)
        odleglosc = numpy.where(odleglosc > 180, 360 - odleglosc, odleglosc)
        assert numpy.array_equal(vis(p, key="wind"), vis(bazowy) & ~(odleglosc < 95))

    def test_sector_wraps_around_north(self, make_gen):
        """Wiatr z 359 deg i obiekt na 5 deg dzieli 6 deg, nie 354 - czyli musi maskowac.

        Bez zawijania odleglosc wyszlaby 354 > 95 i sektor przy polnocy bylby slepy.
        """
        # dec na polnoc od szerokosci stanowiska - obiekt kulminuje na polnocy,
        # wiec jego azymut przechodzi przez 0
        polnocny = obj(dec="+10:00:00")

        bazowy = run_cascade(make_gen([polnocny]), upto="MaskVisibility")
        az = numpy.array(bazowy.ob[0]["calc"]["az"])
        blisko_polnocy = (az < 5) | (az > 355)
        assert blisko_polnocy.any(), "cel testowy mial przechodzic przez azymut 0"

        p = run_cascade(make_gen([polnocny], wind=359.0), upto="MaskWind")

        # kazda minuta o azymucie ~0 musi byc zabita przez wiatr z 359
        assert not vis(p, key="wind")[blisko_polnocy].any()


class TestTwilight:
    def test_inert_by_default(self, make_gen):
        p = run_cascade(make_gen([obj()]), upto="MaskTwilight")

        assert vis(p, key="after_sunset").all()
        assert vis(p, key="before_sunrise").all()

    def test_sunset_dt_delays_the_start(self, make_gen):
        p = run_cascade(make_gen([obj(sunset_dt=2)]), upto="MaskTwilight")

        nt = numpy.array(p.nightTime)
        assert numpy.array_equal(vis(p, key="after_sunset"), nt >= p.sunset + 2 / 24.0)

    def test_sunrise_dt_advances_the_end(self, make_gen):
        p = run_cascade(make_gen([obj(sunrise_dt=2)]), upto="MaskTwilight")

        nt = numpy.array(p.nightTime)
        assert numpy.array_equal(vis(p, key="before_sunrise"), nt <= p.sunrise - 2 / 24.0)

    def test_buffers_are_measured_from_hsun_crossing(self, make_gen):
        """"Zachod" to przejscie Slonca przez hsun, a nie horyzont geometryczny."""
        p = run_cascade(make_gen([obj(sunset_dt=1)]), upto="MaskTwilight")

        pierwsza = numpy.array(p.nightTime)[vis(p, key="after_sunset")][0]
        assert pierwsza >= p.sunset + 1 / 24.0
        assert pierwsza < p.sunset + 1 / 24.0 + ephem.minute * 2


class TestStartEnd:
    def test_absolute_window(self, make_gen):
        p = run_cascade(
            make_gen([obj(t_start="2026/08/13T02:00:00", t_end="2026/08/13T04:00:00")]),
            upto="MaskStartEnd",
        )

        nt = numpy.array(p.nightTime)
        widoczne = nt[vis(p)]
        assert widoczne.min() >= float(ephem.Date("2026/08/13 02:00:00"))
        assert widoczne.max() <= float(ephem.Date("2026/08/13 04:00:00"))

    def test_ut_only_rolls_over_past_midnight(self, make_gen):
        """Plan startuje 12.08 o 23:39 UT; t_start=02:00 znaczy 13.08, nie 12.08."""
        p = run_cascade(make_gen([obj(t_start="02:00:00")]), upto="MaskStartEnd")

        nt = numpy.array(p.nightTime)
        assert nt[vis(p)].min() >= float(ephem.Date("2026/08/13 02:00:00"))

    def test_window_outside_the_night_kills_the_object(self, make_gen):
        p = run_cascade(
            make_gen([obj(t_start="2026/08/20T02:00:00", t_end="2026/08/20T04:00:00")]),
            upto="MaskStartEnd",
        )

        assert minutes(p) == 0


class TestCycle:
    def test_recent_observation_masks_the_night(self, make_gen, write_lc):
        p = make_gen([obj(name="cel", cycle=5)])
        jd_nocy = None
        p.LoadObjects()
        p.MakeTime()
        jd_nocy = night_jd(p)
        write_lc("cel", "Ic", [jd_nocy[0] - 1.0])  # obserwacja sprzed doby, cycle=5

        for step in ("ObjectMask", "CalcObject", "MaskVisibility", "MaskMoon",
                     "MaskWind", "MaskTwilight", "MaskCycle"):
            getattr(p, step)()

        assert minutes(p) == 0

    def test_old_observation_does_not_mask(self, make_gen, write_lc):
        p = make_gen([obj(name="cel", cycle=5)])
        p.LoadObjects()
        p.MakeTime()
        write_lc("cel", "Ic", [night_jd(p)[0] - 30.0])

        for step in ("ObjectMask", "CalcObject", "MaskVisibility", "MaskMoon",
                     "MaskWind", "MaskTwilight", "MaskCycle"):
            getattr(p, step)()

        assert vis(p, key="cycle").all()

    def test_missing_data_file_leaves_object_available(self, make_gen):
        """Brak archiwum = brak wiedzy, obiekt traktowany jak nigdy nieobserwowany."""
        p = run_cascade(make_gen([obj(name="bezdanych", cycle=5)]), upto="MaskCycle")

        assert "cycle" not in p.ob[0]["visibility"]
        assert minutes(p) > 0

    def test_bad_quality_points_are_ignored(self, make_gen, write_lc):
        p = make_gen([obj(name="cel", cycle=5)])
        p.LoadObjects()
        p.MakeTime()
        # swieza obserwacja, ale odrzucona jakosciowo -> nie powinna maskowac
        write_lc("cel", "Ic", [night_jd(p)[0] - 1.0], quality=[2])

        for step in ("ObjectMask", "CalcObject", "MaskVisibility", "MaskMoon",
                     "MaskWind", "MaskTwilight", "MaskCycle"):
            getattr(p, step)()

        assert minutes(p) > 0

    def test_filter_comes_from_seq_not_from_ph_mk(self, make_gen, write_lc):
        """MaskCycle buduje sciezke z pierwszego filtra seq, ph_mk go nie zmienia."""
        p = make_gen([obj(name="cel", seq="1/V/30", cycle=5, ph_mk="2/Ic/0.05", P=1.0)])
        p.LoadObjects()
        p.MakeTime()
        write_lc("cel", "Ic", [night_jd(p)[0] - 1.0])  # dane w Ic - maja zostac zignorowane

        for step in ("ObjectMask", "CalcObject", "MaskVisibility", "MaskMoon",
                     "MaskWind", "MaskTwilight", "MaskCycle"):
            getattr(p, step)()

        assert "cycle" not in p.ob[0]["visibility"]

        # te same dane w V juz maskuja
        q = make_gen([obj(name="cel", seq="1/V/30", cycle=5)])
        q.LoadObjects()
        q.MakeTime()
        write_lc("cel", "V", [night_jd(q)[0] - 1.0])
        for step in ("ObjectMask", "CalcObject", "MaskVisibility", "MaskMoon",
                     "MaskWind", "MaskTwilight", "MaskCycle"):
            getattr(q, step)()

        assert minutes(q) == 0
