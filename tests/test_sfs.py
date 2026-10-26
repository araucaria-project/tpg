"""Testy skyflat selectora (sfs).

Te same zasady co w conftest.py: config podajemy w calosci (bez config.yaml i lokalnego
override'u), a oczekiwania zalezne od efemeryd wyprowadzamy z `field["calc"]` / `s.sun`.
"""

import math

import ephem
import numpy
import pytest

from tpg.sfs import (SkyflatSelector, NoSkyflatField, load_fields, az_dist, ang_dist, NULL_POINT_ZD,
                     select_filters, flat_age, filter_name)

from conftest import SITE, TEL, TEL_CFG, DATE

SKYFLAT = {
    "n": 7,
    "exp": "a",
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
        cfg["skyflat"] = dict(SKYFLAT)
        cfg[TEL] = dict(TEL_CFG, skyflat=dict(fields_file=str(fields), **(skyflat or {})))
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
        assert fields[0]["priority"] == 0
        assert len(errors) == 3

    def test_priority(self, tmp_path):
        f = tmp_path / "fields.txt"
        f.write_text("SKYFLAT a 01:00:00 -30:00:00 priority=5\n"
                     "SKYFLAT b 01:00:00 -30:00:00 priority=x\n")
        errors = []

        fields = load_fields(str(f), log=errors.append)

        assert [(x["name"], x["priority"]) for x in fields] == [("a", 5)]
        assert len(errors) == 1


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

    @pytest.mark.parametrize("session, now, warning", [
        # rano po polnocy z data kalendarzowa zamiast daty nocy -> poranek nastepnego dnia
        ("morning", "2026/08/12 09:00:00", "starts in"),
        ("morning", "2026/08/13 09:00:00", None),
        ("evening", "2026/08/12 15:00:00", None),
        ("evening", "2026/08/13 15:00:00", "already over"),
        ("evening", "2026/08/10 15:00:00", "starts in"),
    ])
    def test_date_warning(self, make_sfs, session, now, warning):
        s = make_sfs(session=session, now=now)
        s.Initiate()
        s.MakeTime()

        if warning:
            assert f"WARNING: the {session} session" in s.msg and warning in s.msg
        else:
            assert "WARNING" not in s.msg

    def test_morning_follows_evening(self, make_sfs):
        ev = make_sfs(session="evening")
        mo = make_sfs(session="morning")
        for s in (ev, mo):
            s.Initiate()
            s.MakeTime()

        assert ev.end_time < mo.start_time < ev.end_time + 1

    def test_fields_file_required_in_telescope_block(self, make_sfs):
        s = make_sfs()
        del s.cfg[TEL]["skyflat"]["fields_file"]

        with pytest.raises(ValueError, match="fields_file"):
            s.Initiate()

    def test_fields_file_argument_replaces_config(self, make_sfs, tmp_path):
        other = tmp_path / "other.txt"
        other.write_text(FIELDS[0] + "\n")
        s = make_sfs(fields_file=str(other))
        del s.cfg[TEL]["skyflat"]["fields_file"]
        s.Initiate()
        s.LoadFields()

        assert [f["name"] for f in s.fields] == ["BF00b"]

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
        assert r["field_line"].split()[1] == r["name"]

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


GROUPS = [["Hb_n"], ["u_s", "b_s", "v_s"], ["B", "V", "Ic"]]


def with_priority(lines, names, priority):
    return [f"{l} priority={priority}" if l.split()[1] in names else l for l in lines]


class TestPriority:
    def run(self, make_sfs, lines, **kw):
        s = make_sfs(lines, **kw)
        r = s.go()
        return s, r

    def test_equal_priorities_rank_by_distance(self, make_sfs):
        _, plain = self.run(make_sfs, FIELDS)
        _, same = self.run(make_sfs, [f"{l} priority=3" for l in FIELDS])

        assert same["name"] == plain["name"]

    def test_higher_priority_wins_over_distance(self, make_sfs):
        # szersze pasmo wysokosci - przy 40-70 tej nocy tylko jedno pole jest widoczne cala sesje
        wide = {"hmin": 30, "hmax": 85}
        s, r = self.run(make_sfs, FIELDS, skyflat=wide)
        full = [f["name"] for f in s.fields if f["fraction"] == 1.0]
        assert len(full) >= 2, "test needs two fields visible the whole session"
        worse = full[-1]

        s, r = self.run(make_sfs, with_priority(FIELDS, [worse], 10), skyflat=wide)

        assert r["name"] == worse
        assert r["priority"] == 10
        assert r["fallback"] is False

    def test_lower_priority_when_higher_not_visible(self, make_sfs):
        s, r = self.run(make_sfs, FIELDS)
        partial = [f["name"] for f in s.fields if f["fraction"] < 1.0]

        # najwyzszy priorytet tylko dla pol niewidocznych przez cala sesje - nie wygrywaja
        s2, r2 = self.run(make_sfs, with_priority(FIELDS, partial, 10))

        assert r2["name"] == r["name"]
        assert r2["priority"] == 0
        assert r2["fallback"] is False

    def test_fallback_prefers_longer_visibility_over_priority(self, make_sfs):
        s, r = self.run(make_sfs, FIELDS, skyflat={"hmin": 30, "hmax": 50})
        assert r["fallback"] is True
        shorter = [f["name"] for f in s.fields if 0 < f["fraction"] < r["fraction"]]
        assert shorter

        s2, r2 = self.run(make_sfs, with_priority(FIELDS, shorter, 10), skyflat={"hmin": 30, "hmax": 50})

        assert r2["name"] == r["name"]


class TestFilters:
    def test_flat_age(self):
        assert flat_age(1.5) == 1.5
        assert flat_age(3) == 3.0
        assert math.isinf(flat_age(None))
        assert flat_age(-0.2) == 0.0

    def test_oldest_from_every_group_evening_order(self):
        history = {"Hb_n": 6, "u_s": 25, "b_s": 55, "v_s": 1, "B": 1, "V": 16, "Ic": 2}

        chosen = select_filters(GROUPS, history, "evening")

        assert [c["filter"] for c in chosen] == ["Hb_n", "b_s", "V"]
        assert [c["age"] for c in chosen] == [6, 55, 16]

    def test_morning_is_reversed(self):
        evening = select_filters(GROUPS, {}, "evening")
        morning = select_filters(GROUPS, {}, "morning")

        assert morning == evening[::-1]

    def test_never_taken_wins_and_ties_keep_group_order(self):
        history = {"u_s": 300, "b_s": None, "B": 6, "V": 6, "Ic": 6}

        chosen = select_filters(GROUPS, history, "evening")

        # Hb_n i v_s nie ma w historii - tez "never"; b_s jest pierwszym "never" w swojej grupie
        assert [c["filter"] for c in chosen] == ["Hb_n", "b_s", "B"]
        assert math.isinf(chosen[1]["age"])

    def test_selector_uses_config_groups(self, make_sfs):
        s = make_sfs(session="morning", skyflat={"groups": [["u"], "B, V"]}, flat_history={"B": 80})
        r = s.go()

        assert [f["filter"] for f in r["filters"]] == ["V", "u"]
        assert s.groups == [["u"], ["B", "V"]]
        assert math.isinf(r["filters"][0]["age"])

    def test_no_groups_no_filters(self, make_sfs):
        r = make_sfs().go()

        assert r["filters"] == []


class TestPlanLine:
    def test_skyflat_line_with_seq(self, make_sfs):
        r = make_sfs(session="morning", skyflat={"groups": [["u"], ["B", "V"]]}).go()

        assert r["seq"] == "7/B/a,7/u/a"
        assert r["line"] == f'SKYFLAT {r["name"]} {r["ra"]} {r["dec"]} seq=7/B/a,7/u/a'

    def test_n_and_exp_from_config(self, make_sfs):
        r = make_sfs(skyflat={"groups": [["u"], ["V"]], "n": 5, "exp": 10}).go()

        assert r["seq"] == "5/u/10,5/V/10"

    def test_full_sequence_in_group_used_as_given(self, make_sfs):
        r = make_sfs(skyflat={"groups": [["u"], ["5/Ic/3", "V"]]}, flat_history={"Ic": 30, "V": 2}).go()

        assert r["seq"] == "7/u/a,5/Ic/3"
        assert r["filters"][1]["filter"] == "Ic"
        assert r["filters"][1]["entry"] == "5/Ic/3"

    def test_full_sequence_history_by_filter_name(self, make_sfs):
        # historia jest po nazwie filtra - Ic ma swiezsze flaty, wiec wygrywa V z domyslna sekwencja
        r = make_sfs(skyflat={"groups": ["5/Ic/3, V"]}, flat_history={"Ic": 1, "V": 20}).go()

        assert r["seq"] == "7/V/a"

    @pytest.mark.parametrize("entry", ["5/Ic", "x/Ic/3", "5//3", "5/Ic/3/1", ""])
    def test_bad_entry_raises(self, make_sfs, entry):
        s = make_sfs(skyflat={"groups": [["u", entry]]})

        with pytest.raises(ValueError, match="skyflat group entry"):
            s.Initiate()

    def test_filter_name(self):
        assert filter_name("Ic") == "Ic"
        assert filter_name("5/Ic/3") == "Ic"
        assert filter_name("7/u_s/a") == "u_s"

    def test_no_groups_line_without_seq(self, make_sfs):
        r = make_sfs().go()

        assert r["line"] == f'SKYFLAT {r["name"]} {r["ra"]} {r["dec"]}'


class TestCli:
    def test_session_without_field_does_not_stop_the_other(self, monkeypatch, capsys):
        import tpg.sfs as sfs

        class FakeSelector:
            def __init__(self, tel, date, session, **kw):
                self.session = session

            def go(self):
                if self.session == "evening":
                    raise NoSkyflatField("no field in the evening")
                return {"line": "SKYFLAT F 01:00:00 -30:00:00 seq=7/V/a"}

        monkeypatch.setattr(sfs, "SkyflatSelector", FakeSelector)
        monkeypatch.setattr("sys.argv", ["sfs", "test", "2026/08/12"])

        with pytest.raises(SystemExit) as e:
            sfs.main()

        out, err = capsys.readouterr()
        assert e.value.code == 1
        assert out == "SKYFLAT F 01:00:00 -30:00:00 seq=7/V/a\n"
        assert "no field in the evening" in err

    def test_all_sessions_ok_exit_0(self, monkeypatch, capsys):
        import tpg.sfs as sfs

        class FakeSelector:
            def __init__(self, tel, date, session, **kw):
                self.session = session

            def go(self):
                return {"line": f"SKYFLAT {self.session}"}

        monkeypatch.setattr(sfs, "SkyflatSelector", FakeSelector)
        monkeypatch.setattr("sys.argv", ["sfs", "test", "2026/08/12"])

        sfs.main()

        assert capsys.readouterr().out == "SKYFLAT evening\nSKYFLAT morning\n"
