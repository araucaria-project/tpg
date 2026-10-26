"""tpg --skyflats: SKYFLAT z sfs wokol nocnej czesci planu.

Konfiguracja skyflat jak w test_sfs.py (pola i limity w tmp_path, nic z maszyny).
Oczekiwania dotyczace sesji wyprowadzamy z samego sfs, a nie z wpisanych godzin.
"""

import pytest

from conftest import TEL, obj
from test_allocate_export import alokuj
from test_sfs import FIELDS, SKYFLAT, NEVER_UP


@pytest.fixture
def make_flat_gen(make_gen, tmp_path):
    """Generator z blokiem skyflat w configu; pola z listy `fields`."""

    def _make(fields=FIELDS, groups=(("V", "Ic"),), tel_cfg=None, **kwargs):
        f = tmp_path / "fields.txt"
        f.write_text("\n".join(fields) + "\n")
        p = make_gen([obj(name="cel", ob_time=10)], tel_cfg=tel_cfg, skyflats=True, **kwargs)
        p.cfg["skyflat"] = dict(SKYFLAT)
        p.cfg[TEL]["skyflat"] = {"fields_file": str(f), "groups": [list(g) for g in groups]}
        return p

    return _make


def skyflats(plan):
    return [(i, l) for i, l in enumerate(plan) if l.startswith("SKYFLAT")]


class TestSkyflats:
    def test_evening_before_night_part_morning_after(self, make_flat_gen):
        p = alokuj(make_flat_gen())
        p.export()

        assert p.plan[0].startswith("SKYFLAT")
        assert p.plan[1] == "WAIT sunset=-18"
        assert p.plan[-1].startswith("SKYFLAT")
        assert len(skyflats(p.plan)) == 2
        assert p.skyflat_info["evening"]["line"] == p.plan[0]
        assert p.skyflat_info["morning"]["line"] == p.plan[-1]

    def test_without_option_no_skyflats(self, make_gen):
        p = alokuj(make_gen([obj(name="cel", ob_time=10)]))
        p.export()

        assert skyflats(p.plan) == []
        assert p.skyflat_info == {}

    def test_makros_wrap_the_skyflats(self, make_flat_gen, tmp_path):
        start = tmp_path / "start.txt"
        start.write_text("DOMEOPEN\n")
        koniec = tmp_path / "koniec.txt"
        koniec.write_text("DOMECLOSE\n")

        p = alokuj(make_flat_gen(tel_cfg={"start_makro": str(start), "end_makro": str(koniec)}))
        p.export()

        assert p.plan[0] == "DOMEOPEN"
        assert p.plan[1].startswith("SKYFLAT")
        assert p.plan[-2].startswith("SKYFLAT")
        assert p.plan[-1] == "DOMECLOSE"

    def test_plan_started_at_night_gets_only_morning(self, make_flat_gen):
        # 03:00 UT 12.08 to noc 11/12.08 na OCM - wieczorna sesja juz minela
        p = alokuj(make_flat_gen(ut="03:00:00"))
        p.export()

        assert len(skyflats(p.plan)) == 1
        assert p.plan[-1].startswith("SKYFLAT")
        assert p.skyflat_info["evening"]["line"] is None
        assert "skipped" in p.skyflat_info["evening"]["note"]

    def test_now_after_evening_session_gets_only_morning(self, make_flat_gen):
        # plan "od zachodu", ale generowany juz po wieczornych flatach (TOI podaje now)
        p = alokuj(make_flat_gen(now="2026/08/13 01:00:00"))
        p.export()

        assert p.plan[0] == "WAIT sunset=-18"
        assert len(skyflats(p.plan)) == 1
        assert p.plan[-1].startswith("SKYFLAT")

    def test_now_before_the_night_keeps_both(self, make_flat_gen):
        p = alokuj(make_flat_gen(now="2026/08/12 15:00:00"))
        p.export()

        assert len(skyflats(p.plan)) == 2

    def test_history_per_session(self, make_flat_gen):
        history = {"evening": {"V": 50, "Ic": 1}, "morning": {"V": 1, "Ic": 50}}
        p = alokuj(make_flat_gen(flat_history=history))
        p.export()

        assert "seq=7/V/a" in p.plan[0]
        assert "seq=7/Ic/a" in p.plan[-1]

    def test_no_field_skips_session_but_keeps_plan(self, make_flat_gen):
        p = alokuj(make_flat_gen(fields=[NEVER_UP]))
        p.export()

        assert skyflats(p.plan) == []
        assert p.plan[0] == "WAIT sunset=-18"
        assert any(l.startswith("OBJECT cel") for l in p.plan)
        for session in ("evening", "morning"):
            assert p.skyflat_info[session]["line"] is None
            assert "WARNING" in p.skyflat_info[session]["note"]

    def test_same_wind_as_the_plan(self, make_flat_gen, monkeypatch):
        import tpg.sfs as sfs
        seen = []
        orig = sfs.SkyflatSelector.__init__

        def spy(self, *a, **kw):
            seen.append(kw.get("wind"))
            orig(self, *a, **kw)

        monkeypatch.setattr(sfs.SkyflatSelector, "__init__", spy)
        p = alokuj(make_flat_gen(wind=123.0))
        p.export()

        assert seen == [123.0, 123.0]
