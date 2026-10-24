"""Alokator i budowa planu wyjsciowego."""

import numpy
import pytest

from conftest import TEL, obj, run_cascade
from tpg.telescope_plan_generator import NoObservableObjects


def alokuj(p, seed=1):
    """Domknij pipeline do konca alokacji."""
    run_cascade(p)
    p.seed = seed
    p.Waga()
    p.RandomizeList()
    p.allocate()
    return p


def kolejnosc(p):
    """Nazwy obiektow w kolejnosci wstawienia do nocy."""
    return [
        p.ob[int(p.nightObject[n])]["name"]
        for n, zmiana in enumerate(p.nightObChange)
        if zmiana
    ]


class TestWeights:
    def test_never_visible_object_gets_zero_weight(self, make_gen):
        p = run_cascade(make_gen([obj(name="widoczny"), obj(name="martwy", h_min=89)]))
        p.Waga()

        wagi = {o["name"]: o["waga"] for o in p.ob}
        assert wagi["martwy"] == 0
        assert wagi["widoczny"] > 0

    def test_zero_weight_object_is_not_drawn(self, make_gen):
        p = run_cascade(make_gen([obj(name="widoczny"), obj(name="martwy", h_min=89)]))
        p.seed = 1
        p.Waga()
        p.RandomizeList()

        nazwy = [p.ob[i]["name"] for i in p.random_index_list]
        assert "martwy" not in nazwy

    def test_scarce_object_gets_higher_weight(self, make_gen):
        """Waga premiuje obiekty o waskim oknie widocznosci."""
        p = run_cascade(make_gen([obj(name="szeroki"), obj(name="waski", h_min=70)]))
        p.Waga()

        wagi = {o["name"]: o["waga"] for o in p.ob}
        assert wagi["waski"] > wagi["szeroki"]

    def test_no_observable_objects_raises(self, make_gen):
        p = run_cascade(make_gen([obj(h_min=89)]))

        with pytest.raises(ValueError, match="No observable objects"):
            p.Waga()

    def test_missing_tag_is_reported_as_such(self, make_gen):
        p = run_cascade(make_gen([obj(name="bez_tagu")], tag="night_monitor"))

        with pytest.raises(NoObservableObjects, match="no object in master file has tag 'night_monitor'"):
            p.Waga()

    def test_tagged_but_invisible_is_not_reported_as_missing_tag(self, make_gen):
        p = run_cascade(make_gen([obj(tag="night_monitor", h_min=89)], tag="night_monitor"))

        with pytest.raises(NoObservableObjects, match="none of 1 objects is observable"):
            p.Waga()

    def test_same_seed_gives_same_order(self, make_gen):
        nazwy = [obj(name=f"c{i}") for i in range(8)]
        a = run_cascade(make_gen(nazwy))
        a.seed = 4242
        a.Waga()
        a.RandomizeList()
        b = run_cascade(make_gen(nazwy))
        b.seed = 4242
        b.Waga()
        b.RandomizeList()

        assert list(a.random_index_list) == list(b.random_index_list)


class TestAllocation:
    def test_higher_priority_goes_first(self, make_gen):
        p = alokuj(make_gen([obj(name="niski", priority=1), obj(name="wysoki", priority=9)]))

        assert kolejnosc(p)[0] == "wysoki"

    def test_negative_priority_is_never_allocated(self, make_gen):
        p = alokuj(make_gen([obj(name="zwykly", priority=5), obj(name="wylaczony", priority=-1)]))

        assert "wylaczony" not in kolejnosc(p)
        assert "zwykly" in kolejnosc(p)

    def test_object_longer_than_the_night_is_skipped(self, make_gen):
        p = alokuj(make_gen([obj(name="olbrzym", ob_time=5000), obj(name="normalny")]))

        assert kolejnosc(p) == ["normalny"]

    def test_allocated_block_is_contiguous_and_marked_busy(self, make_gen):
        p = alokuj(make_gen([obj(name="cel", ob_time=30)]))

        zajete = numpy.array([o == 0 for o in p.nightObject])
        indeksy = numpy.flatnonzero(zajete)
        assert len(indeksy) == 31  # slot rezerwuje ob_time + 1 minute
        assert numpy.array_equal(indeksy, numpy.arange(indeksy[0], indeksy[0] + 31))
        assert not any(p.nightFree[indeksy[0]:indeksy[-1] + 1])

    def test_object_without_cycle_is_scheduled_once(self, make_gen):
        p = alokuj(make_gen([obj(name="jednorazowy", ob_time=10)]))

        assert kolejnosc(p).count("jednorazowy") == 1

    def test_object_with_cycle_returns_later(self, make_gen):
        """cycle=0.1 to 2.4 h - w ciagu ~10 h nocy obiekt musi wrocic kilka razy."""
        p = alokuj(make_gen([obj(name="powracajacy", ob_time=10, cycle=0.1)]))

        wystapienia = kolejnosc(p).count("powracajacy")
        assert wystapienia > 1

    def test_cycle_gap_is_measured_from_slot_start(self, make_gen):
        p = alokuj(make_gen([obj(name="powracajacy", ob_time=10, cycle=0.1)]))

        starty = [n for n, z in enumerate(p.nightObChange) if z]
        odstepy = [b - a for a, b in zip(starty, starty[1:])]
        assert all(d >= int(0.1 * 1440) for d in odstepy)

    def test_uobi_group_is_retired_together(self, make_gen):
        wspolny = [
            obj(name="a", uobi="grupa1", ob_time=10),
            obj(name="b", uobi="grupa1", ob_time=10),
        ]
        p = alokuj(make_gen(wspolny))

        assert len(kolejnosc(p)) == 1

    def test_without_uobi_both_are_scheduled(self, make_gen):
        """Kontrola do testu wyzej - bez uobi noc miesci oba obiekty."""
        p = alokuj(make_gen([obj(name="a", ob_time=10), obj(name="b", ob_time=10)]))

        assert sorted(kolejnosc(p)) == ["a", "b"]


class TestExport:
    def _plan(self, p):
        p.export()
        return p.plan

    def test_starts_with_wait_sunset_when_no_ut(self, make_gen):
        p = alokuj(make_gen([obj(name="cel", ob_time=10)]))

        assert self._plan(p)[0] == "WAIT sunset=-18"

    def test_starts_with_wait_ut_when_ut_given(self, make_gen):
        p = alokuj(make_gen([obj(name="cel", ob_time=10)], ut="03:00:00"))

        assert self._plan(p)[0].startswith("WAIT ut=")

    def test_object_line_carries_observing_attributes(self, make_gen):
        p = alokuj(make_gen([obj(name="cel", ob_time=10, sciprog="test", pi="ocm")]))
        linia = [l for l in self._plan(p) if l.startswith("OBJECT")][0]

        assert "cel" in linia
        assert "seq=1/Ic/30" in linia
        assert "sciprog=test" in linia
        assert "pi=ocm" in linia

    def test_scheduling_attributes_are_stripped(self, make_gen):
        p = alokuj(
            make_gen([obj(name="cel", ob_time=10, cycle=3, P=1.5, tag="x", status="ongoing")])
        )
        linia = [l for l in self._plan(p) if l.startswith("OBJECT")][0]

        for zakazane in ("priority=", "cycle=", "P=", "tag=", "status=", "ob_time="):
            assert zakazane not in linia

    def test_gap_before_an_object_produces_a_wait(self, make_gen):
        """Obiekt z t_start w srodku nocy musi byc poprzedzony WAIT ut=."""
        p = alokuj(make_gen([obj(name="pozny", ob_time=10, t_start="05:00:00")]))
        plan = self._plan(p)

        i = next(n for n, l in enumerate(plan) if l.startswith("OBJECT"))
        assert plan[i - 1].startswith("WAIT ut=")

    def test_makro_replaces_the_object_line(self, make_gen, tmp_path):
        makro = tmp_path / "makro.txt"
        makro.write_text("OBJECT kalibracja 00:00:00 +00:00:00 seq=1/V/1\nSTOP\n")

        p = alokuj(make_gen([obj(name="cel", ob_time=10, makro=str(makro))]))
        plan = self._plan(p)

        assert "OBJECT kalibracja 00:00:00 +00:00:00 seq=1/V/1" in plan
        assert "STOP" in plan
        assert not any(l.startswith("OBJECT cel") for l in plan)

    def test_missing_makro_file_does_not_crash(self, make_gen):
        p = alokuj(make_gen([obj(name="cel", ob_time=10, makro="/nonexistent/makro.txt")]))

        plan = self._plan(p)
        assert not any("nonexistent" in l for l in plan)

    def test_start_and_end_makro_wrap_the_plan(self, make_gen, tmp_path):
        start = tmp_path / "start.txt"
        start.write_text("# POCZATEK\n")
        koniec = tmp_path / "koniec.txt"
        koniec.write_text("# KONIEC\n")

        p = alokuj(
            make_gen(
                [obj(name="cel", ob_time=10)],
                tel_cfg={"start_makro": str(start), "end_makro": str(koniec)},
            )
        )
        plan = self._plan(p)

        assert plan[0] == "# POCZATEK"
        assert plan[-1] == "# KONIEC"

    def test_makros_can_be_disabled(self, make_gen, tmp_path):
        start = tmp_path / "start.txt"
        start.write_text("# POCZATEK\n")

        p = alokuj(
            make_gen(
                [obj(name="cel", ob_time=10)],
                tel_cfg={"start_makro": str(start)},
                add_start_makro=False,
            )
        )

        assert "# POCZATEK" not in self._plan(p)


class TestGoldenPlan:
    """Pelny przebieg na ustalonym ziarnie.

    Ten test celowo przybija konkretny wynik. Jesli padnie po zmianie kodu albo po
    aktualizacji ephem/pyaraucaria, trzeba obejrzec roznice i swiadomie zdecydowac,
    czy nowy plan jest poprawny - a nie odruchowo nadpisac oczekiwanie.
    """

    MASTER = [
        obj(name="alfa", priority=9, ob_time=60),
        obj(name="beta", priority=5, ob_time=45),
        obj(name="gamma", priority=7, ob_time=30, t_start="05:00:00"),
    ]

    def test_plan_is_reproducible(self, make_gen):
        a = alokuj(make_gen(self.MASTER), seed=7)
        a.export()
        b = alokuj(make_gen(self.MASTER), seed=7)
        b.export()

        assert a.plan == b.plan

    def test_plan_matches_expected_shape(self, make_gen):
        p = alokuj(make_gen(self.MASTER), seed=7)
        p.export()

        obiekty = [l.split()[1] for l in p.plan if l.startswith("OBJECT")]
        assert obiekty[0] == "alfa", "najwyzszy priorytet powinien otwierac noc"
        assert set(obiekty) == {"alfa", "beta", "gamma"}
        assert p.plan[0] == "WAIT sunset=-18"

    def test_gamma_respects_its_time_window(self, make_gen):
        import ephem

        p = alokuj(make_gen(self.MASTER), seed=7)

        start = next(
            n for n, z in enumerate(p.nightObChange)
            if z and p.ob[int(p.nightObject[n])]["name"] == "gamma"
        )
        assert p.nightTime[start] >= ephem.Date("2026/08/13 05:00:00")


class TestNightMonitor:
    def nm(self, make_gen, lines, gap=30, seed=1):
        p = make_gen(lines, night_monitor=True, tel_cfg={"night_monitor_gap": gap})
        return alokuj(p, seed=seed)

    def starty(self, p):
        return [n for n, zmiana in enumerate(p.nightObChange) if zmiana]

    def test_only_night_monitor_tag_is_allocated(self, make_gen):
        p = self.nm(make_gen, [obj(name="mon", tag="night_monitor"), obj(name="zwykly")])

        assert set(kolejnosc(p)) == {"mon"}

    def test_object_is_repeated_with_gap(self, make_gen):
        gap = 30
        p = self.nm(make_gen, [obj(name="mon", tag="night_monitor")], gap=gap)

        starty = self.starty(p)
        assert len(starty) > 3
        slots = int(float(p.ob[0]["slotTime"]))
        vis_all = p.ob[0]["visibility"]["all"]
        for a, b in zip(starty, starty[1:]):
            # obserwacja zajmuje slots+1 minut, potem co najmniej gap wolnych minut;
            # dluzej tylko wtedy, gdy obiekt akurat nie byl widoczny (np. hmax przy gorowaniu)
            koniec_przerwy = a + slots + 1 + gap
            assert b >= koniec_przerwy
            if all(vis_all[koniec_przerwy:koniec_przerwy + slots + 1]):
                assert b == koniec_przerwy
            assert all(p.nightFree[a + slots + 1:b])

    def test_gap_becomes_wait_ut_in_plan(self, make_gen):
        p = self.nm(make_gen, [obj(name="mon", tag="night_monitor")])
        p.export()

        obiekty = [i for i, l in enumerate(p.plan) if l.startswith("OBJECT")]
        for i in obiekty[1:]:
            assert p.plan[i - 1].startswith("WAIT ut=")

    def test_higher_priority_wins_when_available(self, make_gen):
        p = self.nm(make_gen, [
            obj(name="niski", tag="night_monitor", priority=1),
            obj(name="wysoki", tag="night_monitor", priority=9),
        ])

        assert set(kolejnosc(p)) == {"wysoki"}

    def test_lower_priority_fills_when_higher_not_visible(self, make_gen):
        p = self.nm(make_gen, [
            obj(name="niski", tag="night_monitor", priority=1),
            obj(name="wysoki", tag="night_monitor", priority=9, h_min=70),
        ])

        assert set(kolejnosc(p)) == {"niski", "wysoki"}

    def test_equal_priority_is_redrawn_each_step(self, make_gen):
        p = self.nm(make_gen, [obj(name=f"m{i}", tag="night_monitor") for i in range(3)])

        assert len(set(kolejnosc(p))) == 3

    def test_cycle_still_blocks_repeat(self, make_gen):
        p = self.nm(make_gen, [obj(name="mon", tag="night_monitor", cycle=1)])

        assert kolejnosc(p) == ["mon"]
