"""Fazowanie: okno obserwacyjne (ph_start/ph_end) i maskowanie pokrycia (ph_mk).

Trzy testy oznaczone REGRESJA pilnuja bledow znalezionych 2026-08-12. Kazdy z nich
padal cicho - plan wygladal poprawnie, a byl liczony ze zlych faz. Jesli ktorykolwiek
zacznie padac, plan znowu klamie, nie wywalajac sie.
"""

import numpy
import pytest

from conftest import minutes, night_jd, obj, run_cascade, vis

# Okres dobrany celowo, nie "na oko".
#
# ephem.Date liczy dni od 1899/12/31 12:00, JD od poludnia 4713 p.n.e. - roznica to
# 2415020 dni. W fazie ta roznica wynosi (2415020 / P) mod 1. Dla okresow okraglych
# - 0.1, 0.5, 1.0, 0.16 - wychodzi dokladnie 0 i pomylenie skal czasu staje sie
# CALKOWICIE NIEWIDOCZNE. Test regresyjny na P=0.1 przechodzilby zawsze, takze na
# zepsutym kodzie. 0.27 daje przesuniecie 0.52, czyli niemal maksymalna czulosc.
P = 0.27
HJD0 = 2460000.0

# przesuniecie fazy, jakie wprowadzilaby zla skala czasu
PRZESUNIECIE_ZLEJ_SKALI = (2415020.0 / P) % 1


def test_okres_testowy_jest_czuly_na_skale_czasu():
    """Straznik testow ponizej: gdyby okres byl 'okragly', nie wykrylyby niczego."""
    czulosc = min(PRZESUNIECIE_ZLEJ_SKALI, 1 - PRZESUNIECIE_ZLEJ_SKALI)
    assert czulosc > 0.4, f"P={P} nie nadaje sie na test skali czasu (czulosc {czulosc:.3f})"


def fazy_nocy(p, period=P, hjd0=HJD0):
    return (night_jd(p) - hjd0) / period % 1


def zle_fazy_nocy(p, period=P, hjd0=HJD0):
    """Fazy tak, jak liczylby je kod uzywajacy surowego ephem.Date zamiast JD."""
    return (numpy.array(p.nightTime) - hjd0) / period % 1


class TestPhaseWindow:
    def test_window_narrows_to_the_requested_phases(self, make_gen):
        p = run_cascade(
            make_gen([obj(P=P, hjd0=HJD0, ph_start=0.1, ph_end=0.3)]),
            upto="MaskPhaseStartEnd",
        )

        faza = fazy_nocy(p)
        w_oknie = (faza >= 0.1) & (faza <= 0.3)
        # kazda widoczna minuta lezy w oknie
        assert not (vis(p) & ~w_oknie).any()
        # i cos w ogole zostalo
        assert minutes(p) > 0

    def test_window_is_cyclic_when_end_precedes_start(self, make_gen):
        p = run_cascade(
            make_gen([obj(P=P, hjd0=HJD0, ph_start=0.8, ph_end=0.1)]),
            upto="MaskPhaseStartEnd",
        )

        faza = fazy_nocy(p)
        w_oknie = (faza >= 0.8) | (faza <= 0.1)
        assert not (vis(p) & ~w_oknie).any()
        assert minutes(p) > 0

    def test_window_uses_julian_date_not_ephem_date(self, make_gen):
        """Gdyby faze liczono z surowego ephem.Date, wyszlaby przesunieta o 2415020 d."""
        p = make_gen([obj(P=P, hjd0=HJD0, ph_start=0.1, ph_end=0.3)])
        p.LoadObjects()
        p.MakeTime()
        for step in ("ObjectMask", "CalcObject", "MaskVisibility", "MaskMoon",
                     "MaskWind", "MaskTwilight", "MaskCycle", "MaskStartEnd"):
            getattr(p, step)()
        przed = vis(p).copy()
        p.MaskPhaseStartEnd()

        dobra = fazy_nocy(p)
        zla = zle_fazy_nocy(p)
        assert numpy.array_equal(vis(p), przed & (dobra >= 0.1) & (dobra <= 0.3))
        assert not numpy.array_equal(vis(p), przed & (zla >= 0.1) & (zla <= 0.3))

    def test_regresja_window_applies_without_ph_mk(self, make_gen):
        """REGRESJA: MaskPhaseStartEnd wczytywalo krzywa blasku, ktorej nie uzywalo,
        i przy jej braku porzucalo cala maske. Okno fazowe zalezy wylacznie od P i hjd0."""
        p = run_cascade(
            make_gen([obj(P=P, hjd0=HJD0, ph_start=0.1, ph_end=0.3)]),
            upto="MaskPhaseStartEnd",
        )

        assert "ph_start" in p.ob[0]["visibility"], "maska fazowa w ogole sie nie uruchomila"
        assert 0 < minutes(p) < len(p.nightTime)

    def test_window_needs_no_light_curve_on_disk(self, make_gen, tmp_path):
        p = run_cascade(
            make_gen([obj(P=P, hjd0=HJD0, ph_start=0.1, ph_end=0.3)]),
            upto="MaskPhaseStartEnd",
        )

        assert not (tmp_path / "data").exists()
        assert minutes(p) > 0

    def test_incomplete_specification_is_ignored(self, make_gen):
        """Bez ph_end okno nie ma sensu - maska musi odpuscic, a nie zgadywac."""
        p = run_cascade(make_gen([obj(P=P, hjd0=HJD0, ph_start=0.1)]), upto="MaskPhaseStartEnd")

        assert "ph_start" not in p.ob[0]["visibility"]


class TestPhaseCoverage:
    def _przygotuj(self, make_gen, write_lc, fazy_danych, ph_mk, **kw):
        """Zbuduj obiekt z krzywa blasku o zadanych fazach i przepusc kaskade."""
        p = make_gen([obj(name="cel", P=P, hjd0=HJD0, ph_mk=ph_mk, **kw)])
        p.LoadObjects()
        p.MakeTime()
        jd = [HJD0 + f * P + k * P for k, f in enumerate(fazy_danych)]
        write_lc("cel", "Ic", jd)
        for step in ("ObjectMask", "CalcObject", "MaskVisibility", "MaskMoon",
                     "MaskWind", "MaskTwilight", "MaskCycle", "MaskStartEnd",
                     "MaskPhaseStartEnd"):
            getattr(p, step)()
        przed = vis(p).copy()
        p.MaskPhase()
        return p, przed

    def test_covered_bin_is_masked(self, make_gen, write_lc):
        p, przed = self._przygotuj(make_gen, write_lc, [0.50, 0.50, 0.50], "2/Ic/0.05")

        zabite = przed & ~vis(p)
        assert zabite.any()
        faza = fazy_nocy(p)
        # zabite minuty leza w oknie wokol 0.5
        assert numpy.abs(faza[zabite] - 0.5).max() <= 0.05 + 1e-9

    def test_below_min_count_does_not_mask(self, make_gen, write_lc):
        """Jeden punkt przy progu 2 nie wystarcza."""
        p, przed = self._przygotuj(make_gen, write_lc, [0.50], "2/Ic/0.05")

        assert numpy.array_equal(przed, vis(p))

    def test_uncovered_phases_survive(self, make_gen, write_lc):
        p, przed = self._przygotuj(make_gen, write_lc, [0.50, 0.50], "2/Ic/0.05")

        faza = fazy_nocy(p)
        daleko = numpy.abs(faza - 0.5) > 0.2
        assert (przed & daleko & vis(p)).sum() == (przed & daleko).sum()

    def test_regresja_phase_computed_from_julian_date(self, make_gen, write_lc):
        """REGRESJA: night_ph liczono z ephem.Date (Dublin JD) zamiast z JD.

        Roznica to 2415020 dni, wiec maska trafiala w losowe fazy. Tu dane siedza
        w fazie 0.5; przy zlej skali czasu zabite minuty nie mialyby z nia zwiazku.
        """
        p, przed = self._przygotuj(make_gen, write_lc, [0.50, 0.50, 0.50], "3/Ic/0.05")

        zabite = przed & ~vis(p)
        assert zabite.any(), "nic nie zamaskowano - maska nie trafila w dane"

        faza = fazy_nocy(p)
        assert numpy.abs(faza[zabite] - 0.5).max() <= 0.05 + 1e-9

        # kontrola: przy blednej skali czasu zabity zbior bylby inny
        zla = zle_fazy_nocy(p)
        zle_zabite = przed & (numpy.abs(zla - 0.5) <= 0.05)
        assert not numpy.array_equal(zabite, zle_zabite)

    def test_regresja_window_wraps_across_the_seam(self, make_gen, write_lc):
        """REGRESJA: okno [f-w, f+w] bylo obcinane przez zakres [0,1).

        Dane w fazie 0.99 nie maskowaly minut w fazie 0.01, mimo ze dzieli je 0.02.
        Obie galezie - przy 1 i przy 0 - musza zostac zabite.
        """
        p, przed = self._przygotuj(make_gen, write_lc, [0.99, 0.99, 0.99], "1/Ic/0.05")

        faza = fazy_nocy(p)
        zabite = przed & ~vis(p)

        przy_jedynce = zabite & (faza >= 0.94)
        przy_zerze = zabite & (faza <= 0.04)
        assert przy_jedynce.any(), "galaz przy fazie 1 niezamaskowana"
        assert przy_zerze.any(), "galaz przy fazie 0 niezamaskowana - okno sie nie zawija"

        # i nic poza oknem cyklicznym
        odleglosc = numpy.minimum(numpy.abs(faza - 0.99), 1.0 - numpy.abs(faza - 0.99))
        assert odleglosc[zabite].max() <= 0.05 + 1e-9

    def test_missing_light_curve_leaves_object_alone(self, make_gen):
        p = run_cascade(make_gen([obj(P=P, hjd0=HJD0, ph_mk="2/Ic/0.05")]))

        assert "ph_mk" not in p.ob[0]["visibility"]
        assert minutes(p) > 0

    @pytest.mark.parametrize("zly", ["bzdura", "2//0.05", "Ic/2/0.05"])
    def test_malformed_ph_mk_is_rejected_by_the_validator(self, make_gen, zly):
        """Bariera stoi przy wczytywaniu, nie w masce.

        Schemat wymusza wzorzec no_of_obs/band/phase_bin, wiec zla specyfikacja
        wywala cala linie master file. try/except w MaskPhase jest dla takich danych
        nieosiagalny - zabezpiecza tylko przed obiektami zbudowanymi z pominieciem
        walidatora.
        """
        p = make_gen([obj(name="cel", P=P, hjd0=HJD0, ph_mk=zly)])
        p.LoadObjects()

        assert p.ob == []
