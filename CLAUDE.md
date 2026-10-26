# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`tpg` (telescope plan generator) builds nightly observing schedules for four telescopes
(`wk06`, `zb08`, `jk15`, `iris`, plus a `test` config) at the Las Campanas / Cerro Armazones
site. Given a night's date, it reads a "master file" of candidate objects, computes each
object's visibility for that night (altitude, moon, wind, phase, timing constraints),
allocates objects into time slots by priority with weighted randomization, and writes an
observing plan in TOI/TEXTUI-readable format.

## Commands

Install / update dependencies (Poetry):
```
poetry install
poetry update
```

Run the generator:
```
poetry run tpg {wk06,zb08,jk15,iris,test} yyyy/mm/dd [hh:mm:ss] [--tag TAG] [--wind DEG] [--fwhm FWHM] \
    [--seed N] [--loud] [--no_save] [--no_ctc] [--no_start_makro] [--no_end_makro] [--master_file PATH] \
    [--night_monitor]
```
Example (see README.md for the full argument table and master-file syntax):
```
poetry run tpg zb08 2023/08/1 >| plans/zb08_230810.txt
poetry run tpg test 2023/08/1 --loud --no_save
```
Use the `test` telescope with `--master_file`/`--no_save` for quick manual iteration against
`test/test_master.txt` without touching the real `plans/` output directory or config-defined paths.

Other entry points defined in `pyproject.toml`:
```
poetry run chexp <factor> [-f FILTER]     # rescale n/filter/exptime tokens read from stdin
poetry run group2master <group>           # dump a TAB.ALL object-database group as master-file lines
poetry run sfs <tel> yyyy/mm/dd [--evening|--morning] [--wind DEG] [--loud]   # sky flat field selector
```
`uobi_generator.py` (not wired as a poetry script) prints random 8-char UOBI ids: run it directly
with `python tpg/uobi_generator.py [n]`.

Run the test suite (pytest, in the `dev` dependency group):
```
poetry install --with dev
poetry run pytest
```
`tests/` covers the mask cascade, the allocator (including `--night_monitor`) and plan export.
Two things to know before adding tests: `tests/conftest.py` replaces the whole config after `Initiate()` so the
developer's `~/.config/tpg/local_config.yaml` cannot leak into results, and expectations that
depend on ephemerides are derived from `ob["calc"]` rather than hardcoded. There is no lint
config. For manual checks, run `tpg` against the `test` telescope config and inspect stdout/the
generated plan and log files.

## Configuration resolution

Two layered YAML configs, merged at runtime by `tpg/config.py:load_config()` (called from
`TelescopePlanGenerator._get_config()` and from sfs):
1. `tpg/config.yaml` — packaged defaults, one block per telescope (`master_file`,
   `start_makro`/`end_makro` paths, `data_file` archive dir, `overhead`, `hmin`/`hmax`, `hsun`,
   `wind_dist`, `moon_dist`, `def_read_mode`, `night_monitor_gap`), plus site coordinates and the
   `ctc` block (readout-mode mapping per telescope).
2. A local override file at the path in `local_config` (default `~/.config/tpg/local_config.yaml`),
   loaded if present and merged over the defaults via `deep_update` (recursive — a local
   telescope block only overrides the fields it lists; e.g. wk06's `dome` geometry lives only in
   `tpg/config.yaml`). `local_config.yaml` at the repo root is a personal example of this override file, not
   something the script reads directly.

`master_file` paths that don't resolve as given are retried relative to the package directory
(`tpg/`), so both absolute production paths and repo-relative paths work.

## Master files and observing blocks

`master_files/*.txt` are the tracked, hand-edited source of truth for observing targets (one
per telescope, e.g. `zb08_master.txt`, `jk15_master.txt`); commits like "updated exposure times"
or activating/deactivating targets are routine edits to these files, not code changes.
`master_files/{tel}_night_monitor.txt` hold the objects for `--night_monitor` mode and are passed
with `--master_file`; they are deliberately kept out of the regular master files, where they would
leak into normal plans as gap fillers even at `priority=0`. `makro/*.txt` are reusable
TOI-readable macro blocks a master-file line can splice in via `makro=path`. `plans/` is git-ignored — it holds generated nightly plans and their `_log.txt`
siblings, numbered `{tel}_{yyyymmdd}_{n}.txt`.

Each master-file line is parsed by `pyaraucaria.obs_plan.obs_plan_parser.ObsPlanParser` and
validated against schemas loaded via `pyaraucaria.ob_validator.ObsValidator` (`base_schema.yaml`
merged with a local `tpg_schema.yaml`, plus `base_rules.yaml`). The object-attribute vocabulary
(`seq`, `priority`, `uobi`, `cycle`, `ph_mk`, `t_start`/`t_end`, phase-fold keys, `makro`,
`ob_time`, etc.) is documented in README.md and is authoritative for what a line can contain —
read it before touching parsing/validation/masking code.

`pyaraucaria` and `ctc` (github.com/araucaria-project/{pyaraucaria,ctc}) are separate git
repositories pulled in as Poetry git dependencies; a checkout of `pyaraucaria` may also exist as
a sibling project directory during local development. `ctc` (`CycleTimeCalc`) is optional at
runtime — if its base folder / package isn't available, `telescope_plan_generator.py` falls back
to computing slot times manually from `seq` (see `ObsValidator.calc_seq_time`). `--no_ctc` forces
the manual path.

## Sky flat selector (sfs)

`tpg/sfs.py` (`SkyflatSelector`) picks a twilight flat field from `skyflat.fields_file` for the
evening or morning session (Sun between `skyflat.hsun` limits). It mirrors the tpg structure in
miniature — `go()` runs `Initiate` / `LoadFields` / `MakeTime` / `CalcFields`, then masks
(`MaskAltitude`, `MaskSun`, `MaskMoon`, `MaskWind`) narrowing per-field `visibility["all"]`
(numpy bool arrays over the session minutes), then `Rank` (whole-session visibility first, then
mean distance to the gradient null point 13 deg from zenith opposite the Sun). It deliberately
does not import `telescope_plan_generator.py` (astropy/pyaraucaria/ctc): shared config loading
(`load_config`, `deep_update`, `resolve_path`) lives in `tpg/config.py`, the GEM dome geometry in
`tpg/dome_eq.py`. `SelectFilters` (after `Rank`) picks one filter per `skyflat.groups` group —
oldest flats in the given history first — via the pure function `select_filters`; evening in group
order, morning reversed.

## Generation pipeline

All logic lives in `tpg/telescope_plan_generator.py` (`TelescopePlanGenerator`); `tpg/tpg.py` is
just the argparse CLI wrapper. `go()` runs a fixed, order-dependent sequence of steps, each
narrowing a per-object `visibility` array (one bool per minute of the night, from evening to
morning `hsun` crossing) stored under a named key plus a merged `"all"` key:

1. `Initiate` / `LoadObjects` / `parseObjects` — load config, parse/validate master-file lines,
   compute each object's `slotTime` (via `ob_time`, CTC, or manual `seq` calculation).
2. `MakeTime` — build the night's minute-resolution time grid (`ephem`-based sunset/sunrise).
3. `ObjectMask` → `CalcObject` → `MaskVisibility` → `MaskMoon` → `MaskWind` → `MaskTwilight` →
   `MaskCycle` → `MaskStartEnd` → `MaskPhaseStartEnd` → `MaskPhase` — each ANDs a new condition
   (tag/status/fwhm, altitude limits, moon distance/phase, wind vs telescope or (GEM `dome` block) dome azimuth, sunset/sunrise buffers,
   revisit `cycle` from archival light-curve data, absolute time windows, phase-fold windows,
   phase-coverage-based masking) into `visibility["all"]`.
4. `Waga` / `RandomizeList` — assign each object a weight favoring scarce visibility windows,
   then draw a random allocation order (seeded, reproducible via `--seed`). `Waga` raises
   `NoObservableObjects` (a `ValueError`) when nothing is observable; the CLI turns it into a
   one-line stderr message and exit code 1.
5. `allocate` — walk the night minute-by-minute; at each free minute, scan the random order for
   the highest-priority object that both fits contiguously in its remaining slot time and stays
   visible throughout, schedule it, and either remove it (no `cycle`) or re-mask it until
   `cycle` days later. Objects sharing a `uobi` are allocated/removed as one unit.
   With `--night_monitor` (implies tag `night_monitor`) `allocate` delegates to
   `_allocate_night_monitor`: objects stay in the pool, at each free minute one of the fitting
   highest-priority objects is drawn at random, and `night_monitor_gap` minutes after each
   observation are left free (export turns them into `WAIT ut=...`).
6. `export` — emit `WAIT ut=...` gaps and OBJECT lines (or spliced `makro=` file contents,
   filtered to schema-allowed keys) in plan order, then prepend/append `start_makro`/`end_makro`.
7. `SavePlan` — write the plan and its log to `plan_catalog`, auto-numbering per telescope/date.

When modifying the pipeline, preserve step order — later masks and the allocator assume earlier
steps have already populated `visibility["all"]` and `slotTime`.
