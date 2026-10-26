# tpg - telescope plan generator
version: 1.1.1

Script for generating observing plans for **wk06**, **zb08**, **jk15** and **iris** telescopes. 
Script allocates objects in time based on priority, visibility and many other specifications. 
The slot duration time is calculated based on specification ``seq`` of each object. 
Script creates files readable by *TOI* and *TEXTUI*.

### RUNNING THE SCRIPT
```
poetry run tpg {wk06,zb08,jk15,iris,test} date [ut] [--tag TAG] [--wind WIND] [--fwhm FWHM] [--seed SEED] [--loud] [--no_save] [--no_ctc] [--no_start_makro] [--no_end_makro] [--master_file MASTER_FILE] [--night_monitor] [--skyflats]
date: YYYY/MM/DD [hh:mm:ss] 
```

```example
poetry run tpg zb08 2023/08/1
poetry run tpg wk06 2023/08/1 03:45:00 --no_start_makro
poetry run tpg wk06 2023/08/1 --seed 101010  --wind 124
poetry run tpg wk06 2023/08/1 --loud --tag stromgren
poetry run tpg zb08 2023/08/1 --night_monitor --master_file master_files/zb08_night_monitor.txt
```
If `ut` is specified, the plan will **start at UT**.
If it's not specified, it will **start at sunset** (defined by `hsun` in ```config.yaml```)

By default the generated plan is saved automatically to `plan_catalog` (see `config.yaml`) as
`{tel}_{yyyymmdd}_{n}.txt`, with a matching `{tel}_{yyyymmdd}_{n}_log.txt` log file; use `--no_save`
to skip this. stdout is **not** a clean plan file — it interleaves the plan's `WAIT`/`OBJECT` lines
with diagnostic messages (config loaded, sunset/sunrise, random seed, makro status, etc.), so
redirecting it (e.g. `>| plan.txt`) will not produce a valid TOI/TEXTUI plan.

If nothing can be allocated, tpg prints a one-line reason to stderr and exits with code 1, without
saving a plan: either no object in the master file has the required tag (`--tag` / `--night_monitor`),
or no object survives the visibility masks this night (`--loud` shows which mask drops them).

| Argument                    | desription                                                                             |
| --------------------------- |----------------------------------------------------------------------------------------|
| `tel`                       | telescope name: `wk06`, `zb08`, `jk15`, `iris`, `test`.                                |
| `date`                      | observing date and ut `YYYY/MM/DD` or `YYYY/MM/DD hh:mm:ss`. If no ut, start at sunset |
| `--tag TAG`                 | **tag** reqired for objects to allocate.                                               |
| `--wind WIND`               | wind **direction** to avoid (degrees).                                                 |
| `--fwhm FWHM`               | current **fwhm**, used to apply objects' `fwhm_min`/`fwhm_max` limits.                |
| `--seed SEED`               | random seed for object radmoization.                                                   |
| `--loud`                    | enables “loud” mode — prints additional processing steps.                              |
| `--no_save`                 | do not save the generated plan to `plan_catalog`.                                      |
| `--no_ctc`                  | do not use **CTC** for slot-time calculation, even if available.                       |
| `--no_start_makro`          | disables adding **makro** at the beginning of the plan.                                |
| `--no_end_makro`            | disables adding **makro** at the end of the plan.                                      |
| `--master_file MASTER_FILE` | override the master file configured for `tel`.                                         |
| `--night_monitor`           | monitoring mode: only `tag=night_monitor` objects, repeated, with `night_monitor_gap` minutes break (see below). |
| `--skyflats`                | add the evening/morning `SKYFLAT` chosen by sfs (see "SKYFLAT SELECTOR" below). |



### INSTALLATION
Requires Python >=3.10,<3.14 (Python 3.14 is not yet supported: `scipy` does not ship a
prebuilt wheel for it, which forces a from-source build requiring a Fortran compiler and can
fail on a fresh machine). Once `scipy` publishes 3.14 wheels, the cap in `pyproject.toml` can
be raised.

pull from github
```commandline
cd tpg
poetry install
```
upgrading
```commandline
poetry update
```

### TESTS
```commandline
poetry install --with dev
poetry run pytest
```
The suite covers the mask cascade, the allocator and plan export. It is self-contained: every
test builds its own master file and light-curve data under a temporary directory, and the
config is replaced after `Initiate()` so a local `~/.config/tpg/local_config.yaml` cannot
change the results.

**Planned:** `pyproject.toml` will be moved to the standard `[project]` layout (PEP 621) and the
project will switch from **poetry** to **uv**. `poetry check` already reports the current
`[tool.poetry]` metadata and `[tool.poetry.scripts]` as deprecated. Until that happens the commands
above stay valid; afterwards they become `uv sync` / `uv run tpg ...`.

### CONFIGURATION
`config.yaml`

if local configuration file exist (path given by `local_config` in config.yaml, by default:
`~/.config/tpg/local_config.yaml`), specified parameters will be overwritten. The merge is
key by key at every level — a local file only needs the fields it changes (typically the
paths); everything else, e.g. `dome`, comes from `config.yaml`.

```plainttext
plan_catalog: ~/tpg_plans/                                      # where plans and logs are saved

zb08:                                                           # telescope
  master_file: /home/user/tpg/master_files/zb08_master.txt      # file with OB definitions
  start_makro: /home/user/tpg/master_files/zb08_start_makro.txt # start_makro
  end_makro: /home/user/tpg/master_files/zb08_end_makro.txt     # end makro
  data_file: /data/fits/zb08/processed-ofp/targets/             # observations root, also root for obs_data
  overhead: 15                                                  # SECONDS added to each exposure
  hmin: 35                                                      # deg, min altitude to observe
  hmax: 80                                                      # deg, max altitude to observe
  hsun: -18                                                     # deg, Sun alt to start and stop obs
  wind_dist: 95                                                 # deg, wind sector to avoid in azimuth
  moon_dist: 40                                                 # deg, default distance to avoid moon
  def_read_mode: 2                                              # for slot time calculation
  night_monitor_gap: 30                                         # min, break between observations in --night_monitor mode (default 30)
```

Telescopes on a German equatorial mount (currently only wk06) also get a `dome` block. The dome
slit is then not at the telescope azimuth, and the wind mask (`--wind`, `wind_dist`) is applied
to the **dome** azimuth instead. The side of pier is not known in advance, so a minute is masked
when the dome azimuth for *either* side falls into the wind sector. Without a `dome` block the
dome is assumed to follow the telescope azimuth (alt-az telescopes).

```plainttext
wk06:
  dome:
    mount: gem          # dome azimuth from GEM geometry, both sides of pier
    r_dome: 2050        # mm, dome radius
    spx: -110           # mm, dome center -> ra/dec axes intersection, (+) E / (-) W
    spy: -110           # mm, dome center -> ra/dec axes intersection, (+) N / (-) S
    gem: 670            # mm, ra/dec axes intersection -> optical axis
```

The geometry is a temporary copy of `pyaraucaria.dome_eq` in `tpg/dome_eq.py` (the pyaraucaria
version only computes for the current moment); once pyaraucaria accepts a time it should be
used instead and the copy removed.

`overhead` is **seconds per exposure** (readout and similar), not minutes per observing block —
it is passed to `ObsValidator.calc_seq_time` as its per-exposure term. The only fixed per-block
cost in the manual path is a hardcoded 5 s; slew between targets is not modelled (see ROADMAP).

### MASTER FILE (format examples)

`OBJECT NAME RA DEC seq=<N/filter/exp_time>  priority=<priority> [optional_parameters]` 

`example_master.txt`
```
OBJECT example1  21:24:00.2421   18:16:43.7   seq=1/Ic/30,1/V/30      priority=4    uobi=24j65   comment="binary and dusty envelope"
OBJECT example2  07:24:14       -00:31:41     seq=3x(1/B/60,1/V/30)   priority=6    t_start=08:20:00    t_end=10:00:00
OBJECT example3  07:24:14       -00:31:41     seq=10x(1/B/60)         priority=6    t_start=2025/11/23T08:20:00  t_end=2025/11/25T10:00:00
OBJECT example4  07:24:14       -00:31:41     seq=2x(1/B/36)          priority=7    P=0.1   ph_start=0.8  ph_end=0.1  hjd0=2460310.5
OBJECT example5  07:45:11.79    -36:58:29.1   seq=1/V/36              priority=7    ph_mk=2/Ic/0.05  P=13.7677045 
OBJECT example6  07:24:14       -00:31:41     seq=2/Ic/3              priority=10   cycle=5 obs_data=example6/Ic/example6_Ic.txt  
OBJECT example7  07:24:14       -00:31:41     seq=2/Ic/3.4            priority=11   makro=/home/user/makro/makro_ex10.txt  ob_time=5    
```

#### Specification:
* `seq=1/Ic/30,2/V/45` - (N/filter/exp_time) both sequences will be executed.  
  `seq=10x(1/Ic/30,2/V/45)` - (N/filter/exp_time) sequences will be repeated 10 times.   
* `priority=5`  - script allocates higher priority first. If priority is < 0, script will never allocate object
* `tag=photometric,standard` - if **tpg** is executed with `--tag TAG` option, only objects with this **TAG** will be allocated
* `uobi=34dfg` - unique observing block identifier. If some objects have **the same uobi**, allocating one of them retires all of them: those without `cycle` are removed from the observing list, those with `cycle` are masked for their `cycle` time (as if they had been observed).
* `comment="some comment"` - comment visible to observer
* `cycle=3` - object will be observed if last observation was more than 3 days ago. Based on archival data. 

  `cycle=0.3` - object can be repeated during night 0.3*24h after last planed scheduling. accounts for archival data.
  The in-night interval is counted from the **start** of the previous slot, so for `cycle` shorter than
  the slot time the object can be scheduled again right after itself.
*  `ph_mk=2/Ic/0.05` - parameters for phase masking:  no_of_obs/band_light_curve/phase_bin.
   `phase_bin` is a **half-width**: the window checked around the current phase is `phase ± phase_bin`,
   and it **wraps across the 0/1 seam** — phases 0.99 and 0.01 are 0.02 apart, not 0.98
*  `obs_data=example6/Ic/example6_Ic.txt` - explicite specified data file for cycle or phase masking.
   The path is always resolved **relative to `data_file`** of the telescope (see `config.yaml`);
   a leading `/` does not make it absolute. Without `obs_data` the path is built automatically as
   `<data_file>/<name>/<filter>/light-curve/<name>_<filter>_diff_light_curve.txt` (lowercased name)
* `t_start=2024/01/09T1:20:05` - date and ut for start

  `t_start=2024/01/09` - alternative date only 

  `t_start=1:20:05` - alternative ut only. The time is first attached to the **calendar date of the plan start**;
  if the result falls before the plan start, one day is added. So an ut smaller than the start ut always
  means "the coming night", never "yesterday"
* `ph_start=0.1 ph_end=0.3 P=7.5 hjd0=2460310.2` - phase to observe, all keys are required   
  `ph_start=0.8 ph_end=0.1 P=7.5 hjd0=2460310.2` - yes, it is cyclic

  The phase window needs no observational data — it is computed from `P` and `hjd0` alone.
  If `hjd0` is omitted while `P` is given, `hjd0=2460000` is assumed
* `h_min=20` - object can be observed if more than 20 deg above horizon
* `h_max=70` - object can be observed if less than 70 deg above horizon
* `min_moon_dist=30` - object can be observed if more than 30 deg from the moon
* `max_moon_phase=70` - object can be observed if moon illumination is less than 70 percent.
  Only the illuminated fraction is checked — the Moon's altitude is **not** taken into account,
  so a bright but already set Moon still masks the object. Use `min_moon_dist` for the geometry.
* `sunset_dt=1` - object can be observed only 1h or more after sunset
* `sunrise_dt=1` - object can be observed only until 1h before sunrise

  "sunset"/"sunrise" here mean the Sun crossing `hsun` (by default -18 deg), i.e. the same limits that
  define the beginning and the end of the night
* `fwhm_min=1.5` - object can be observed if current fwhm (given with `--fwhm`) is more than 1.5
* `fwhm_max=3` - object can be observed if current fwhm (given with `--fwhm`) is less than 3
* `status=finished` or `status=inactive` - object will never be allocated, regardless of visibility
* `ob_time=12` - slot time will be set as 12 minutes (not calculated by tpg)
* `makro=/home/user/makro_ex10.txt` - if scheduled file makro_ex10.txt will be included in the plan, instead of ob line
  `makro=/home/user/makro_ex10.txt priority=10 cycle=5 ob_time=5 obs_data=./test/obs_data/example5.txt` - should work
* `read_mod=2` - readout mode, passed to **CTC** for slot-time calculation and written to the plan
* `sciprog=cal`   - just for plan management    
* `pi=ocm` - just for plan management 

### REMARKS
* no white spaces in *NAME* *RA* *DEC* *seq* 
* If there is no object to observe, it will allocate:
```WAIT ut=12:45:34``` with the `ut` for the next object
* comment has to be enclosed with **" "** characters 
* `min_moon_dist`, `max_moon_phase` are evaluated **for every minute** of the night, so the Moon
distance mask follows the Moon across the sky instead of being fixed for the whole night
* phase masking: object is masked if in defined `phase_bin` are already `no_of_obs` points for `band_light_curve`
* `cycle` works on two independent levels:
  * **against archival data** — needs **obs data**. The light-curve path is built from the **first filter
  of `seq`** (`ph_mk` is *not* consulted here), unless `obs_data` points at the file explicitly.
  If the file is missing or unreadable, this part is skipped and the object is treated as never observed.
  * **within the night** — after the object is scheduled, its visibility is masked for `cycle * 24h`
  counted from the **start** of the slot. This always works, no data file needed.
* since version 0.1.2 there is `makro` functionality. This replaces the OB entry with the contents of the file, allowing  
custom complex observing blocks, calibrations, or engineering sequences. The makro file should be **TOI readable** (OBJECT at the beginning), as it will be just inserted in the final plan. Previous function of splitting long (>10) sequences to multiple shorter sequences was removed. 
* since version 1.0.0 there are `start_makro` and `end_makro` that are inserted at the beginning and end of the plan. Both files should be **TOI readable**.
* plan lines contain only the attributes allowed by the OBJECT schema (`name`, `ra`, `dec`, `seq`,
`read_mod`, `uobi`, `sciprog`, `pi`, `comment`, ...). Scheduling attributes — `priority`, `cycle`,
`t_start`/`t_end`, `ph_*`, `P`, `hjd0`, `h_min`/`h_max`, `min_moon_dist`, `max_moon_phase`,
`fwhm_*`, `status`, `tag`, `obs_data`, `ob_time` — are used by tpg but are **not** written to the plan.
* Objects sharing the same `uobi` are treated as a single observing block.
If any of them is allocated, all others are considered as observed according to theirs specification. 



### How it works

Each object in the master file contains a set of attributes such as coordinates, 
sequence definition (`seq`), priority, timing constraints, cycle limits, phase 
constraints, moon/wind limits, and optional macro substitution.

#### 1. Visibility computation
The night is a grid of minutes, from the evening to the morning crossing of the Sun altitude
limit `hsun`. Every object carries a **visibility array** over that grid — one boolean per minute —
which starts out all `True`.

The whole selection stage is a **cascade of masks**: each step ANDs one more condition into the
array, so **every step can only narrow it, never widen it**. In the order in which they are applied:

| step | narrows the array by |
| ---- | -------------------- |
| `tag` / `status` / `fwhm_min` / `fwhm_max` | drops the object for the whole night if it does not qualify |
| altitude | object must be above `h_min` and below `h_max` |
| `min_moon_dist`, `max_moon_phase` | Moon distance (per minute) and Moon illumination |
| wind | object (or, with a `dome` block, dome for either side of pier) inside the forbidden azimuth sector around `--wind` |
| `sunset_dt`, `sunrise_dt` | buffer after the evening / before the morning `hsun` crossing |
| `cycle` | too soon after the last archival observation |
| `t_start`, `t_end` | absolute time window |
| `ph_start`, `ph_end`, `P`, `hjd0` | phase window (cyclic-aware) |
| `ph_mk` | phase bins already covered by existing data |

Each condition is also kept under its own key, next to the merged `"all"` array that the scheduler
actually uses. With `--loud` the log prints, per object, how many hours are left under each key —
the quickest way to see *which* mask killed a target. (Most keys hold that one condition alone;
`t_start`/`t_end`/`ph_*` hold the array as already narrowed up to that step.)

```
OBJECT example1 [28m]  all: 0.0  h_min: 0.0  h_max: 10.2  min_moon_dist: 4.6  ...
```

The merged `"all"` array is the **final visibility array** of the object.

#### 2. Slot-time calculation
For every object, the **predicted execution time** (slot length) is computed
from the sequence specification (`seq`), using [CTC](https://github.com/araucaria-project/ctc) if available. 
If CTC is not available, it tpg will calculate it ''manually''.  
This includes:

* exposure times  
* number of cycles  
* readout modes  
* telescope overheads  

Optional fixed override via `ob_time`.

Objects with defined `makro=` do not generate an OB entry — instead, the macro
file is inserted directly into the plan, and `ob_time` is strongly suggested. 

#### 3. Scheduling and allocation
The main allocation loop advances minute-by-minute through the night.

1. Before the loop starts, a **random order** of candidates is drawn **once** (reproducible with
`--seed`), with weights favouring objects with short total visibility. Objects that are never
visible during the night get weight 0 and are dropped.
2. At each free minute, that fixed order is scanned and the **highest-priority** object that still
fits is chosen (ties resolved by the drawn order).
3. If the object fits continuously in the visibility window for its full predicted slot time, it is scheduled.
4. After successful allocation:
   * Objects **without** `cycle` are removed from the candidate pool.
   * Objects **with** `cycle` remain, but their visibility array is masked until  
     `allocation_time + cycle * 24h`, where `allocation_time` is the **start** of the slot.
5. The allocated minutes are marked as busy, so the loop simply walks on to the first minute that is
still free.
6. If no object fits the current time slot, the scheduler inserts:
``WAIT ut=HH:MM:SS``
and jumps to the next possible time.

The process repeats until the morning Sun altitude limit is reached.

**`--night_monitor` mode** replaces steps 2–5: only objects tagged `night_monitor` are considered
(the mode implies `--tag night_monitor`), and they are **not** removed after allocation. At each
step the candidates that fit are narrowed to the highest priority and one of them is drawn
at random (re-drawn every step, reproducible with `--seed`). After each observation the next
`night_monitor_gap` minutes (per-telescope config key, default 30) are left free, which the
export turns into `WAIT ut=...`. `cycle` still applies: an object with `cycle` drops out for
`cycle` days after its first allocation.

Keep night-monitor objects in a **separate master file** (`master_files/{tel}_night_monitor.txt`) and
pass it with `--master_file`. In the regular master file they would be ordinary candidates in normal
runs too — even at `priority=0` they fill slots where nothing else fits, and their mere presence
changes the random draw for a given `--seed`.

#### 4. Final plan construction
After allocation:

* the plan opens with `WAIT sunset=<hsun>` when it starts at sunset, or with `WAIT ut=HH:MM:SS`
when an explicit ut was given
* optional `start_makro` and `end_makro` files are inserted  
(unless disabled via `--no_start_makro` or `--no_end_makro`)
* scheduled OB entries and macro expansions are written in TOI / TEXTUI format
* attributes not allowed by the OBJECT schema (`priority`, `cycle`, `t_start`, ...) are dropped

#### Summary flowchart

Objects → Predicted slot time 
→ Compute altitude and visibility
→ Masking (cycle/phase/moon/wind/time)   
→ Priority-based randomized scheduling  
→ Insert macros → Output plan



### SKYFLAT SELECTOR (sfs)

`sfs` picks the field for twilight sky flats, separately for the evening and the morning session:
```
poetry run sfs {wk06,zb08,jk15,iris,test} yyyy/mm/dd [--evening | --morning] [--wind DEG] [--fields_file PATH] [--history FILTER=DAYS ...] [--loud]
```
The date is the night (evening) date, as in tpg. Without `--evening`/`--morning` both sessions are
printed. stdout gets only the plan line `SKYFLAT name ra dec seq=...` per session (see below), ready to
paste or pipe; warnings go to stderr. `--loud` adds a `#` summary line (session window, null point
distance) and a `# filters:` line per session, plus every field with its masks and score.
If a session has no visible field, sfs reports it on stderr, still computes the other session and
exits with code 1 at the end.

* **Session**: the Sun between the `skyflat.hsun` limits (default +1 -> -15 deg in the evening,
-15 -> +1 in the morning), 1-minute grid.
* **Date**: the night (evening) date — `2026/10/26 --morning` is the morning of 27 Oct. sfs warns
when the requested session is already over or starts more than 24 h from now (typically today's
calendar date passed after midnight, which plans the next day's morning).
* **Masks** (every minute of the session, a field must pass all of them for the whole session):
field altitude within `hmin`-`hmax` (currently 40-78 deg), azimuth at least `sun_az_dist`
(default 100 deg) from the Sun azimuth, `moon_dist` from the Moon, and with `--wind` the same
`wind_dist` check as tpg (dome azimuth for both sides of pier on a GEM with a `dome` block).
* **Ranking**: mean angular distance from the twilight sky gradient null point, 13 deg from the
zenith on the anti-solar vertical (Chromey & Hasselbacher 1996, PASP 108, 944), i.e. ~77 deg
altitude opposite the Sun; with `hmax` below that it favours the highest field opposite the Sun.
* **Priorities**: a field line may have `priority=N` (higher first, default 0; for now all fields
are equal). Among the fields visible for the whole session the highest priority wins, the null
point distance decides only within one priority; a lower priority is used only when no field of
a higher one passes the masks.
* If no field passes the masks for the whole session, the one visible for the longest part of it
is chosen (priority, then distance break ties) and marked as fallback (`WARNING` in the output).
* **Filters**: `skyflat.groups` (per telescope) is an ordered list of filter groups, from the
filter needing the most light (darkest flats) to the least, e.g. for jk15
`[[u_s, Ha_n, Hb_n], [u, z, v_s, b_s, y_s, Hb_w, Ha_w], [B, V, Ic, g, r, i]]`; the number of groups is the number of
filters one session can do on that telescope. From every group sfs takes the filter with the
oldest flats (never taken = oldest; ties keep the group order). The evening goes in group order
(the sky gets darker), the morning in reverse.
* **Flat history** (`--history u=3.5 V=12 ...`): days since the last flats per filter, for the
session being planned (evening and morning histories are kept apart by the caller — run
`--evening`/`--morning` separately with their own history). A filter missing from the history, or
given as `none`, is treated as never taken; without `--history` the first filter of every group is
taken. In Python the history is a plain dict `{filter: days | None}`; converting from other
sources (e.g. flat_overwatch timestamps) is up to the caller. With `--loud` the output gets a
`# filters: ...` line in observing order with the age of each filter's flats.
* **Full sequences in groups**: a group entry may be a whole sequence for its filter instead of a
bare name, e.g. `groups: [[u_s, Ha_n], [5/Ic/3, V]]` — if `Ic` is chosen, `5/Ic/3` goes into the
plan as given; bare names get the defaults below. Flat history is still looked up by the filter
name (`Ic`). A malformed entry (not `filter` or `N/filter/exp`) stops sfs with an error.
* **Plan line**: `SKYFLAT <field> <ra> <dec> seq=7/u/a,7/z/a,...` — `skyflat.n` frames per filter
(default 7) with `skyflat.exp` exposure (default `a`, auto exposure: TOI waits for the right sky
level and computes the exposure itself). Without groups the line has no `seq`.

Configuration: defaults in the top-level `skyflat:` block of `config.yaml` (`hsun`, `hmin`,
`hmax`, `sun_az_dist`, `n`, `exp`, optionally `moon_dist`/`wind_dist`, which otherwise come from the
telescope block); a telescope block may override single fields in its own `skyflat:` block, and
holds its `groups` (no groups = field only, `SKYFLAT` line without `seq`). `fields_file`
has no default and **must** be set in every telescope's `skyflat:` block (or given with
`--fields_file`). It may be absolute or relative to the package / repo root; all telescopes
currently use `master_files/sky_flat_fields.txt`, one `SKYFLAT name ra dec [priority=N]` line per field.

**In tpg plans** (`--skyflats`, or `TelescopePlanGenerator(..., skyflats=True, flat_history={"evening": {...},
"morning": {...}}, now=None)`): the evening `SKYFLAT` goes before the night part of the plan (`WAIT sunset=...` /
`WAIT ut=...`), the morning one after it, and `start_makro`/`end_makro` still wrap everything (a makro may open or
close the dome). A session that is over at the moment the plan should start from is skipped — that moment is the
`date`/`ut` argument (with no ut: the day, so a plan "from sunset" gets the evening flats), or `now` if given and
later (TOI passes the current time). A session without a field is skipped with a `WARNING`, the plan is still made; with no `skyflat.groups`
for the telescope no `SKYFLAT` is added (it would have no `seq`).
The outcome per session is in `p.skyflat_info` (`{"evening": {"line": ... | None, "note": ...}, ...}`); tpg uses the
same `--wind` for sfs. The night date for sfs comes from the plan's own sunrise, not from the argument.

**Flat session makros** (`skyflat.evening_makro` / `skyflat.morning_makro` in the telescope's `skyflat` block,
paths like `fields_file`): with `--skyflats` the whole session comes from the makro instead of a single line — e.g.
`makro/wk06_evening_makro.txt`:
```
# lines starting with # are not copied to the plan
WAIT sunset=-1
BELL
SKYFLAT
BELL
WAIT sunset=-16
BELL
ZERO seq=11/u/0
BELL
```
The bare `SKYFLAT` line is replaced by the sfs choice; a skipped session leaves out the whole makro (no `WAIT
sunset=-1` in a plan made at night). A missing makro file or one without a bare `SKYFLAT` line falls back to the single
`SKYFLAT` line with a `WARNING` in `skyflat_info`. Order: `start_makro`, evening makro, night part, morning makro,
`end_makro`. Current makros: wk06/zb08 wait for `sunset=-16`, jk15 for `-14` before the zeros.

From Python (e.g. TOI), one call per session:
```python
from tpg.sfs import SkyflatSelector, NoSkyflatField

r = SkyflatSelector("zb08", "2026/10/26", "evening", wind=None,
                    flat_history={"u": 3.5, "V": 12, "Ic": None}).go()
r["line"]       # 'SKYFLAT BF22b 22:23:50.0 -30:58:00.0 seq=7/u/a,...' - plan line
r["filters"]    # [{"filter": "u", "entry": "u", "age": 3.5}, ...] in observing order
r["fallback"]   # True if no field is visible for the whole session
```
Also in the result: `name`, `ra`, `dec`, `priority`, `field_line`, `seq`, `session`, `start`/`end`
(`ephem.Date`), `score`, `fraction`. Warnings (fallback, missing history, suspicious date) are
collected in the selector's `msg`. `NoSkyflatField` (a `ValueError`) is raised when no field is
visible at any minute of the session; a malformed config raises `ValueError`.

### ROADMAP

Ideas and known limitations, not commitments. Ordered roughly by return on effort.
For scale: of the 917 active OBJECT lines in `master_files/`, 596 use `cycle`, 540 use `ph_mk`
and 305 use `ph_start` — the phase and cadence machinery is the core of this catalogue, so that
is where most of the value sits.

#### Observing model

* **Slew and dome rotation are not modelled.** `overhead` is per-exposure readout, and the only
fixed per-block cost is a hardcoded 5 s, so the scheduler effectively assumes the telescope moves
between any two points of the sky instantly. Over a night the plan drifts late.
`CycleTimeCalc` already exposes `set_telescope_start_az_alt`, `set_dome_start_az`,
`set_start_time` and `set_observatory_location`, and tpg calls none of them — it calls
`reset_time()` before every block and computes each one from an unspecified position.
Wiring this up is mostly plumbing, but it has an architectural cost: `slotTime` stops being a
property of an object and becomes a function of the *pair* (previous target, next target), so it
can no longer be computed once at parse time — the allocator has to evaluate it as it goes.
* **The Moon's altitude is ignored.** `max_moon_phase` looks only at the illuminated fraction, so
a bright Moon that has already set still masks the target for the rest of the night. The cheap fix
is to apply the mask only while the Moon is above the horizon. The proper fix is to replace both
`min_moon_dist` and `max_moon_phase` with one physical quantity — sky brightness at the target
(e.g. Krisciunas & Schaefer 1991), which combines separation, phase, Moon altitude and target
altitude. Telling: only 3 objects set `min_moon_dist` explicitly, which suggests the present knobs
are not the ones anyone wants to turn.
* **Phase coverage is binary.** `ph_mk` either masks a minute or not. Weighting by *coverage gain*
— a point in an empty phase bin being worth more than one in a nearly full bin — would reuse the
existing `Waga` mechanism rather than introduce a new concept.
* **No multi-night planning.** Every night is scheduled independently and `cycle` is a coarse
substitute for cadence. A period of 13.7 d cannot be covered in one night. A "plan the next N
nights, maximising phase coverage" mode would be the largest gain in scientific value here, and
also the largest piece of work.
* **Bright targets cannot use twilight.** The minute grid *starts* at the `hsun` crossing, so a
per-object `hsun` would change nothing until the grid boundaries themselves are derived from the
earliest `hsun` among the objects. Worth 30–60 min a night for standards and calibrators.
* **Altitude limits stand in for airmass.** `h_min`/`h_max` are a proxy for what actually matters.
`pyaraucaria` already ships `airmass.py`, so an `airmass_max` attribute would be wiring, not physics.
* **Time scales are mixed.** `hjd0` is named heliocentric and `jd_obs` is geocentric, and they are
subtracted directly. Harmless for one site and long periods; the heliocentric correction reaches
±8.3 min, which is 0.06 in phase at `P = 0.1 d` — larger than a typical `phase_bin`. This becomes
a correctness issue as soon as archival data from several sources is mixed.

#### Architecture

* **An observation-history provider.** After the 2026-08-12 clean-up the whole dependency on
archival data is two masks asking one question: *give me the JDs of good past observations of this
target in this band*. Behind a small protocol that could be the current OFP directory layout, a
database, a TAP query, or an external survey — which is what would make `cycle` mean "since anyone
observed it" rather than "since we did, and only if the file is where we expect". It also makes
`cycle` and `ph_mk` testable without files on disk. A provider must declare which time scale it
returns (see above).
* **The cascade as a list of functions.** Each mask is essentially `(ob, grid, ctx) -> ndarray`.
A registry instead of a hardcoded `go()` would let a new mask be added without touching the
pipeline, and is the natural extension point for anyone reusing this.
* **No `sys.exit()` in library code.** `LoadObjects` and `_parse_date` kill the process. The class
is already used as a library (by the test suite), so these should raise and the CLI wrapper should
translate to exit codes — which also fixes `tpg` always exiting with status 1.
* **Smaller, but real:** two bare `except:` clauses remain — a third, in `MaskPhaseStartEnd`, was
what silently disabled the phase window until it was removed. The one left in `MaskPhase` is
unreachable for schema-valid input, since a malformed `ph_mk` is already rejected by the validator.
Per-key visibility arrays are inconsistent (some hold the isolated condition, some the cumulative
one, which muddies `--loud` diagnostics); `visibility["all"]` starts as a list and silently becomes
an `ndarray` halfway through the cascade; `d_ver` treats `0` as a missing key, which is currently
harmless only because the parser returns strings.

#### Housekeeping

* CI running `pytest` on push — the suite exists, the workflow file does not.
* `pyproject.toml` to the `[project]` layout and poetry to uv (see INSTALLATION).
* `poetry.lock` is gitignored, so two machines can resolve different versions. Worth deciding
deliberately when the packaging is touched.

#### Deliberately not planned

Rewriting the allocator as a solver (ILP/CP-SAT). The greedy walk with weighted randomisation is
**explainable**: when an observer asks why a target did not make it, `--loud` answers with hours
remaining per mask key. A solver would buy a few percent of packing efficiency and take that
answer away. With 243 objects sharing `priority=5`, the weighted draw already does real work
inside each priority class.





