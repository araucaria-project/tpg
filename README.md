# tpg - telescope plan generator
version: 1.1.1

Script for generating observing plans for **wk06**, **zb08**, **jk15** and **iris** telescopes. 
Script allocates objects in time based on priority, visibility and many other specifications. 
The slot duration time is calculated based on specification ``seq`` of each object. 
Script creates files readable by *TOI* and *TEXTUI*.

### RUNNING THE SCRIPT
```
poetry run tpg {wk06,zb08,jk15,iris,test} date [ut] [--tag TAG] [--wind WIND] [--fwhm FWHM] [--seed SEED] [--loud] [--no_save] [--no_ctc] [--no_start_makro] [--no_end_makro] [--master_file MASTER_FILE]
date: YYYY/MM/DD [hh:mm:ss] 
```

```example
poetry run tpg zb08 2023/08/1
poetry run tpg wk06 2023/08/1 03:45:00 --no_start_makro
poetry run tpg wk06 2023/08/1 --seed 101010  --wind 124
poetry run tpg wk06 2023/08/1 --loud --tag stromgren
```
If `ut` is specified, the plan will **start at UT**.
If it's not specified, it will **start at sunset** (defined by `hsun` in ```config.yaml```)

By default the generated plan is saved automatically to `plan_catalog` (see `config.yaml`) as
`{tel}_{yyyymmdd}_{n}.txt`, with a matching `{tel}_{yyyymmdd}_{n}_log.txt` log file; use `--no_save`
to skip this. stdout is **not** a clean plan file — it interleaves the plan's `WAIT`/`OBJECT` lines
with diagnostic messages (config loaded, sunset/sunrise, random seed, makro status, etc.), so
redirecting it (e.g. `>| plan.txt`) will not produce a valid TOI/TEXTUI plan.

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
shallow — a telescope key present locally replaces that telescope's **whole** block.

```plainttext
plan_catalog: ~/tpg_plans/                                      # where plans and logs are saved

zb08:                                                           # telescope
  master_file: /home/user/tpg/master_files/zb08_master.txt      # file with OB definitions
  start_makro: /home/user/tpg/master_files/zb08_start_makro.txt # start_makro
  end_makro: /home/user/tpg/master_files/zb08_end_makro.txt     # end makro
  data_file: /data/fits/zb08/processed-ofp/targets/             # observations root, also root for obs_data
  overhead: 15                                                  # overhead for time slot calculation
  hmin: 35                                                      # min altitude to observe
  hmax: 80                                                      # max altitude to observe
  hsun: -18                                                     # Sun alt to start and stop obs
  wind_dist: 95                                                 # wind distance to avoid in azimuth
  moon_dist: 40                                                 # distance to avoid moon
  def_read_mode: 2                                              # for slot time calculation
```      
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
| wind | object inside the forbidden azimuth sector around `--wind` |
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





