# tpg - telescope plan generator
version: 1.0.0

Script for generating observing plans for **wk06**, **zb08** and **jk15** telescopes. 
Script allocates objects in time based on priority, visibility and many other specifications. 
The slot duration time is calculated based on specification ``seq`` of each object. 
Script creates files readable by *TOI* and *TEXTUI*.

### RUNNING THE SCRIPT
```
poetry run tpg {wk06,zb08,jk15} date [ut] [--tag TAG] [--wind WIND] [--loud] [--seed] [--no_start_makro] [--no_end_makro]
date: YYYY/MM/DD [hh:mm:ss] 
```

```example
poetry run tpg zb08 2023/08/1  >| zb08_230810.txt
poetry run tpg wko6 2023/08/1 03:45:00 --no_start_makro
poetry run tpg wko6 2023/08/1 --seed 101010  --wind 124
poetry run tpg wko6 2023/08/1 --loud --tag stromgren
```
If `ut` is specified, the plan will **start at UT**.
If it's not specified, it will **start at sunset** (defined by `hsun` in ```config.yaml```)

| Argument           | desription                                                                             |
| ------------------ |----------------------------------------------------------------------------------------|
| `tel`              | telescope name: `wk06`, `zb08`, `jk15`.                                                |
| `date`             | observing date and ut `YYYY/MM/DD` or `YYYY/MM/DD hh:mm:ss`. If no ut, start at sunset |
| `--tag TAG`        | **tag** reqired for objects to allocate.                                               |
| `--wind WIND`      | wind **direction** to avoid (degrees).                                                 |
| `--seed SEED`      | random seed for object radmoization.                                                   |
| `--loud`           | enables “loud” mode — prints additional processing steps.                              |
| `--no_start_makro` | disables adding **makro** at the beginning of the plan.                                |
| `--no_end_makro`   | disables adding **makro** at the end of the plan.                                      |



### INSTALLATION
pull from github
```commandline
cd tpg
poetry install
```
upgrading
```commandline
poetry update
```

### CONFIGURATION
`config.yaml`

if local configuration file exist (in config.yaml by default: `~/.config/tpg/config.yaml`), 
specified parameters will be overwritten. 

```plainttext
wk06:                                                           # telescope
  master_file: /home/user/tpg/zb0806_master.txt                 # file with OB definitions
  start_makro: /home/user/tpg/master_files/zb08_start_makro.txt # start_makro
  end_makro: /home/user/tpg/master_files/zb08_end_makro.txt     # end makro
  data_file: /data/fits/wk06/processed-ofp/targets/             # directory with observations
  overhead: 15                                                  # overhead for time slot calculation
  hmin: 35                                                      # min altitude to observe
  hmax: 80                                                      # max altitude to observe
  hsun: -18                                                     # Sun alt to start and stop obs
  wind_dist: 95                                                 # wind distance to avoid in azimuth
  moon_dist: 40                                                 # distance to avoid moon
  def_read_mode: 2                                              # for slot time calculation
```      
### MASTER FILE (format examples)

`NAME RA DEC seq=<N/filter/exp_time>  priority=<priority> [optional_parameters]` 

`example_master.txt`
```
example1  21:24:00.2421   18:16:43.7   seq=1/Ic/30,1/V/30      priority=4    uobi=24j65   comment="binary and dusty envelope"
example2  07:24:14       -00:31:41     seq=3x(1/B/60,1/V/30)   priority=6    t_start=08:20:00    t_end=10:00:00
example3  07:24:14       -00:31:41     seq=10x(1/B/60)         priority=6    t_start=2025/11/23T08:20:00  t_end=2025/11/25T10:00:00
example4  07:24:14       -00:31:41     seq=2x(1/B/36)          priority=7    P=0.1   ph_start=0.8  ph_end=0.1  hjd0=2460310.5
example5  07:45:11.79    -36:58:29.1   seq=1/V/36              priority=7    ph_mk=2/Ic/0.05  P=13.7677045 
example6  07:24:14       -00:31:41     seq=2/Ic/3              priority=10   cycle=5 obs_data=home/user/obs_data/example5.txt  
example7  07:24:14       -00:31:41     seq=2/Ic/3.4            priority=11   makro=/home/user/makro/makro_ex10.txt  ob_time=300    
```

#### Specification:
* `seq=1/Ic/30,2/V/45` - (N/filter/exp_time) both sequences will be executed.  
  `seq=10x(1/Ic/30,2/V/45)` - (N/filter/exp_time) sequences will be repeated 10 times.   
* `priority=5`  - script allocates higher priority first. If priority is < 0, script will never allocate object
* `tag=photometric,standard` - if **tpg** is executed with `--tag TAG` option, only objects with this **TAG** will be allocated
* `uobi=34dfg` - unique observing block identifier. If some objects have **the same uobi**, all of them will be removed from the observing list if one of them will be allocated.
* `comment="some comment"` - comment visible to observer
* `cycle=3` - object will be observed if last observation was more than 3 days ago. Based on archival data. 

  `cycle=0.3` - object can be repeated during night 0.3*24h after last planed scheduling. accounts for archival data. 
*  `ph_mk=2/Ic/0.05` - parameters for phase masking:  no_of_obs/band_light_curve/phase_bin
*  `obs_data=/home/user/obs_data/example5.txt` - explicite specified data file for cycle or phase masking
* `t_start=2024/01/09T1:20:05` - date and ut for start

  `t_start=2024/01/09` - alternative date only 

  `t_start=1:20:05` - alternative ut only (if tpg start time is before 24:00 and object specified before 16:00, it will be considered +1 day)
* `ph_star=0.1 ph_end=0.3 P=7.5 hjd0=2460310.2` - phase to observe, all keys are required   
  `ph_star=0.8 ph_end=0.1 P=7.5 hjd0=2460310.2` - yes, it is cyclic
* `h_min=20` - object can be observed if more than 20 deg above horizon
* `h_max=70` - object can be observed if less than 70 deg above horizon
* `min_moon_dist=30` - object can be observed if more than 30 deg from the moon
* `max_moon_phase=70` - object can be observed if moon phase is less than 70 percent, and the moon is up.
* `ob_time=120` - slot time will be set as 120 seconds (not calculated by tpg)
* `makro=/home/user/makro_ex10.txt` - if scheduled file makro_ex10.txt will be included in the plan, instead of ob line
  `makro=/home/user/makro_ex10.txt priority=10 cycle=5 ob_time=300 obs_data=./test/obs_data/example5.txt` - should work
* `sciprog=cal`   - just for plan management    
* `pi=ocm` - just for plan management 

### REMARKS
* no white spaces in *NAME* *RA* *DEC* *seq* 
* If there is no object to observe, it will allocate:
```WAIT ut=12:45:34``` with the `ut` for the next object
* comment has to be enclosed with **" "** characters 
* `min_moon_dist`, `max_moon_phase` are calculated at the begining of the visibility of the object - only once per all night!
* phase masking: object is masked if in defined `phase_bin` are already `no_of_obs` points for `band_light_curve`
* for `cycle` full functionality, **obs data** are required, and for this the filter has to be specified. If `ph_mk` is not specified, it will take first filter defined in `seq`.
* since version 0.1.2 there is `makro` functionality. This replaces the OB entry with the contents of the file, allowing  
custom complex observing blocks, calibrations, or engineering sequences. The makro file should be **TOI readable** (OBJECT at the beginning), as it will be just inserted in the final plan. Previous function of splitting long (>10) sequences to multiple shorter sequences was removed. 
* since version 1.0.0 there are `start_makro` and `end_makro` that are inserted at the beginning and end of the plan. Both files should be **TOI readable**.
* since version 0.1.4 `uobi` parameter will not be exported if there is also `cycle` parameter for this ob.
* Objects sharing the same `uobi` are treated as a single observing block.
If any of them is allocated, all others are considered as observed according to theirs specification. 



### How it works

Each object in the master file contains a set of attributes such as coordinates, 
sequence definition (`seq`), priority, timing constraints, cycle limits, phase 
constraints, moon/wind limits, and optional macro substitution.

#### 1. Visibility computation
For the selected night and telescope location, the altitude of each object is
computed for every minute between the evening and morning Sun altitude limits
(`hsun`).  
From this, an initial **visibility array** is created:

* `True` — object is above `h_min` and below `h_max`
* `False` — outside altitude limits

This base visibility is then further **masked** using the object’s additional attributes:

* `cycle` — masked according to archival/plan data
* `t_start`, `t_end` — time window restriction (absolute or relative)
* `ph_start`, `ph_end`, `P`, `hjd0` — phase restrictions (cyclic-aware)
* `ph_mk` — phase-masking of already scheduled/archival data
* `min_moon_dist`, `max_moon_phase` — moon restrictions
* wind direction — masked if object lies within forbidden azimuth sector

The result is a **final visibility array** for each object.

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

1. At each moment, all **visible and eligible** objects are filtered.
2. Objects are grouped by **priority**.
3. Within the highest non-empty priority group, one object is selected **randomly**, with weights proportional to total visibility time (shorter visibility → higher probability).
4. If the object fits continuously in the visibility window for its full predicted slot time, it is scheduled.
5. After successful allocation:
   * Objects **without** `cycle` are removed from the candidate pool.
   * Objects **with** `cycle` remain, but their visibility array is masked until  
     `allocation_time + cycle * 24h`.
6. The time pointer moves to the end of the allocated slot.
7. If no object fits the current time slot, the scheduler inserts:
``WAIT ut=HH:MM:SS``
and jumps to the next possible time.

The process repeats until the morning Sun altitude limit is reached.

#### 4. Final plan construction
After allocation:

* optional `start_makro` and `end_makro` files are inserted  
(unless disabled via `--no_start_makro` or `--no_end_makro`)
* scheduled OB entries and macro expansions are written in TOI / TEXTUI format
* optional `uobi` fields are suppressed when `cycle` is present

#### Summary flowchart

Objects → Compute altitude → Base visibility  
→ Masking (cycle/phase/moon/wind/time)  
→ Predicted slot time  
→ Priority-based randomized scheduling  
→ Insert macros → Output plan





