# tpg - Telescope Plan Generator

Script for generating observing plans for **wk06**, **zb08** and **jk15** telescopes. Script alocates time based on priority and visibility. Calculates slot time based on specification ``seq`` of each object. Creates files readable by *TOI* and *TEXTUI*.

### RUNNING SCRIPT
```
poetry run tpg {wk06,zb08,jk15} date [ut] [--wind WIND] [-h]
date: YYY/MM/DD [hh:mm:ss] 
```

```example
poetry run tpg zb08 2023/08/1  >| zb08_230810.txt
poetry run tpg wko6 2023/08/1 03:45  --wind 124
```
If `ut` is specified, the plan will **start at UT**.
If it's not specified, it will **start at sunset** (defined by `hsun` in ```config.yaml```)

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

```plainttext
wk06:                                                        # telescope
  master_file: /home/observer/src/tpg/wk06_master.txt        # file with OB definitions
  data_file: /data/fits/wk06/processed-ofp/targets/          # directory with observations
  overhead: 15                                               # overhead for time slot calculation
  hmin: 35                                                   # min altitude to observe
  hmax: 91                                                   # max altitude to observe
  hsun: -18                                                  # Sun alt to start and stop obs
  wind_dist: 95                                              # wind distance to avoid in azimuth
  moon_dist: 40                                              # distance to avoid moon
```      
### MASTER FILE (format)
`zb08_master.txt`
```
# NAME RA DEC seq=<N/filter/exp_time>  priority=<priority> [optional: comment="<comment for observer>" t_start=2024/01/09T1:20:05 t_end=2024/01/09T10:10:05]

AU_Peg        21:24:00.24    +18:16:43.77     seq=1/Ic/30     priority=4     comment="adjust time to seeing"
BL_Her        18:01:09.22    +19:14:56.69     seq=10x(1/Ic/30)     priority=4        
BX_Del        20:21:18.97    +18:26:16.29     seq=1/Ic/30     priority=4     t_start=2024/01/09T1:20:05 t_end=2024/01/11           
SW_Tau        04:24:32.97     +04:07:24.08  seq=2/Ic/3.4,2/V/1.8   ph_mk=2/Ic/0.05    P=1.583535  priority=4

```
#### Specification:
* `seq=1/Ic/30,2/V/45` - (N/filter/exp_time) both sequences will be executed.  
  `seq=10x(1/Ic/30,2/V/45)` - (N/filter/exp_time) sequences will be repeated 10 times.   
  `seq=INFx(1/Ic/30,2/V/45)` - (N/filter/exp_time) sequences will be repeated as long it's visible.
* `priority=5`  - script allocates higher priority first. If priority is < -50, script will never allocate object
* `comment="some comment"` - comment visible to observer
*  `ph_mk=2/Ic/0.05` - parameters for phase masking:  no_of_obs/band_light_curve/phase_bin
* `t_start=2024/01/09T1:20:05` - date and ut  for start
  `t_start=2024/01/09` - alternative date only 
  `t_start=1:20:05` - alternative ut only (before 16:00 will be considered next day)
* `t_end=2024/01/09T10:10:05` - date and ut  for start  
  `t_end=2024/01/09` - alternative date only  
  `t_end=10:10:05` - alternative ut only (before 16:00 will be considered next day)
* `ph_star=0.1 ph_end=0.3 P=7.5 hjd0=2460310.2` - phase to observe, all keys are required   
  `ph_star=0.8 ph_end=0.1 P=7.5 hjd0=2460310.2` - yes, it is cyclic

### REMARKS
* no white spaces in *NAME* *RA* *DEC* *seq* 
* If there is no object to observe, it will allocate:
```WAIT ut=12:45:34``` with the `ut` for the next object
* comment has to be enclosed with **" "** characters 
* cycle and phase coverage are not implemented here, because observation database is not yet configured.   
* phase masking: object is masked if in defined `phase_bin` are already `no_of_obs` points for `band_light_curve`






