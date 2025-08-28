# tpg - Telescope Plan Generator

Script for generating observing plans for **wk06**, **zb08** and **jk15** telescopes. Script alocates time based on priority and visibility. Calculates slot time based on specification ``seq`` of each object. Creates files readable by *TOI* and *TEXTUI*.

### RUNNING SCRIPT
```
poetry run tpg {wk06,zb08,jk15} date [ut] [--wind WIND] [-h] [--loud]
date: YYYY/MM/DD [hh:mm:ss] 
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
example20   07:24:14         -00:31:41     makro=./test/makro_ex10.txt   priority=10 cycle=5    ob_time=300 obs_data=./test/obs_data/example5.txt


```
#### Specification:
* `seq=1/Ic/30,2/V/45` - (N/filter/exp_time) both sequences will be executed.  
  `seq=10x(1/Ic/30,2/V/45)` - (N/filter/exp_time) sequences will be repeated 10 times.   
  `seq=INFx(1/Ic/30,2/V/45)` - (N/filter/exp_time) sequences will be repeated as long it's visible.
* `priority=5`  - script allocates higher priority first. If priority is < -50, script will never allocate object
* `comment="some comment"` - comment visible to observer
* `cycle=3` - object will be observed if last observation was more than 3 days ago.

  `cycle=0.3` - object can be repeted during night 0.3*24h after last planed scheduling.
*  `ph_mk=2/Ic/0.05` - parameters for phase masking:  no_of_obs/band_light_curve/phase_bin
*  `obs_data=./test/obs_data/example5.txt` - data file if cycle or phase data file is different then object name
* `t_start=2024/01/09T1:20:05` - date and ut  for start

  `t_start=2024/01/09` - alternative date only 

  `t_start=1:20:05` - alternative ut only (before 16:00 will be considered next day)
* `t_end=2024/01/09T10:10:05` - date and ut  for start  
  `t_end=2024/01/09` - alternative date only  
  `t_end=10:10:05` - alternative ut only (before 16:00 will be considered next day)
* `ph_star=0.1 ph_end=0.3 P=7.5 hjd0=2460310.2` - phase to observe, all keys are required   
  `ph_star=0.8 ph_end=0.1 P=7.5 hjd0=2460310.2` - yes, it is cyclic
* `h_min=20` - object can be observed if more than 20 deg above horizon
* `h_max=70` - object can be observed if less than 70 deg above horizon
* `min_moon_dist=30` - object can be observed if more than 30 deg from the moon
* `max_moon_phase=0.7` - object can be observed if moon phase is less than 0.7 
* `ob_time=120` - slot time will be set as 120 seconds (not calculated by tpg)
* `makro=./test/makro_ex10.txt` - if scheduled file makro_ex10.txt will be included in the plan, instead of ob line
  `makro=./test/makro_ex10.txt priority=10 cycle=5 ob_time=300 obs_data=./test/obs_data/example5.txt` - should work


### REMARKS
* no white spaces in *NAME* *RA* *DEC* *seq* 
* If there is no object to observe, it will allocate:
```WAIT ut=12:45:34``` with the `ut` for the next object
* comment has to be enclosed with **" "** characters 
* phase masking: object is masked if in defined `phase_bin` are already `no_of_obs` points for `band_light_curve`
* for `cycle` full funcionality, obs data are required, and for this the filter has to be specified. If `ph_mk` is not specified, it will take first filter defined in `seq`.

### How it works
Each object specified in the master file has multiple atributes. Objects are divided in to different *priority* groups, and each object has specific visibility and scheduling conditions. For each object for specified night and geographical location it's altitude is calculated for each minute of specified night. Based on this array, the *visibility* array is created, where *True* value is assigned if object altitude is above *hmin* and below *hmax* atributes. This base visibility array then is masked with additinal atributes, like *t_star*, *t_end*, *ph_start*, *ph_end* and *ph_mk*.  

Additionaly, *visibility* can be masked for wind direction (based on calculated azimuth for each night minute), and distance from the moon.
The reuslt is the finall *visibility* array.

In the next step, for each object the *predicted execution time* is calculated with *ctc*, based on specified *sequence*. 
Object allocation is sequential, that is for each minute of the night, starting from specified sun altitude, one object is alocated for observing slot equal to its *predicted execution time*. Next object is allocated afterwards, repeatong untill specified end of the night - sun altitude for sunrise. Objects alocation works like this: From all objects, objects with highest priorities are selected, and one of them is selected randomly if it's visibility conditions are true for all it's *predicted execution time*. Random selection is weighted with it's total voisibility time for particular night - there is higer chance for selecting object with shorter visibility. If object is allocated, it's removced from the list of potential objects to allocate, and the time marker is moved to the end of it's observing slot (*predicted execution time*). Hovewer, if for this particular object, *cycle* atribute is specified, than instead of removing the object from the potential objects for allocation, it's *visibility* array is masked *False* from actual time marker (allocation time), untill *cycle* part of 24h. If no object of particular *priority* is avialable for the slot, the priority is lowered and precess is repeted for objects with lower priority. 

#### Simplified schema:
-> Loading objects from master file 


-> visibility is calclated for all nigh, based on it's alititude 


-> visibility is masked based on additional atributes 


-> objects are allocated, starting from the highest priority, randomizing selection within priority group. 


-> time marker is shifted and object for next slot is selected. 





