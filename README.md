# tpg - Telescope Plan Generator

Script for generating observing plans for **zb08** and **jk15** telescopes. Script alocates time based on priority and visibility. Calculates slot time based on specification ``seq`` of each object. Creates files readable by *TOI* and *TEXTUI*.

### RUNNING SCRIPT
```
tpg.py <zb08|jk15> <date> (YYYY/MM/DD) [optional: master_catalog_file]
```
```
tpg.py zb08 2023/08/1  >| zb08_230810.txt
```

### INSTALLATION
pull from github
```commandline
cd tpg
pipx install .
```
upgrading
```commandline
pipx upgrade tpg
```

### CONFIGURATION
`tpg.cfg`

```plainttext
zb08_FILE=~/tpg/zb08_master.txt    # file with objects with specification for  
HMIN=35                        # default visibility based on min altitude condition for all objects
HMIN=80                        # default visibility based on max altitude condition for all objects
HSUN=-18                        # default sun altitude, as the begining of the obserwations
```      
### MASTER FILE (format)
`zb08_master.txt`
```
# NAME RA DEC seq=<N/filter/exp_time>  priority=<priority> [optional: comment="<comment for observer>" t_start=2024/01/09T1:20:05 t_end=2024/01/09T10:10:05]

AU_Peg	    21:24:00.24    +18:16:43.77 	seq=1/Ic/30     priority=4	 comment="adjust time to seeing"
BL_Her	    18:01:09.22    +19:14:56.69 	seq=10x(1/Ic/30)     priority=4	    
BX_Del	    20:21:18.97    +18:26:16.29 	seq=1/Ic/30     priority=4	 t_start=2024/01/09T1:20:05 t_end=2024/01/11           
```
#### Specification:
* `seq=1/Ic/30,2/V/45` - (N/filter/exp_time) both sequences will be executed.  
  `seq=10x(1/Ic/30,2/V/45)` - (N/filter/exp_time) sequences will be repeated 10 times.
* `priority=5`  - script allocates higher priority first. If priority is < -50, script will never allocate object
* `comment="some comment"` - comment visible to observer
* `t_start=2024/01/09T1:20:05` - date and ut  for start
  `t_start=2024/01/09` - alternative date only 
  `t_start=1:20:05` - alternative ut only (before 16:00 will be considered next day)
* `t_end=2024/01/09T10:10:05` - date and ut  for start  
  `t_end=2024/01/09` - alternative date only  
  `t_end=10:10:05` - alternative ut only (before 16:00 will be considered next day)

### REMARKS
* no white spaces in *NAME* *RA* *DEC* *seq* 
* If there is no object to observe, it will allocate:
```WAIT ut=12:45:34``` with the `ut` for the next object
* comment has to be enclosed with **" "** characters 
* cycle and phase coverage are not implemented here, because observation database is not yet configured.   







