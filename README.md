# tpg - Telescope Plan Generator

Script for generating observing plans for **zb08** and **jk15** telescopes. Script alocates time based on priority and visibility. Calculates slot time based on specification ``seq`` of each object. Creates files readable by *TOI* and *TEXTUI*.

### RUNNING SCRIPT:
```
tpg.py <zb08|jk15> <date> (YYYY/MM/DD)
```
```
tpg.py zb08 2023/08/1  >| zb08_230810.txt
```


### CONFIGURATION
`scheduler.cfg`

```plainttext
zb08_FILE=./zb08_master.txt    # file with objects with specification for  
HMIN=40                        # default visibility altitude condition for all objects
HSUN=-9                        # default sun altitude, as the begining of the obserwations
```      
### MASTER FILE (format)
`./zb08_master.cfg`
```
# NAME RA DEC seq=<N/filter/exp_time>  priority=<priority>

AU_Peg	    21:24:00.24    +18:16:43.77 	seq=1/Ic/30     priority=4	 #binary and dusty envelope
BL_Her	    18:01:09.22    +19:14:56.69 	seq=1/Ic/30     priority=4	    
BX_Del	    20:21:18.97    +18:26:16.29 	seq=1/Ic/30     priority=4	            
```
#### Specification:
* `seq=1/Ic/30,2/V/45` - (N/filter/exp_time) both sequences will be executed.            
* `priority=5`  - script allocates higher priority first. If priority is < -50, script will never allocate object

### REMARKS
* no white spaces in *specification*
* If there is no object to observe, it will allocate:
```WAIT ut=12:45:34``` with the `ut` for the next object
* cycle and phase coverage are not implemented here, because observation database is not yet configured.   







