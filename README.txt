SCHEDULER README
06.12.2019
Marek Gorski

prepares observing plans for V16 and IRIS telescopes.


Running script:

python27 scheduler.py telescope (v16 or iris)  2019/08/20 [--nv --ph --loud]  >| v16_plan_nv.txt

--nv : print plan in the NV file format; readable with NV
--ph : reject object based on phase caverage"
--loud : write to std err additional informations 


script uses files in the format of NV catalogs. 
file paths are specified as the variables in the config file 

CONFIGURATION file : scheduler.cfg

IRISFILE=/Users/mgorski/DEVELOPER/OCA/OCANV/CATALOGS/iris_master_nv.txt
V16FILE=/Users/mgorski/DEVELOPER/OCA/OCANV/CATALOGS/v16_master_nv.txt

OBSDATA=/home/observer/scripts/OCANV/ObsData            # file with observed jd of objects
TABALL=/home/observer/scripts/OCANV/CATALOGS/TAB.ALL    # file with periods for phase calculations
LOGDIR=/home/observer/Documents/v16/                    # file with log of observed objects for cycle 
      
OBSNUM=3               # minimum total number of observations to reject
PHCOV=0.05             # phase distance to reject
PHNUM=2                # minimum number of observations in phase to reject

HMIN=40                # default visibility altitude condition for all objects
HSUN=-9                # default sun altitude, as the begining of the obserwations
      

Script alocate time based on priority and visibility. Calculates slot time based on DIT and NDIT for each object 

# NAME RA DEC BAND DIT NDIT   nv_overhead=50  nv_total_time=800  priority=10 flag=ignore_phs flag=ignore_vis

Additional markers are:

   priority=5                  - program allocastes higher priority first. If priority is < -50, script will never allocate object
   nv_overhead=50              - overhead in seconds added to each DIT
   nv_total_time=800           - totl time in seconds (not calculating DIT and NDIT)
   t_start=22:30:00            - visibility only after preset time each night
   t_start=2019/11/21d22:30:00  - visibility only after preset day time
   t_end=05:00:00              - visibility only before preset time each night
   t_end=2019/11/21d05:00:00    - visibility only after preset day time
   flag=ignore_phs             - with --ph option, will ignore phase coverage. Object will be observed, even if it's phase is already covered 
   flag=ignore_vis             - will ignore visibility. Exapmple - objects lower than 40 deg, or FLAT for IRIS. 
   cycle=4                     - change priority based on log files and cycle 


Script reads ecah object from files NV catalogs, calculate visibility, randomize list, mask phase and cycle and allocate objects. 

If there is no object to observe, it will allocate object 

WASTING_TIME object with ra=sideral time and dec=-10 deg

if priority < -50 : no allocation for object

phase coverage : object will be rejected to observe if there is more than 4 observation already done, 
and in range of 0.1 phase there is at least 2 observations done


cycle: if object was obserwed (scrit looks in log files) within last cycle days, priority will be changed to 1
if object was not observed since for more than cycle: priority=priority+1
