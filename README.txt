TELESCOPE PLAN GENERATOR README
08.08.2023
Marek Gorski

prepares observing plans for zb08 and jk15 telescopes.

Running script:

tpg.py telescope (zb08 or jk15)  2023/08/10  >| zb08_230810.txt


file paths are specified in the config file

CONFIGURATION file : scheduler.cfg

zb08_FILE=./zb08_master.txt

# **NOT IMPLEMENTED** OBSDATA=/home/observer/scripts/OCANV/ObsData            # file with observed jd of objects
# **NOT IMPLEMENTED** TABALL=/home/observer/scripts/OCANV/CATALOGS/TAB.ALL    # file with periods for phase calculations
# **NOT IMPLEMENTED** LOGDIR=/home/observer/Documents/v16/                    # file with log of observed objects for cycle
      
# **NOT IMPLEMENTED**  OBSNUM=3               # minimum total number of observations to reject
# **NOT IMPLEMENTED** PHCOV=0.05             # phase distance to reject
# **NOT IMPLEMENTED** PHNUM=2                # minimum number of observations in phase to reject

HMIN=40                # default visibility min altitude condition for all objects
HMAX=85                # default visibility max altitude condition for all objects
HSUN=-9                # default sun altitude, as the begining of the obserwations
      

Script alocate time based on priority and visibility. Calculates slot time based on specification ''seq'' each object

# NAME RA DEC seq=1/Ic/30  priority=10          # **NOT IMPLEMENTED** nv_overhead=50  nv_total_time=800   flag=ignore_phs flag=ignore_vis

Additional markers are:

   seq=1/Ic/30                 -
   priority=5                  - script allocates higher priority first. If priority is < -50, script will never allocate object
   comment="some comment"
# **NOT IMPLEMENTED**   nv_overhead=50              - overhead in seconds added to each DIT
# **NOT IMPLEMENTED**   nv_total_time=800           - totl time in seconds (not calculating DIT and NDIT)
# **NOT IMPLEMENTED**   t_start=22:30:00            - visibility only after preset time each night
# **NOT IMPLEMENTED**   t_start=2019/11/21d22:30:00  - visibility only after preset day time
# **NOT IMPLEMENTED**   t_end=05:00:00              - visibility only before preset time each night
# **NOT IMPLEMENTED**   t_end=2019/11/21d05:00:00    - visibility only after preset day time
# **NOT IMPLEMENTED**   flag=ignore_phs             - with --ph option, will ignore phase coverage. Object will be observed, even if it's phase is already covered
# **NOT IMPLEMENTED**   flag=ignore_vis             - will ignore visibility. Exapmple - objects lower than 40 deg, or FLAT for IRIS.
# **NOT IMPLEMENTED**   cycle=4                     - change priority based on log files and cycle


Script reads each object from master file for specific telescope, calculate visibility and allocate objects.

If there is no object to observe, it will allocate:

WAIT ut=(time for next object)

if priority < -50 : no allocation for object

# **NOT IMPLEMENTED** phase coverage : object will be rejected to observe if there is more than 4 observation already done,
# **NOT IMPLEMENTED** and in range of 0.1 phase there is at least 2 observations done


# **NOT IMPLEMENTED** cycle: if object was obserwed (scrit looks in log files) within last cycle days, priority will be changed to 1
# **NOT IMPLEMENTED** if object was not observed since for more than cycle: priority=priority+1
