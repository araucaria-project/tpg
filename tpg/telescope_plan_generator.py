import sys,os,yaml,numpy,ephem
import random
import time

from astropy.table import Table

try:
    from ctc import CycleTimeCalc
    CTC_AVAILABLE = True
except ImportError:
    CTC_AVAILABLE = False


class TelescopePlanGenerator():
    def __init__(self,tel,date,tag=None,wind=None,loud=False,seed=None,add_start_makro=True,add_end_makro=True,done_uobi=[]):
        self.cwd = os.getcwd()
        self.script_dir = os.path.dirname(os.path.abspath(__file__))
        self.flags = []
        self.tag = tag
        self.wind = wind
        self.done = done_uobi
        self.loud = loud
        self.seed = seed
        self.add_start_makro = add_start_makro
        self.add_end_makro = add_end_makro

        self.msg = ""

        with open(self.script_dir+'/config.yaml', 'r') as cfg_file:
            self.cfg = yaml.safe_load(cfg_file)
            txt = f'config loaded from {self.script_dir}/config.yaml\n'
            if self.loud: sys.stderr.write(txt)
            self.msg = self.msg + txt

        local_config_file = self.cfg["local_config"]
        local_config_file = os.path.expanduser(local_config_file)
        if os.path.exists(local_config_file):
            with open(local_config_file) as f:
                local_cfg = yaml.safe_load(f) or {}
                self.cfg.update(local_cfg)
                txt = f'local config loaded from {local_config_file}\n'
                if self.loud: sys.stderr.write(txt)
                self.msg = self.msg + txt

        self.tel = tel
        if len(date) == 1:          # bo to lista
            self.flags.append("start_at_sunset")
            self.date = date[0] + " 12:00:00"
            txt = f'start at sunset\n'
        elif len(date) == 2:
            self.date = date[0] + " " + date[1]
            txt = f'start at {self.date}\n'
        else:
            print("date missing: yyyy/mm/dd")
            sys.exit(0)

        if "/" in date[0]:
            if len(date[0].split("/"))<3:
                print("wrong date format: yyyy/mm/dd")
                sys.exit(0)
            if len(date[0].split("/")[0])<4:
                print("wrong date format: yyyy/mm/dd")
                sys.exit(0)
        elif "-" in date[0]:
            if len(date[0].split("-")) < 3:
                print("wrong date format: yyyy/mm/dd")
                sys.exit(0)
            if len(date[0].split("-")[0])<4:
                print("wrong date format: yyyy/mm/dd")
                sys.exit(0)
        else:
            print("wrong date format: yyyy/mm/dd")
            sys.exit(0)

        self.plan = []

        if self.loud: sys.stderr.write(txt)
        self.msg = self.msg + txt

        # --------------- MAIN SCRIPT --------------
    def go(self):
        self.Initiate()  # initiate script variables
        self.LoadObjects()  # loads objects to observe
        self.MakeTime()  # prepares time variable
        self.CalcObject()  # calculates visibilities for all objects
        self.MaskMoon()
        self.MaskWind()
        self.MaskCycle()
        self.MaskStartEnd()  # mask visibilite base on t_start and t_end
        self.MaskPhaseStartEnd()  # mask visibilite base on t_start and t_end
        self.MaskPhase()
        self.Waga()  # calculate probability for randomization order of objects to observe
        self.RandomizeList()  # randomize listo of objects to observev
        self.allocate()  # SCHEDULING OBJECTS

        # wypisuje log
        vis_list = ["tag", "h_min", "h_max", "min_moon_dist", "max_moon_phase", "wind","cycle", "t_start", "t_end",
                    "ph_start","ph_end","ph_mk","all"]
        txt = "name    "+"    ".join(vis_list) + "\n"
        for n,ob in enumerate(self.ob):
            txt = txt + self.ob[n]["name"]
            txt = txt + f'  ({self.ob[n]["slotTime"]:.0f})m'
            for k in vis_list:
                if k in self.ob[n]["visibility"].keys():
                    vis = numpy.array(self.ob[n]["visibility"][k])
                    minutes = len(vis[vis])
                    txt = txt+ "    "+str(minutes)
                else:
                    txt = txt + "    --"

            txt = txt + "\n"

        if self.loud: sys.stderr.write(txt)
        self.msg = self.msg + txt

        self.export()  # export plan

    def Initiate(self):

        self.ctc_active = True

        # plan                      # 1 minute resolution
        self.nightTime = []  # time in ephem.Date format
        self.nightObject = []  # index of objects to observe
        self.nightFree = []  # Bool marker of allocation
        self.nightObChange = []    # Bool marker to indicate ob change

        self.ob = []

        # additional
        self.random_index_list = []  # list with random order of object indices to oberve
        self.pliki_hjd = []  # files with phase coverage for objects

    def LoadObjects(self):  # Load Object file
        # first try patch as it is, if file not found, take base path form __FILE__ and append config file path using os.path.join
        if not os.path.isfile(self.cfg[self.tel]["master_file"]):
            self.cfg[self.tel]["master_file"] = os.path.join(self.script_dir, self.cfg[self.tel]["master_file"])
            if not os.path.isfile(self.cfg[self.tel]["master_file"]):
                sys.stderr.write(f'No Objects Catalog file {self.cfg[self.tel]["master_file"]}\n')
                sys.stderr.write('Check config.yaml file\n')
                sys.exit()

        if os.path.exists(self.cfg["ctc"]["ctc_base_folder"]) and CTC_AVAILABLE:
            self.ctc = CycleTimeCalc(telescope=self.tel,base_folder=self.cfg["ctc"]["ctc_base_folder"], tpg=True)
            self.ctc.set_rm_modes(self.cfg["ctc"]["rm_modes_mhz"])
            self.ctc.set_start_rmode(self.cfg[self.tel]["def_read_mode"])
            self.ctc_active = True
        else:
            self.ctc_active = False

        txt = f'ctc active: {self.ctc_active}\n'
        if self.loud: sys.stderr.write(txt)
        self.msg = self.msg + txt

        file = open(self.cfg[self.tel]["master_file"], 'r')
        for line in file:
            if line.strip():
                if len(line.split()) > 0:  # no empty lines
                    if line.strip()[0] != "#":
                        tmp_ob = self.parseObjects(line)
                        self.ob.append(tmp_ob)


    def parseObjects(self,line):
        tmp_ob = dict()
        line = line.strip()
        if line:  # no empty lines

            if line.strip()[0] != "#":
                flag = []
                seq = "1/Ic/30"
                priority = 0
                comment = ""

                overhead = self.cfg[self.tel]["overhead"]
                tags = ""
                ts = False
                te = False
                hmax = False
                hmin = False
                ps = False
                pe = False
                period = False
                hjd0 = 2460000
                ph_mk = False
                cycle = False
                read_mod = None
                uobi = False
                makro = False
                obs_data = False
                max_moon_phase = False
                min_moon_dist = False



                if "comment=" in line:
                    tmp = line.split("comment=")[1]
                    if "\"" in line:
                        comment = tmp.split("\"")[1]
                    if "\'" in line:
                        comment = tmp.split("\'")[1]

                if "uobi=" in line:
                    uobi = line.split("uobi=")[1].split()[0]

                if "seq=" in line:
                    seq = line.split("seq=")[1].split()[0]

                if "priority=" in line:
                    priority = line.split("priority=")[1].split()[0]

                if "tag=" in line:
                    tags = line.split("tag=")[1].split()[0]

                if "t_start=" in line:
                    ts = line.split("t_start=")[1].split()[0]

                if "t_end=" in line:
                    te = line.split("t_end=")[1].split()[0]

                if "h_min=" in line:
                    hmin = line.split("h_min=")[1].split()[0]

                if "h_max=" in line:
                    hmax = line.split("h_max=")[1].split()[0]

                if "ph_start=" in line and "ph_end=" in line:
                    ps = line.split("ph_start=")[1].split()[0]
                    pe = line.split("ph_end=")[1].split()[0]

                if "P=" in line:
                    period = line.split("P=")[1].split()[0]

                if "hjd0=" in line:
                    hjd0 = line.split("hjd0=")[1].split()[0]

                if "ph_mk=" in line:
                    ph_mk = line.split("ph_mk=")[1].split()[0]

                if "cycle=" in line:
                    cycle = line.split("cycle=")[1].split()[0]

                if "read_mod=" in line:
                    read_mod = line.split("read_mod=")[1].split()[0]

                if "min_moon_dist" in line:
                    min_moon_dist = float(line.split("min_moon_dist=")[1].split()[0])

                if "max_moon_phase" in line:
                    max_moon_phase = float(line.split("max_moon_phase=")[1].split()[0])

                if "obs_data" in line:
                    obs_data = (line.split("obs_data=")[1].split()[0])

                if "makro" in line:
                    makro = (line.split("makro=")[1].split()[0])

                for i, tmp in enumerate(line.split()):

                    if i == 0:
                        name = line.split()[0]
                    elif i == 1:
                        ra = line.split()[1]
                    elif i == 2:
                        dec = line.split()[2]

                if "ob_time" in line:
                    ob_time = float(line.split("ob_time=")[1].split()[0])
                    slotTime = ob_time
                else:
                    slotTime = 10
                    k = 1
                    if "x(" in seq and ")" in seq:
                        k = seq.split("x")[0]
                        if k == "INF":
                            k = 1
                            flag = "inf"
                        k = int(k)

                        seq0 = seq.split("x(")[1].split(")")[0]
                    else:
                        seq0 = seq
                    for x_seq in seq0.split(","):
                        slotTime = slotTime + (
                                float(x_seq.split("/")[0]) * (float(x_seq.split("/")[2]) + float(overhead)))
                    slotTime = (k * slotTime) / 60.

                    if "INF" in seq:
                        tmp = seq.replace("INF", "1")
                    else:
                        tmp = seq

                    program = f"OBJECT  {name} {ra} {dec} seq={tmp}"
                    if read_mod is not None:
                        program = program + f" read_mod={read_mod}"

                    if self.ctc_active:
                        self.ctc.reset_time()
                        ctc_ob_time = self.ctc.calc_time(program) / 60.
                        slotTime = ctc_ob_time
                    else:
                        slotTime = calc_slot_time(program, overhead) / 60.


                if uobi and uobi in self.done:
                    pass
                else:

                    tmp_ob["name"] = name
                    tmp_ob["ra"] = ra
                    tmp_ob["dec"] = dec
                    tmp_ob["seq"] = seq
                    tmp_ob["comment"] = comment
                    tmp_ob["slotTime"] = slotTime
                    tmp_ob["line"] = line
                    tmp_ob["priority"] = priority
                    tmp_ob["tag"] = tags
                    tmp_ob["t_start"] = ts
                    tmp_ob["t_end"] = te
                    tmp_ob["h_min"] = hmin
                    tmp_ob["h_max"] = hmax
                    tmp_ob["ph_start"] = ps
                    tmp_ob["ph_end"] = pe
                    tmp_ob["period"] = period
                    tmp_ob["hjd0"] = hjd0
                    tmp_ob["ph_mk"] = ph_mk
                    tmp_ob["cycle"] = cycle
                    tmp_ob["uobi"] = uobi
                    tmp_ob["makro"] = makro
                    tmp_ob["obs_data"] = obs_data
                    tmp_ob["max_moon_phase"] = max_moon_phase
                    tmp_ob["min_moon_dist"] = min_moon_dist
                    tmp_ob["flag"] = flag

        return tmp_ob


    def MakeTime(self):  # Make observing time
        self.oca = ephem.Observer()
        self.oca.date = self.date
        self.oca.lon = str(-70.1963)
        self.oca.lat = str(-24.598)
        self.oca.elev = 2817
        self.oca.horizon = str(self.cfg[self.tel]["hsun"])  # altitude of sun to start and stop observations

        if "start_at_sunset" in self.flags:
            self.sunset = self.oca.next_setting(ephem.Sun(), use_center=True)

            txt = f'sunset: {self.sunset}\n'
            if self.loud: sys.stderr.write(txt)
            self.msg = self.msg + txt

            t = self.sunset
            self.oca.date = self.sunset
            self.sunrise = self.oca.next_rising(ephem.Sun(), use_center=True)

            txt = f'sunrise: {self.sunrise}\n'
            if self.loud: sys.stderr.write(txt)
            self.msg = self.msg + txt

        else:
            t = ephem.Date(self.date)
            self.oca.date = t
            self.sunrise = self.oca.next_rising(ephem.Sun(), use_center=True)

            txt = f'sunrise: {self.sunrise}\n'
            if self.loud: sys.stderr.write(txt)
            self.msg = self.msg + txt

        self.start_time = t
        while t < self.sunrise:
            t = t + 1 * ephem.minute  # time resolution : 1 minute
            self.nightTime.append(t)
            self.nightObject.append(-100)
            self.nightObChange.append(False)
            self.nightFree.append(True)


    def CalcObject(self):  # Colculate visibility for all objects
        for n, tmp in enumerate(self.ob):

            if self.ob[n].get("ra", False) and self.ob[n].get("dec", False) and self.ob[n].get("name", False):
                self.ob[n]["visibility"] = dict()

                # sprawdzanie opcji tag
                ok = True
                if self.tag:
                    ok = False
                    if d_ver(self.ob[n],["tag"]):
                        if self.tag in self.ob[n]["tag"]:
                            ok = True
                if ok:
                    vis = [True for _ in range(len(self.nightTime))]

                else:
                    vis = [False for _ in range(len(self.nightTime))]

                self.ob[n]["visibility"]["all"] = vis.copy()
                self.ob[n]["visibility"]["tag"] = vis.copy()

                # sprawdzanie widocznosci
                if all(self.ob[n]["visibility"]["all"]):
                    s = ephem.FixedBody()
                    s._ra = self.ob[n]["ra"]
                    s._dec = self.ob[n]["dec"]
                    self.oca.date = self.nightTime[0]

                    if self.ob[n].get("h_min",False):
                        self.oca.horizon = float(self.ob[n]["h_min"])
                    else:
                        self.oca.horizon = str(self.cfg[self.tel]["hmin"])
                    try:
                        tr = self.oca.next_rising(s)
                        ts = self.oca.next_setting(s)
                        if tr < ts:
                            vis1 = [x > tr and x < ts for x in self.nightTime]
                        else:
                            vis1 = [x < ts or x > tr for x in self.nightTime]
                    except ephem.NeverUpError:
                        vis1 = [False for _ in range(len(self.nightTime))]
                    except ephem.AlwaysUpError:
                        vis1 = [True for _ in range(len(self.nightTime))]

                    if self.ob[n].get("h_max",False):
                        self.oca.horizon = float(self.ob[n]["h_max"])
                    else:
                        self.oca.horizon = str(self.cfg[self.tel]["hmax"])

                    try:
                        tr = self.oca.next_rising(s)
                        ts = self.oca.next_setting(s)
                        if tr < ts:
                            vis2 = [not (x > tr and x < ts) for x in self.nightTime]
                        else:
                            vis2 = [not (x < ts or x > tr) for x in self.nightTime]
                    except ephem.NeverUpError:
                        vis2 = [True for _ in range(len(self.nightTime))]
                    except ephem.AlwaysUpError:
                        vis2 = [False for _ in range(len(self.nightTime))]

                    self.ob[n]["visibility"]["h_min"] = vis1
                    vis = [a and b for a,b in zip(vis1,vis2)]
                    self.ob[n]["visibility"]["h_max"] = vis.copy()
                    self.ob[n]["visibility"]["all"] = vis

    def MaskMoon(self):
        moon = ephem.Moon()
        self.oca.date = self.nightTime[0]
        self.oca.horizon = 0
        tr = self.oca.next_rising(moon)
        ts = self.oca.next_setting(moon)
        if tr < ts:
            moonvis = [not (x > tr and x < ts) for x in self.nightTime]
        else:
            moonvis = [(x > tr or x < ts) for x in self.nightTime]

        for n, tmp in enumerate(self.ob):
            if d_ver(self.ob[n],["name","ra","dec"]):
                s = ephem.FixedBody()
                s._ra = self.ob[n]["ra"]
                s._dec = self.ob[n]["dec"]
                tmp = numpy.array(self.nightTime)
                time = tmp[numpy.array(self.ob[n]["visibility"]["all"])]
                if len(time)>0:
                    self.oca.date = time[0]
                    moon.compute(self.oca)
                    s.compute(self.oca)
                    phase = float(moon.phase)
                    sep = ephem.separation(s,moon) * 180/ephem.pi

                    if d_ver(self.ob[n], ["min_moon_dist"]):
                        moon_dist = float(self.ob[n]["min_moon_dist"])
                    else:
                        moon_dist = float(self.cfg[self.tel]["moon_dist"])

                    if float(sep) < float(moon_dist):
                        # mask everything for faster calculation
                        self.ob[n]["visibility"]["all"] = [False for _ in range(len(self.nightTime))]
                        self.ob[n]["visibility"]["min_moon_dist"] = [False for _ in range(len(self.nightTime))]

                    if d_ver(self.ob[n], ["max_moon_phase"]):
                        if float(phase) > float(self.ob[n]["max_moon_phase"]):
                            self.ob[n]["visibility"]["all"] = moonvis
                            self.ob[n]["visibility"]["max_moon_phase"] = moonvis

    def MaskWind(self):  # cala noc, jeden kierunek
        if self.wind:
            for n, tmp in enumerate(self.ob):
                if d_ver(self.ob[n], ["name", "ra", "dec","visibility"]):
                    self.ob[n]["visibility"]["wind"] = self.ob[n]["visibility"]["all"].copy()
                    s = ephem.FixedBody()
                    s._ra = self.ob[n]["ra"]
                    s._dec = self.ob[n]["dec"]

                    time = numpy.array(self.nightTime)

                    for t,tmp in enumerate(time):
                        if self.ob[n]["visibility"]["all"][t]:
                            self.oca.date = time[t]
                            s.compute(self.oca)
                            az = s.az * 180/ephem.pi
                            dist = numpy.abs(self.wind - az)
                            if dist > 180:
                                dist = 360. - dist
                            if dist < float(self.cfg[self.tel]["wind_dist"]):
                                self.ob[n]["visibility"]["all"][t] = False
                                self.ob[n]["visibility"]["wind"][t] = False

    def MaskCycle(self):  # mask cycle covarage
        for n, tmp in enumerate(self.ob):
            if d_ver(self.ob[n], ["name","cycle", "visibility"]):
                ZNACZNIK = False
                if self.ob[n]["cycle"]:
                    jd = []
                    try:
                        if d_ver(self.ob[n], ["obs_data"]):
                            file = self.ob[n]["obs_data"]
                            txt = f'{self.ob[n]["name"]} data source: {file}\n'
                            if self.loud: sys.stderr.write(txt)
                            self.msg = self.msg + txt

                        else:
                            filter = self.ob[n]["seq"].split("/")[1]
                            file = self.cfg[self.tel]["data_file"] + self.ob[n]["name"].lower() + "/" + filter + "/light-curve/" + self.ob[n]["name"].lower() + "_" + filter + "_diff_light_curve.txt"

                        lc_tab = Table.read(file, format="ascii")
                        jd = lc_tab["jd_obs"]



                    except (FileNotFoundError) as e:
                        pass
                    except (ValueError) as e:
                        txt = f'ERROR loading {file}\n'
                        if self.loud: sys.stderr.write(txt)
                        self.msg = self.msg + txt


                    last_jd = False
                    try:
                        last_jd = max(jd)
                    except:
                        pass
                    if last_jd:
                        vis_tmp = numpy.array(self.ob[n]["visibility"]["all"])
                        for i, t in enumerate(self.nightTime):
                            if self.ob[n]["visibility"]["all"][i]:
                                jd = ephem.julian_date(ephem.Date(t))
                                if (jd - last_jd) < float(self.ob[n]["cycle"]):
                                    self.ob[n]["visibility"]["all"][i] = False
                                    ZNACZNIK = True
                                else:
                                    pass

                if ZNACZNIK:
                    vis = self.ob[n]["visibility"]["all"]
                    self.ob[n]['visibility']['cycle'] = vis
                    vis = numpy.array(vis)
                    t1 = len(vis_tmp[vis_tmp])
                    t2 = len(vis[vis])

    def MaskStartEnd(self):  # Masks visibility based on t_start t_end
        for n, tmp in enumerate(self.ob):
            ZNACZNIK = False
            if d_ver(self.ob[n], ["name", "t_start", "visibility"]):
                ts = self.ob[n]["t_start"]
                if ts:
                    if "T" in ts:
                        ts = ts.replace("T", " ")
                    if "/" not in ts:
                        date = str(self.start_time).split()[0]

                        if int(str(self.start_time).split()[1].split(":")[0]) > 16 and int(ts.split(":")[0]) < 16:
                            date = ephem.Date(date) + 1
                            date = str(ephem.Date(date)).split()[0]
                        ts = date + " " + ts
                    ts = ephem.Date(ts)
                    time_s = ts

                    if time_s:
                        for i, tmp in enumerate(self.nightTime):
                            if float(self.nightTime[i]) < float(time_s) and self.ob[n]["visibility"]["all"][i] and time_s:
                                self.ob[n]["visibility"]["all"][i] = False
                                ZNACZNIK = True
                self.ob[n]["visibility"]["t_start"] = self.ob[n]["visibility"]["all"].copy()

            if d_ver(self.ob[n], ["name", "t_end", "visibility"]):
                te = self.ob[n]["t_end"]
                if te:
                    if "T" in te:
                        te = te.replace("T", " ")
                    if "/" not in te:
                        date = str(self.start_time).split()[0]

                        if int(str(self.start_time).split()[1].split(":")[0]) > 16 and int(te.split(":")[0]) < 16:
                            date = ephem.Date(date) + 1
                            date = str(ephem.Date(date)).split()[0]
                        te = date + " " + te
                    te = ephem.Date(te)
                    time_e = te

                    if time_e:
                        for i, tmp in enumerate(self.nightTime):
                            if float(self.nightTime[i]) > float(time_e) and self.ob[n]["visibility"]["all"][i] and time_e:
                                self.ob[n]["visibility"]["all"][i] = False
                                ZNACZNIK = True
                self.ob[n]["visibility"]["t_end"] = self.ob[n]["visibility"]["all"].copy()
                print(ts,te)


    def MaskPhaseStartEnd(self):  # Masks visibility based on pstart and pend
        for n, tmp in enumerate(self.ob):
            ZNACZNIK = False
            if d_ver(self.ob[n], ["name", "ph_start","ph_end", "visibility"]):
                ps = self.ob[n]["ph_start"]
                pe = self.ob[n]["ph_end"]
                period = self.ob[n]["period"]
                hjd0 = self.ob[n]["hjd0"]

                if True:
                    for i, t in enumerate(self.nightTime):
                        if self.ob[n]["visibility"]["all"][i]:
                            jd = ephem.julian_date(ephem.Date(t))
                            p = (jd - float(hjd0)) / float(period) % 1
                            if pe < ps:
                                if p > float(pe) and p < float(ps):
                                    self.ob[n]["visibility"]["all"][i] = False
                                    ZNACZNIK = True
                                else:
                                    pass
                            else:
                                if p > float(ps) and p < float(pe):
                                    pass
                                else:
                                    self.ob[n]["visibility"]["all"][i] = False
                                    ZNACZNIK = True

                    if ZNACZNIK:
                        self.ob[n]["visibility"]["ph_start"] = self.ob[n]["visibility"]["all"].copy()
                        self.ob[n]["visibility"]["ph_end"] = self.ob[n]["visibility"]["all"].copy()


    def MaskPhase(self):  # mask phase covarage
        for n, tmp in enumerate(self.ob):
            if d_ver(self.ob[n], ["name", "ph_mk", "visibility"]):
                ZNACZNIK = False
                if self.ob[n]["ph_mk"]:
                    jd = []
                    try:
                        filter = self.ob[n]["ph_mk"].split("/")[1]
                        if d_ver(self.ob[n], ["obs_data"]):
                            file = self.ob[n]["obs_data"]
                        else:
                            file = self.cfg[self.tel]["data_file"] + self.ob[n]["name"].lower() + "/" + filter + "/light-curve/" + self.ob[n]["name"].lower() + "_" + filter + "_diff_light_curve.txt"

                        lc_tab = Table.read(file, format="ascii")
                        jd = lc_tab["jd_obs"]

                    except (FileNotFoundError) as e:
                        pass
                    except (ValueError) as e:
                        txt = f'ERROR loading {file}\n'
                        if self.loud: sys.stderr.write(txt)
                        self.msg = self.msg + txt


                    P = self.ob[n]["period"]
                    hjd0 = self.ob[n]["hjd0"]
                    if P and len(jd)>0:
                        ph = (numpy.array(jd) - float(hjd0)) / float(P) % 1

                        night_ph = (numpy.array(self.nightTime) - float(hjd0)) / float(P) % 1

                        for i,tmp in enumerate(night_ph):
                            np1 = night_ph[i] - float(self.ob[n]["ph_mk"].split("/")[2])
                            np2 = night_ph[i] + float(self.ob[n]["ph_mk"].split("/")[2])
                            maska1 = ph > np1
                            maska2 = ph < np2
                            maska = [a and b for a,b in zip(maska2,maska1)]
                            if len(ph[maska]) > int(self.ob[n]["ph_mk"].split("/")[0]):
                                self.ob[n]["visibility"]["all"][i] = False
                                ZNACZNIK = True

                if ZNACZNIK:
                    self.ob[n]["visibility"]["ph_mk"] = self.ob[n]["visibility"]["all"].copy()


    def Waga(self):  #
        no_obj = False
        if True:
            for n, tmp in enumerate(self.ob):
                if d_ver(self.ob[n], ["name", "visibility"]):
                    maska = numpy.array(self.ob[n]["visibility"]["all"])
                    maska = maska[maska]
                    allminutes = len(self.nightTime)
                    if len(maska) == 0:
                        waga = 0
                    else:
                        no_obj = True
                        waga = 1. - float(len(maska)) / float(allminutes)
                    if waga < 0.3 and len(maska) != 0: waga = 0.3
                    self.ob[n]["waga"] = waga
        if not no_obj:
            raise ValueError(f"No observable objects!")


    # tutaj mozna dodac opcje wylaczenia losowania
    def RandomizeList(self):  # randomize listo of objects to observev
        mozliwosci = numpy.array([i for i, _ in enumerate(self.ob)])
        probability = numpy.array([ob["waga"] for ob in self.ob])
        maska = (probability != 0)
        mozliwosci = mozliwosci[maska]
        probability = probability[maska]

        suma = probability.sum()
        if suma > 0:
            probability = probability / suma

        if not self.seed:
            self.seed = random.randrange(0, 2 ** 32)

        txt = f"random seed: {self.seed}\n"
        if self.loud: sys.stderr.write(txt)
        self.msg = self.msg + txt

        numpy.random.seed(self.seed)
        self.random_index_list = numpy.random.choice(mozliwosci, len(mozliwosci), replace=False, p=probability)

    def allocate(self):
        for nn, tmp in enumerate(self.nightTime):
            if self.nightFree[nn]:

                OK = False
                done = -1
                priormax = -1

                for nk, tmp in enumerate(self.random_index_list):
                    k = self.random_index_list[nk]

                    # wbrew pozorom ta funkcja przeszukuje caly zbior, ale wybiera pierwsz obiekt
                    # z listy random_index_list ( bo jest > priormax, a nie >= )
                    # moze to nie jest najbardziej optymalny kod
                    if int(self.ob[k]["priority"]) > priormax and int(self.ob[k]["priority"]) > -1:
                        visibility = self.ob[k]["visibility"]["all"]
                        slots = float((self.ob[k]["slotTime"])) / 1.
                        if "inf" in self.ob[k]["flag"]:
                            loops = 1
                            go = True
                            while go:
                                if all(visibility[nn:(nn + int(loops * slots)) + 2]) and len(visibility[nn:]) > int(
                                        loops * slots) + 1:
                                    loops = loops + 1
                                else:
                                    go = False
                            loops = loops
                            slots = int(slots * loops)
                        slots = int(slots)
                        bounds = (nn, nn + slots)

                        if all(self.nightFree[bounds[0]:bounds[1] + 1]) and all(visibility[bounds[0]:bounds[1] + 1]):
                            priormax = int(self.ob[k]["priority"])
                            bfin = bounds
                            kfin = k
                            OK = True
                            done = nk


                if OK:
                    self.nightFree[bfin[0]:bfin[1] + 1] = [False] * ((bfin[1] + 1) - bfin[0])
                    self.nightObject[bfin[0]:bfin[1] + 1] = [int(kfin)] * ((bfin[1] + 1) - bfin[0])
                    self.nightObChange[bfin[0]] = True

                    if "inf" in self.ob[kfin]["flag"]:
                        self.ob[kfin]["seq"] = self.ob[kfin]["seq"].replace("INFx", f"{loops}x")

                    # usuwanie juz zaalocowanych
                    q = self.random_index_list[done]

                    if d_ver(self.ob[q], ["uobi"]):
                        uobi_tmp = self.ob[q]["uobi"]
                        if len(uobi_tmp.strip()) > 0:
                            indeksy = [i for i, item in enumerate(self.ob) if item["uobi"] == uobi_tmp]
                            ni = [j for j, val in enumerate(self.random_index_list) if val in indeksy]
                    else:
                        ni = [done]

                    del_list = []
                    for ki in ni:
                        q = self.random_index_list[ki]
                        if self.ob[q]["cycle"] != False:
                            dn = int(1440 * float(self.ob[q]["cycle"]))  # 1440 minut ma doba
                            if dn > len(self.ob[q]["visibility"]["all"][nn:]):
                                dn = len(self.ob[q]["visibility"]["all"][nn:])
                            self.ob[q]["visibility"]["all"][nn:nn+dn] = [False] * dn
                        else:
                            del_list.append(ki)

                    self.random_index_list = numpy.delete(self.random_index_list, del_list)



    def export(self):
        if True:
            ZNACZNIK = False
            if "start_at_sunset" in self.flags:
                txt = f'WAIT sunset={str(self.cfg[self.tel]["hsun"])}'
            else:
                txt = f'WAIT ut={str(ephem.Date(self.nightTime[0])).split()[1]}'
            self.plan.append(txt)
            for n, tmp in enumerate(self.nightTime):
                if self.nightFree[n]:
                    ZNACZNIK = True
                if self.nightObChange[n]:  # tutaj sprawdzamy ze jest nowy obiekt alokowany
                    if ZNACZNIK:   # tutaj wstawiamy wait, jak nie byla przerwa w alokacji
                        ZNACZNIK = False
                        self.oca.date = ephem.Date(self.nightTime[n])
                        txt = "WAIT ut=" + str(self.oca.date).split()[1]
                        self.plan.append(txt)
                    if True:
                        k = int(self.nightObject[n])
                        if d_ver(self.ob[k], ["name", "makro"]):
                                with open(self.ob[k]["makro"], "r", encoding="utf-8") as makro_file:
                                    for line in makro_file:
                                        if line.strip():
                                            self.plan.append(line.strip())
                        else:
                            txt = self.ob[k]["line"]
                            par2remove = self.cfg["params_2_remove"]
                            if d_ver(self.ob[k], ["uobi","cycle"]):
                                par2remove.append("uobi")  # bo nie powinno sie powtarzac uobi
                            for x in par2remove:
                                if x+"=" in txt:
                                    ab = txt.split(x+"=")
                                    a = ab[0]
                                    b = ab[1]
                                    r = b.split()
                                    w = " ".join(r[1:])
                                    q = a.strip() + " " + w.strip()
                                    txt = "   ".join(q.split())
                            self.plan.append(f'OBJECT {txt}')


        if self.add_start_makro and d_ver(self.cfg[self.tel],["start_makro"]):
            start_plan = []
            ok = False
            if os.path.exists(self.cfg[self.tel]["start_makro"]):
                with open(self.cfg[self.tel]["start_makro"], "r", encoding="utf-8") as makro_file:
                    for line in makro_file:
                        if line.strip():
                            start_plan.append(line.strip())
                            ok = True
            if ok:
                self.plan = start_plan + self.plan
                txt = f"start_makro added\n"
                if self.loud: sys.stderr.write(txt)
                self.msg = self.msg + txt
            else:
                txt = f"NO start_makro FILE\n"
                if self.loud: sys.stderr.write(txt)
                self.msg = self.msg + txt

        if self.add_end_makro and d_ver(self.cfg[self.tel],["end_makro"]):
            end_plan = []
            ok = False
            if os.path.exists(self.cfg[self.tel]["end_makro"]):
                with open(self.cfg[self.tel]["end_makro"], "r", encoding="utf-8") as makro_file:
                    for line in makro_file:
                        if line.strip():
                            end_plan.append(line.strip())
                            ok = True
            if ok:
                self.plan = self.plan + end_plan
                txt = f"end_makro added\n"
                if self.loud: sys.stderr.write(txt)
                self.msg = self.msg + txt
            else:
                txt = f"NO end_makro FILE\n"
                if self.loud: sys.stderr.write(txt)
                self.msg = self.msg + txt


def d_ver(slownik,klucze):
    if not isinstance(klucze, (list, tuple)):
        raise TypeError(f'klucze should be a list or tuple, not {type(klucze)}')
    result = True
    for k in klucze:
        if k not in slownik.keys():
            result = False
        else:
            if slownik[k]:
                pass
            else:
                result = False
    return result


def calc_slot_time(program, overhed):
    slotTime = 10
    k = 1
    seq = program.split("seq=")[1].split()[0]
    try:
        if "x(" in seq and ")" in seq:
            k = int(seq.split("x")[0])
            seq = seq.split("x(")[1].split(")")[0]
        for x_seq in seq.split(","):
            if "a" in x_seq:
                x_seq = x_seq.replace("/a", "/5")
            slotTime = slotTime + (float(x_seq.split("/")[0]) * (float(x_seq.split("/")[2]) + float(overhed)))
        slotTime = k * slotTime
    except:
        slotTime = 0
    return slotTime
