import sys,os,yaml,numpy,ephem
import time

from ctc import CycleTimeCalc


class telescope_plan_generator():
    def __init__(self,tel,date,wind=None,loud=False,done_uobi=[]):

        self.cwd = os.getcwd()
        self.script_dir = os.path.dirname(os.path.abspath(__file__))
        self.flags = []
        self.done = done_uobi
        self.loud = loud
        self.msg = ""

        with open(self.script_dir+'/config.yaml', 'r') as cfg_file:
            self.cfg = yaml.safe_load(cfg_file)
            txt = f'config loaded from config.yaml\n'
            if self.loud: sys.stderr.write(txt)
            self.msg = self.msg + txt

        self.tel = tel
        if len(date) == 1:
            self.flags.append("start_at_sunset")
            self.date = date[0] + " 16:00:00"
            txt = f'start at sunset\n'
        elif len(date) == 2:
            self.date = date[0] + " " + date[1]
            txt = f'start at {self.date}\n'
        else:
            txt = f'date missing\n'

        if len(date[0].split("/")[0])<3:
            print("wrong date format: yyyy/mm/dd")
            sys.exit(0)

        self.wind = wind
        self.plan = []

        if self.loud: sys.stderr.write(txt)
        self.msg = self.msg + txt

        # --------------- MAIN SCRIPT --------------


        self.Initiate()  # initiate script variables
        self.LoadObjects()  # loads objects to observe
        self.MakeTime()  # prepares time variable
        self.CalcObject()  # calculates visibilities for all objects
        self.MaskMoon()
        self.MaskWind()
        self.MaskPhase()
        #self.Waga()  # calculate propability for randomization order of objects to observe
        #self.RandomizeList()  # randomize listo of objects to observev
        self.MaskStartEnd()  # mask visibilite base on t_start and t_end
        self.MaskPhaseStartEnd()  # mask visibilite base on t_start and t_end
        self.Waga()  # calculate propability for randomization order of objects to observe
        self.RandomizeList()  # randomize listo of objects to observev
        self.allocate()  # SCHEDULING OBJECTS
        self.export()  # export plan


    def Initiate(self):

        self.ctc_active = True

        # plan                      # 1 minute resolution
        self.nightTime = []  # time in ephem.Date format
        self.nightObject = []  # index of objects to observe
        self.nightFree = []  # Bool marker of allocation

        # objects to obserrve
        self.name = []
        self.ra = []
        self.dec = []
        self.seq = []
        self.comment = []
        self.slotTime = []  # minutes of allocation for object
        self.spec = []  # specification for NV to observe
        self.tstart = []  #
        self.tend = []  #
        self.obflag = []  #
        self.uobi = []

        self.pstart = []
        self.pend = []
        self.period = []
        self.hjd0 = []

        self.cycle = []

        self.ph_mask = []

        # additional
        self.random_index_list = []  # list with random roder of object indices to oberve
        self.pliki_hjd = []  # files with phase coverage for objects

        self.probability = []  # propability for object list randomization
        self.visibility = []  # visibility of each object (self.nightTime indicies)
        self.priority = []  # priority for observations
        self.cycle = [] # specifcation how often to observe, e.g. cycle=3 - every three days

    def LoadObjects(self):  # Load Object file
        # first try patch as it is, if file not found, take base path form __FILE__ and append config file path using os.path.join
        if not os.path.isfile(self.cfg[self.tel]["master_file"]):
            self.cfg[self.tel]["master_file"] = os.path.join(self.script_dir, self.cfg[self.tel]["master_file"])
            if not os.path.isfile(self.cfg[self.tel]["master_file"]):
                sys.stderr.write(f'No Objects Catalog file {self.cfg[self.tel]["master_file"]}\n')
                sys.stderr.write('Check config.yaml file\n')
                sys.exit()

        if os.path.exists(self.cfg["ctc"]["ctc_base_folder"]):
            self.ctc = CycleTimeCalc(telescope=self.tel,base_folder=self.cfg["ctc"]["ctc_base_folder"], tpg=True)
            self.ctc.set_rm_modes(self.cfg["ctc"]["rm_modes_mhz"])
            self.ctc.set_start_rmode(self.cfg[self.tel]["def_read_mode"])
            self.ctc_active = True
        else:
            self.ctc_active = False

        plik = open(self.cfg[self.tel]["master_file"], 'r')
        for line in plik:
            line = line.strip()
            if line:  # no empty lines
                if line.strip()[0] != "#":
                    flag = []
                    seq = "1/Ic/30"
                    priority = 0
                    comment = ""

                    overhead = self.cfg[self.tel]["overhead"]
                    ts = False
                    te = False
                    ps = False
                    pe = False
                    period = False
                    hjd0 = False
                    ph_mk = False
                    cycle = False
                    read_mod = None
                    uobi = False
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

                    if "t_start=" in line:
                        ts = line.split("t_start=")[1].split()[0]

                    if "t_end=" in line:
                        te = line.split("t_end=")[1].split()[0]

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


                    for i, tmp in enumerate(line.split()):

                        if i == 0:
                            name = line.split()[0]
                        elif i == 1:
                            ra = line.split()[1]
                        elif i == 2:
                            dec = line.split()[2]

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
                        tmp = seq.replace("INF","1")
                    else:
                        tmp = seq

                    program = f"OBJECT  {name} {ra} {dec} seq={tmp}"
                    if read_mod is not None:
                        program = program + f" read_mod={read_mod}"

                    if self.ctc_active:
                        self.ctc.reset_time()
                        ctc_ob_time = self.ctc.calc_time(program)/60.
                        slotTime = ctc_ob_time
                        txt = f'{name} CTC slot time: {slotTime}'
                    else:
                        slotTime = calc_slot_time(program,overhead)/60.
                        txt = f'{name} TPG slot time: {slotTime}\n'

                    if self.loud: sys.stderr.write(txt)
                    self.msg = self.msg + txt

                    if uobi and uobi in self.done:
                        pass
                    else:
                        self.name.append(name)
                        self.ra.append(ra)
                        self.dec.append(dec)
                        self.seq.append(seq)
                        self.comment.append(comment)
                        self.slotTime.append(slotTime)
                        self.spec.append(line)
                        self.priority.append(priority)
                        self.tstart.append(ts)
                        self.tend.append(te)
                        self.pstart.append(ps)
                        self.pend.append(pe)
                        self.period.append(period)
                        self.hjd0.append(hjd0)
                        self.ph_mask.append(ph_mk)
                        self.cycle.append(cycle)
                        self.uobi.append(uobi)

                        self.obflag.append(flag)

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

        while t < self.sunrise:
            t = t + 1 * ephem.minute  # time resolution : 1 minute
            self.nightTime.append(t)
            self.nightObject.append(-100)
            self.nightFree.append(True)


    def CalcObject(self):  # Colculate visibility for all objects
        for n, tmp in enumerate(self.name):

            s = ephem.FixedBody()
            s._ra = self.ra[n]
            s._dec = self.dec[n]
            self.oca.date = self.nightTime[0]

            self.oca.horizon = str(self.cfg[self.tel]["hmin"])
            try:
                tr = self.oca.next_rising(s)
                ts = self.oca.next_setting(s)
                if tr < ts:
                    vis1 = [x > tr and x < ts for x in self.nightTime]
                else:
                    vis1 = [x < ts or x > tr for x in self.nightTime]
            except ephem.NeverUpError:
                #print("newer up", self.name[n])
                vis1 = [False for _ in range(len(self.nightTime))]
            except ephem.AlwaysUpError:
                vis1 = [True for _ in range(len(self.nightTime))]

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

            vis = [a and b for a,b in zip(vis1,vis2)]
            self.visibility.append(vis)

            vis = numpy.array(vis)
            txt = f'{self.name[n]} visible for {len(vis[vis])} minutes\n'
            if self.loud: sys.stderr.write(txt)
            self.msg = self.msg + txt


    def MaskMoon(self):
        moon = ephem.Moon()
        for n, tmp in enumerate(self.name):
            s = ephem.FixedBody()
            s._ra = self.ra[n]
            s._dec = self.dec[n]
            tmp = numpy.array(self.nightTime)
            time = tmp[numpy.array(self.visibility[n])]
            if len(time)>0:
                self.oca.date = time[0]
                moon.compute(self.oca)
                s.compute(self.oca)
                sep = ephem.separation(s,moon) * 180/ephem.pi
                if sep < float(self.cfg[self.tel]["moon_dist"]):
                    # mask everything for faster calculation
                    self.visibility[n] = [False for _ in range(len(self.nightTime))]
                    txt = f'{self.name[n]} masked because of MOON\n'
                    if self.loud: sys.stderr.write(txt)
                    self.msg = self.msg + txt



    def MaskWind(self):
        if self.wind:
            for n, tmp in enumerate(self.name):
                s = ephem.FixedBody()
                s._ra = self.ra[n]
                s._dec = self.dec[n]
                tmp = numpy.array(self.nightTime)
                time = tmp[numpy.array(self.visibility[n])]
                if len(time)>0:
                    self.oca.date = time[0]
                    s.compute(self.oca)
                    az = s.az * 180/ephem.pi
                    dist = numpy.abs(self.wind - az)
                    if dist > 180:
                        dist = 360. - dist
                    if dist < float(self.cfg[self.tel]["wind_dist"]):
                        self.visibility[n] = [False for _ in range(len(self.nightTime))]
                        # masked for all night because of calc time
                        txt = f'{self.name[n]} masked because of WIND\n'
                        if self.loud: sys.stderr.write(txt)
                        self.msg = self.msg + txt


    def Waga(self):  #
        for n, tmp in enumerate(self.name):
            maska = numpy.array(self.visibility[n])
            maska = maska[maska]
            allminutes = len(self.nightTime)
            if len(maska) == 0:
                waga = 0
            else:
                waga = 1. - float(len(maska)) / float(allminutes)
            if waga < 0.3 and len(maska) != 0: waga = 0.3
            self.probability.append(waga)

    def RandomizeList(self):  # randomize listo of objects to observev
        mozliwosci = numpy.array([i for i, _ in enumerate(self.name)])
        probability = numpy.array(self.probability)
        maska = (probability != 0)
        mozliwosci = mozliwosci[maska]
        probability = probability[maska]
        for n, tmp in enumerate(mozliwosci):
            suma = sum(probability)
            probability = [x / suma for x in probability]
        self.random_index_list = numpy.random.choice(mozliwosci, len(mozliwosci), replace=False, p=probability)

    def allocate(self):
        for nn, tmp in enumerate(self.nightTime):
            if self.nightFree[nn]:

                OK = False
                done = -1
                priormax = -1

                for nk, tmp in enumerate(self.random_index_list):
                    k = self.random_index_list[nk]
                    if int(self.priority[k]) > priormax and int(self.priority[k]) > -50:

                        visibility = self.visibility[k]
                        slots = float((self.slotTime[k])) / 1.
                        if "inf" in self.obflag[k]:
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
                            priormax = int(self.priority[k])
                            bfin = bounds
                            kfin = k
                            OK = True
                            done = nk




                if OK:
                    self.nightFree[bfin[0]:bfin[1] + 1] = [False] * ((bfin[1] + 1) - bfin[0])
                    self.nightObject[bfin[0]:bfin[1] + 1] = [int(kfin)] * ((bfin[1] + 1) - bfin[0])
                    q = self.random_index_list[done]
                    if self.cycle[q]:
                        dn = int(1440 * float(self.cycle[q]))
                        if dn > len(self.visibility[q][nn:]) : dn = len(self.visibility[q][nn:])
                        self.visibility[q][nn:nn+dn] = [False] * dn
                    else:
                        self.random_index_list = numpy.delete(self.random_index_list, done)
                    if "inf" in self.obflag[kfin]:
                        self.seq[kfin] = self.seq[kfin].replace("INFx", f"{loops}x")
                        if "loud" in self.flags:
                            txt = f" {self.name[kfin]} INF changed to {loops}"
                            sys.stderr.write(txt + "\n")


    def export(self):
        if True:
            prev = -100
            ZNACZNIK = False
            if "start_at_sunset" in self.flags:
                txt = f'WAIT sunset={str(self.cfg[self.tel]["hsun"])}'
            else:
                txt = f'WAIT ut={str(ephem.Date(self.nightTime[0])).split()[1]}'
            self.plan.append(txt)
            for n, tmp in enumerate(self.nightTime):
                #print(self.nightTime[n], self.nightObject[n])
                if int(self.nightObject[n]) != prev:
                    if ZNACZNIK:
                        ZNACZNIK = False
                        empty_time = float(self.nightTime[n] - t0) * 24. * 60. * 60.
                        self.oca.date = ephem.Date(self.nightTime[n])
                        txt = "WAIT ut=" + str(self.oca.date).split()[1]
                        self.plan.append(txt)
                    prev = int(self.nightObject[n])
                    if prev > -50:
                        k = int(self.nightObject[n])

                        if "x" in self.seq[k]:
                            q = self.seq[k].split("x")[0]
                            r = self.seq[k].split("x")[1]
                            if int(q) > 10:
                                self.seq[k] = f"10x{r}"
                                for _ in range(int(int(q)/10)):
                                    txt = f"OBJECT  {self.name[k]} {self.ra[k]} {self.dec[k]}      seq={self.seq[k]}"
                                    if len(self.comment[k]) > 0: txt = txt + f"    comment=\"{self.comment[k]}\""
                                    self.plan.append(txt)
                                z = int(q)-10*int(int(q)/10)
                                if z != 0:
                                    self.seq[k] = f"{z}x{r}"
                                    txt = f"OBJECT  {self.name[k]} {self.ra[k]} {self.dec[k]}      seq={self.seq[k]}"
                                    if len(self.comment[k]) > 0: txt = txt + f"    comment=\"{self.comment[k]}\""
                                    self.plan.append(txt)
                            else:
                                txt = f"OBJECT  {self.name[k]} {self.ra[k]} {self.dec[k]}      seq={self.seq[k]}"
                                if len(self.comment[k]) > 0: txt = txt + f"    comment=\"{self.comment[k]}\""
                                if self.uobi[k]:
                                    txt = txt + f" uobi={self.uobi[k]}"
                                self.plan.append(txt)
                        else:
                            txt = f"OBJECT  {self.name[k]} {self.ra[k]} {self.dec[k]}      seq={self.seq[k]}"
                            if self.uobi[k]:
                                txt = txt + f" uobi={self.uobi[k]}"
                            if len(self.comment[k]) > 0: txt = txt + f"    comment=\"{self.comment[k]}\""
                            self.plan.append(txt)

                    else:
                        t0 = self.nightTime[n]
                        ZNACZNIK = True

    def MaskStartEnd(self):  # Masks visibility based on t_start t_end
        for n, tmp in enumerate(self.name):
            ZNACZNIK = False

            ts = self.tstart[n]
            if ts:
                if "T" in ts:
                    ts = ts.replace("T", " ")
                if "/" not in ts:
                    date = self.date.split()[0]
                    if int(ts.split(":")[0]) < 16:
                        date = ephem.Date(date) + 1
                        date = str(ephem.Date(date)).split()[0]
                    ts = date + " " + ts
                ts = ephem.Date(ts)
                time_s = ts
                if time_s or time_e:
                    for i, tmp in enumerate(self.nightTime):
                        if float(self.nightTime[i]) < float(time_s) and self.visibility[n][i] and time_s:
                            self.visibility[n][i] = False
                            ZNACZNIK = True


            te = self.tend[n]
            if te:
                if "T" in te:
                    te = te.replace("T", " ")
                if "/" not in te:
                    date = self.date.split()[0]
                    if int(te.split(":")[0]) < 16:
                        date = ephem.Date(date) + 1
                        date = str(ephem.Date(date)).split()[0]
                    te = date + " " + te
                te = ephem.Date(te)
                time_e = te
                if time_e < time_s:
                    time_e = time_e + ephem.day

                if time_s or time_e:
                    for i, tmp in enumerate(self.nightTime):
                        if float(self.nightTime[i]) > float(time_e) and self.visibility[n][i] and time_e:
                            self.visibility[n][i] = False
                            ZNACZNIK = True

                if ZNACZNIK:
                    txt = f'{self.name[n]} time masked {time_s} to {time_e}\n'
                    if self.loud: sys.stderr.write(txt)
                    self.msg = self.msg + txt

    def MaskPhaseStartEnd(self):  # Masks visibility based on pstart and pend
        for n, tmp in enumerate(self.name):
            ZNACZNIK = False

            ps = self.pstart[n]
            pe = self.pend[n]
            period = self.period[n]
            hjd0 = self.hjd0[n]

            if ps:
                for i, t in enumerate(self.nightTime):
                    if self.visibility[n][i]:
                        jd = ephem.julian_date(ephem.Date(t))
                        p = (jd - float(hjd0)) / float(period) % 1
                        if pe < ps:
                            if p > float(pe) and p < float(ps):
                                self.visibility[n][i] = False
                                ZNACZNIK = True
                            else:
                                pass
                        else:
                            if p > float(ps) and p < float(pe):
                                pass
                            else:
                                self.visibility[n][i] = False
                                ZNACZNIK = True

                if ZNACZNIK:
                    txt = f'{self.name[n]} phase masked {ps} to {pe}\n'
                    if self.loud: sys.stderr.write(txt)
                    self.msg = self.msg + txt

    def MaskPhase(self):  # mask phase covarage
        if True:
            for n, tmp in enumerate(self.name):
                ZNACZNIK = False
                if self.ph_mask[n]:
                    jd = []
                    try:
                        filter = self.ph_mask[n].split("/")[1]
                        file = self.cfg[self.tel]["data_file"] + self.name[n].lower() + "/" + filter + "/light-curve/" + self.name[n].lower() + "_" + filter + "_diff_light_curve.txt"
                        with open(file, "r") as plik:
                            if plik != None:
                                for line in plik:
                                    if len(line.strip()) > 0:
                                        try:
                                            if len(line.split())>9:
                                                if int(line.split()[9]) < 2:
                                                    jd.append(float(line.split()[3]))
                                            else:
                                                jd.append(float(line.split()[3]))
                                        except ValueError:
                                            pass

                    except (FileNotFoundError, ValueError) as e:
                        pass
                    P = self.period[n]
                    hjd0 = self.hjd0[n]
                    if P and len(jd)>0:
                        if not hjd0:
                            hjd0 = 2460000
                        ph = (numpy.array(jd) - float(hjd0)) / float(P) % 1

                        night_ph = (numpy.array(self.nightTime) - float(hjd0)) / float(P) % 1

                        for i,tmp in enumerate(night_ph):
                            np1 = night_ph[i] - float(self.ph_mask[n].split("/")[2])
                            np2 = night_ph[i] + float(self.ph_mask[n].split("/")[2])
                            maska1 = ph > np1
                            maska2 = ph < np2
                            maska = [a and b for a,b in zip(maska2,maska1)]
                            if len(ph[maska]) > int(self.ph_mask[n].split("/")[0]):
                                self.visibility[n][i] = False
                                ZNACZNIK = True
                if ZNACZNIK:
                    txt = f'{self.name[n]} masked for observed data\n'
                    if self.loud: sys.stderr.write(txt)
                    self.msg = self.msg + txt


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