import sys,os,yaml,numpy,ephem
import random
import time
import copy

from astropy.table import Table
from tqdm import tqdm

from pyaraucaria.obs_plan.obs_plan_parser import ObsPlanParser
from pyaraucaria.ob_validator import ObsValidator

try:
    from ctc import CycleTimeCalc
    CTC_AVAILABLE = True
except ImportError:
    CTC_AVAILABLE = False


class TelescopePlanGenerator():
    def __init__(self,tel,date,tag=None,wind=None,fwhm=None,loud=False,seed=None,add_start_makro=True,add_end_makro=True,save_plan=True, ignore_ctc=False, master_file=None, done_uobi=[]):
        self.cwd = os.getcwd()
        self.script_dir = os.path.dirname(os.path.abspath(__file__))
        self.arg_date = date
        self.flags = []
        self.tel = tel
        self.tag = tag
        self.wind = wind
        self.fwhm = fwhm
        self.done = done_uobi
        self.loud = loud
        self.seed = seed
        self.add_start_makro = add_start_makro
        self.add_end_makro = add_end_makro
        self.save_plan = bool(save_plan)
        self.ignore_ctc = bool(ignore_ctc)
        self.master_file = master_file

        self.msg = ""

        # --------------- MAIN SCRIPT --------------
    def go(self):
        self.Initiate()  # initiate script variables
        self.LoadObjects()  # loads objects to observe
        self.MakeTime()  # prepares time variable
        self.ObjectMask() # maskowanie tag, fwhm, status
        self.CalcObject()  # calculates visibilities for all objects
        self.MaskVisibility()
        self.MaskMoon()
        self.MaskWind()
        self.MaskTwilight()
        self.MaskCycle()
        self.MaskStartEnd()  # mask visibilite base on t_start and t_end
        self.MaskPhaseStartEnd()  # mask visibilite base on t_start and t_end
        self.MaskPhase()
        self.Waga()  # calculate probability for randomization order of objects to observe
        self.RandomizeList()  # randomize listo of objects to observev
        self.allocate()  # SCHEDULING OBJECTS
        self._printObjectLog()
        self.export()  # export plan
        self.SavePlan()


    def log_msg(self,txt,always_loud=True):
        if always_loud:
            print(txt)
        else:
            if self.loud:
                print(txt)
        self.msg = self.msg + txt + "\n"

    def _printObjectLog(self):
        txt = ""
        for ob in self.ob:
            cm = ob["command_name"]
            name = ob["name"]
            slot = ob.get("slotTime", "--")
            try:
                slot = float(slot)
                slot = f'{slot:.0f}'
            except:
                pass
            txt = txt + f'{cm} {name} [{slot}m]'

            for kw, vis in ob["visibility"].items():
                vis_array = numpy.array(vis)
                hours_visible = vis_array.sum()/60.
                txt = txt + f'  {kw}: {hours_visible:.1f}'

            txt = txt + "\n"
        self.log_msg(txt,always_loud=False)

    def _get_config(self):

        with open(self.script_dir+'/config.yaml', 'r') as cfg_file:
            self.cfg = yaml.safe_load(cfg_file)
            txt = f'config loaded from {self.script_dir}/config.yaml'
            self.log_msg(txt)

        local_config_file = self.cfg["local_config"]
        local_config_file = os.path.expanduser(local_config_file)
        if os.path.exists(local_config_file):
            with open(local_config_file) as f:
                local_cfg = yaml.safe_load(f) or {}
                self.cfg.update(local_cfg)
                txt = f'local config loaded from {local_config_file}'
                self.log_msg(txt)

    def _parse_date(self):

        if len(self.arg_date) == 1:          # bo to lista
            self.flags.append("start_at_sunset")
            self.date = self.arg_date[0] + " 12:00:00"
            txt = f'start at sunset\n'
        elif len(self.arg_date) == 2:
            self.date = self.arg_date[0] + " " + self.arg_date[1]
            txt = f'start at {self.date}\n'
        else:
            print("date missing: yyyy/mm/dd")
            sys.exit(0)

        if "/" in self.arg_date[0]:
            if len(self.arg_date[0].split("/"))<3:
                print("wrong date format: yyyy/mm/dd")
                sys.exit(0)
            if len(self.arg_date[0].split("/")[0])<4:
                print("wrong date format: yyyy/mm/dd")
                sys.exit(0)
        elif "-" in self.arg_date[0]:
            if len(self.arg_date[0].split("-")) < 3:
                print("wrong date format: yyyy/mm/dd")
                sys.exit(0)
            if len(self.arg_date[0].split("-")[0])<4:
                print("wrong date format: yyyy/mm/dd")
                sys.exit(0)
        else:
            print("wrong date format: yyyy/mm/dd")
            sys.exit(0)

        self.log_msg(txt)

    def Initiate(self):

        self._get_config()
        self._parse_date()

        self.ctc_active = True
        self.plan_saved = False

        # plan                      # 1 minute resolution
        self.nightTime = []  # time in ephem.Date format
        self.nightObject = []  # index of objects to observe
        self.nightFree = []  # Bool marker of allocation
        self.nightObChange = []    # Bool marker to indicate ob change

        self.ob = []
        self.plan = []

        # additional
        self.random_index_list = []  # list with random order of object indices to observe
        self.pliki_hjd = []  # files with phase coverage for objects

        BASE_SCHEMA = ObsValidator.load_schema("base_schema.yaml")
        TPG_SCHEMA = ObsValidator.load_schema("tpg_schema.yaml")

        self.schema = merge_schemas(BASE_SCHEMA, TPG_SCHEMA)

        self.plan_rules = ObsValidator.load_schema("base_rules.yaml")
        self.validator = ObsValidator(self.schema, self.plan_rules)

    def init_ctc(self):
        if self.ignore_ctc:
            self.ctc_active = False
            return
        if os.path.exists(self.cfg["ctc"]["ctc_base_folder"]) and CTC_AVAILABLE:
            self.ctc = CycleTimeCalc(telescope=self.tel,base_folder=self.cfg["ctc"]["ctc_base_folder"], tpg=True)
            self.ctc.set_rm_modes(self.cfg["ctc"]["rm_modes_mhz"])
            self.ctc.set_start_rmode(self.cfg[self.tel]["def_read_mode"])
            self.ctc_active = True
        else:
            self.ctc_active = False

        txt = f'ctc active: {self.ctc_active}\n'
        self.log_msg(txt)

    def LoadObjects(self):  # Load Object file
        # first try patch as it is, if file not found, take base path form __FILE__ and append config file path using os.path.join
        if not os.path.isfile(self.cfg[self.tel]["master_file"]):
            self.cfg[self.tel]["master_file"] = os.path.join(self.script_dir, self.cfg[self.tel]["master_file"])
            if not os.path.isfile(self.cfg[self.tel]["master_file"]):
                self.log_msg(f'No Objects Catalog file {self.cfg[self.tel]["master_file"]}')
                self.log_msg('Check config.yaml file')
                sys.exit()

        self.init_ctc()

        if self.master_file:
            fname = self.master_file
        else:
            fname = self.cfg[self.tel]["master_file"]
        try:
            with open(fname, 'r') as file:
                for line in tqdm(file, desc="Loading objects"):
                    if not line or line.strip().startswith("#"):
                        continue
                    tmp_ob = self.parseObjects(line)
                    if tmp_ob:
                        self.ob.append(tmp_ob)


        except FileNotFoundError:
            self.log_msg(f"File not found: {fname}")
            sys.exit(1)



    def parseObjects(self,line):

        ob = None
        ob_ok = None
        if len(line.strip()) > 0:
            if len(line.split()) > 0:
                if "#" not in line.split()[0]:
                    ob_tmp = ObsPlanParser.convert_from_string(line)

                    if ob_tmp is None:
                        txt = f'** parser error \n * {line}'
                        self.log_msg(txt)
                    else:
                        ob = ObsValidator.convert_to_obdict(ob_tmp)
                        validator_result = self.validator.validate_ob(ob)
                        if bool(validator_result["valid"]):
                            ob_ok = True
                        else:
                            ob_ok = False
                            txt = f'** validator error:  '
                            for key, value in validator_result["result"].items():
                                if not value:
                                    txt = txt + f'{key} '
                            txt = txt + "\n" + f'* {line}'
                            self.log_msg(txt)


        if ob != None and ob_ok:

            # uzupelnienie do liczenia fazowania
            if "P" in ob.keys() and "hjd0" not in ob.keys():
                ob["hjd0"] = 2460000

            # slotTime calculation
            if "ob_time" in ob.keys():
                ob["slotTime"] = ob["ob_time"]
            elif "seq" in ob.keys():
                overhead = self.cfg[self.tel]["overhead"]
                seq = ob["seq"]
                slotTime = 10
                k = 1
                if "x(" in seq and ")" in seq:
                    k = seq.split("x")[0]
                    k = int(k)
                    seq0 = seq.split("x(")[1].split(")")[0]
                else:
                    seq0 = seq
                for x_seq in seq0.split(","):
                    slotTime = slotTime + (
                            float(x_seq.split("/")[0]) * (float(x_seq.split("/")[2]) + float(overhead)))
                slotTime = (k * slotTime) / 60.

                program = f'{ob["command_name"]} {ob["name"]} {ob["ra"]} {ob["dec"]} seq={seq}'
                if "read_mod" in ob.keys():
                    read_mod = ob["read_mod"]
                    program = program + f" read_mod={read_mod}"
                if self.ctc_active:
                    self.ctc.reset_time()
                    ctc_ob_time = self.ctc.calc_time(program) / 60.
                    slotTime = ctc_ob_time
                else:
                    slotTime = calc_slot_time(program, overhead) / 60.
                ob["slotTime"] = slotTime

            else:
                ob["slotTime"] = 10

        else:
            ob = False

        return ob

    def MakeTime(self):
        self.oca = ephem.Observer()
        self.oca.lon = str(self.cfg["obs_lon"])
        self.oca.lat = str(self.cfg["obs_lat"])
        self.oca.elev = float(self.cfg["obs_elev"])
        self.oca.horizon = str(self.cfg[self.tel]["hsun"])

        sun = ephem.Sun()
        t0 = ephem.Date(self.date)
        self.oca.date = t0

        # zachód i wschód tej nocy
        next_sunset = self.oca.next_setting(sun, use_center=True)
        next_sunrise = self.oca.next_rising(sun, use_center=True)
        mid_night = (next_sunset + next_sunrise) / 2

        now = ephem.now()

        if now < next_sunset:
            # przed zachodem → bierzemy dzisiejszy zachód
            self.sunset = next_sunset
            self.oca.date = self.sunset
            self.sunrise = next_sunrise
        elif now < mid_night:
            # po zachodzie, ale przed środkiem nocy → bierzemy poprzedni zachód
            self.sunset = self.oca.previous_setting(sun, use_center=True)
            self.oca.date = self.sunset
            self.sunrise = self.oca.next_rising(sun, use_center=True)
        else:
            # po środku nocy → bierzemy kolejny zachód
            self.oca.date = now
            self.sunset = self.oca.next_setting(sun, use_center=True)
            self.oca.date = self.sunset
            self.sunrise = self.oca.next_rising(sun, use_center=True)

        # ustaw start nocy
        if "start_at_sunset" in self.flags:
            t = self.sunset
        else:
            t = t0

        self.start_time = t

        # generowanie listy minut nocy
        self.nightTime = []
        self.nightObject = []
        self.nightObChange = []
        self.nightFree = []

        current = t
        while current < self.sunrise:
            current = current + ephem.minute
            self.nightTime.append(current)
            self.nightObject.append(-100)
            self.nightObChange.append(False)
            self.nightFree.append(True)

    def ObjectMask(self):

        for n, tmp in enumerate(tqdm(self.ob, desc="Selecting objects")):

            self.ob[n]["visibility"] = dict()
            self.ob[n]["visibility"]["all"] = [True for _ in range(len(self.nightTime))]

            # sprawdzanie statusu
            if "status" in self.ob[n].keys():
                if self.ob[n]["status"] in ["finished", "waiting", "halted", "inactive"]:
                    vis = [False for _ in range(len(self.nightTime))]
                    self.ob[n]["visibility"]["all"] = vis.copy()
                    self.ob[n]["visibility"]["status"] = vis.copy()

            # sprawdzanie opcji tag
            ok = True
            if self.tag:
                ok = False
                if d_ver(self.ob[n], ["tag"]):
                    if self.tag in self.ob[n]["tag"]:
                        ok = True
            if ok:
                pass
            else:
                vis = [False for _ in range(len(self.nightTime))]
                self.ob[n]["visibility"]["tag"] = vis.copy()
                self.ob[n]["visibility"]["all"] = self.ob[n]["visibility"]["tag"].copy()

            # sprawdzanie fwhm
            if "max_fwhm" in self.ob[n].keys() and self.fwhm is not None:
                if float(self.ob[n]["max_fwhm"]) > float(self.fwhm):
                    vis = [False for _ in range(len(self.nightTime))]
                    self.ob[n]["visibility"]["max_fwhm"] = vis.copy()
                    self.ob[n]["visibility"]["all"] = self.ob[n]["visibility"]["max_fwhm"].copy()

            if "min_fwhm" in self.ob[n].keys() and self.fwhm is not None:
                if float(self.ob[n]["min_fwhm"]) < float(self.fwhm):
                    vis = [False for _ in range(len(self.nightTime))]
                    self.ob[n]["visibility"]["min_fwhm"] = vis.copy()
                    self.ob[n]["visibility"]["all"] = self.ob[n]["visibility"]["min_fwhm"].copy()

    def CalcObject(self):
        moon = ephem.Moon()

        for ob in tqdm(self.ob, desc="calculating objects"):
            if not d_ver(ob, ["name", "ra", "dec"]):
                continue

            s = ephem.FixedBody()
            s._ra = ob["ra"]
            s._dec = ob["dec"]

            alt = []
            az = []
            moon_sep = []
            moon_phase = []

            for x in self.nightTime:
                self.oca.date = x

                s.compute(self.oca)
                moon.compute(self.oca)

                alt.append(float(s.alt) * 180.0 / ephem.pi)
                az.append(float(s.az) * 180.0 / ephem.pi)
                moon_sep.append(float(ephem.separation(s, moon)) * 180.0 / ephem.pi)
                moon_phase.append(float(moon.phase))

            ob["calc"] = {"alt": alt,"az": az,"moon_sep": moon_sep,"moon_phase": moon_phase}

    def MaskVisibility(self):
        for ob in tqdm(self.ob, desc="Masking visibility"):
            if "calc" not in ob.keys():
                continue

            alt = ob["calc"]["alt"]

            hmin = float(ob.get("h_min", self.cfg[self.tel]["hmin"]))
            hmax = float(ob.get("h_max", self.cfg[self.tel]["hmax"]))

            hmin_mask = [a >= hmin for a in alt]
            hmax_mask = [a <= hmax for a in alt]
            all_mask = [m1 and m2 for m1, m2 in zip(hmin_mask, hmax_mask)]

            ob["visibility"]["h_min"] = hmin_mask
            ob["visibility"]["h_max"] = hmax_mask
            ob["visibility"]["all"] = all_mask

    def MaskMoon(self):
        for ob in tqdm(self.ob, desc="Masking Moon"):
            if "calc" not in ob.keys():
                continue

            sep = ob["calc"]["moon_sep"]
            phase = ob["calc"]["moon_phase"]

            moon_dist = float(ob.get("min_moon_dist", self.cfg[self.tel]["moon_dist"]))
            max_phase = ob.get("max_moon_phase", None)

            base_vis = ob["visibility"]["all"]

            min_moon_mask = []
            max_phase_mask = []
            all_mask = []

            for i in range(len(base_vis)):
                ok_base = base_vis[i]

                cond_dist = sep[i] >= moon_dist
                cond_phase = True if max_phase is None else phase[i] <= float(max_phase)

                min_moon_mask.append(cond_dist)
                max_phase_mask.append(cond_phase)

                all_mask.append(ok_base and cond_dist and cond_phase)

            ob["visibility"]["min_moon_dist"] = min_moon_mask
            ob["visibility"]["max_moon_phase"] = max_phase_mask
            ob["visibility"]["all"] = all_mask

    def MaskWind(self):
        if not self.wind:
            return

        for ob in tqdm(self.ob, desc="Masking Wind"):
            if "calc" not in ob.keys():
                continue

            az = ob["calc"]["az"]
            vis = ob["visibility"]["all"]

            wind_dist = float(self.cfg[self.tel]["wind_dist"])

            new_vis = []

            for i in range(len(vis)):
                ok = vis[i]

                if ok:
                    dist = abs(self.wind - az[i])
                    if dist > 180:
                        dist = 360 - dist

                    if dist < wind_dist:
                        ok = False

                new_vis.append(ok)

            ob["visibility"]["wind"] = new_vis.copy()
            ob["visibility"]["all"] = new_vis

    def MaskTwilight(self):
        for ob in tqdm(self.ob, desc="masking sunset_dt and sunrise_dt"):
            if "visibility" not in ob.keys():
                continue

            base_vis = ob["visibility"]["all"]

            # parametry (w godzinach → dni)
            after_sunset = float(ob.get("sunset_dt", 0.0)) / 24.0
            before_sunrise = float(ob.get("sunrise_dt", 0.0)) / 24.0

            after_sunset_mask = []
            before_sunrise_mask = []
            all_mask = []

            for i, t in enumerate(self.nightTime):
                ok_base = base_vis[i]

                cond_after = True
                if after_sunset > 0:
                    cond_after = t >= self.sunset + after_sunset

                cond_before = True
                if before_sunrise > 0:
                    cond_before = t <= self.sunrise - before_sunrise

                after_sunset_mask.append(cond_after)
                before_sunrise_mask.append(cond_before)

                all_mask.append(ok_base and cond_after and cond_before)

            ob["visibility"]["after_sunset"] = after_sunset_mask
            ob["visibility"]["before_sunrise"] = before_sunrise_mask
            ob["visibility"]["all"] = all_mask


    def MaskCycle(self):
        for ob in tqdm(self.ob, desc="masking cycle"):
            last_jd = None
            if not d_ver(ob, ["name", "cycle", "visibility"]):
                continue

            if not ob["cycle"]:
                continue

            jd_list = []
            try:
                if d_ver(ob, ["obs_data"]):
                    if "/" == self.cfg[self.tel]["data_file"][-1] or ob["obs_data"][0]:
                        separator = ""
                    else:
                        separator = "/"
                    file = self.cfg[self.tel]["data_file"]+separator+ob["obs_data"]
                else:
                    filtr = ob["seq"].split("/")[1]
                    file = (self.cfg[self.tel]["data_file"] + ob["name"].lower() + "/" + filtr + "/light-curve/" + ob["name"].lower() + "_" + filtr + "_diff_light_curve.txt")

                lc_tab = Table.read(file, format="ascii")
                jd = numpy.array(lc_tab["jd_obs"])
                if "quality" in lc_tab.colnames:
                    quality = numpy.array(lc_tab["quality"])
                    mask_quality = quality < 2
                    jd = jd[mask_quality]
                if len(jd) == 0:
                    continue
                last_jd = max(jd)

            except FileNotFoundError:
                txt = f'[{ob["name"]}] no data file {file}'
                self.log_msg(txt,always_loud=False)
            except ValueError:
                txt = f'[{ob["name"]}] ERROR loading {file}\n'
                self.log_msg(txt,always_loud=False)

            base_vis = ob["visibility"]["all"]

            cycle_mask = []
            all_mask = []

            if last_jd is None:
                continue

            for i, t in enumerate(self.nightTime):
                jd_t = ephem.julian_date(ephem.Date(t))

                cond_cycle = (jd_t - last_jd) >= float(ob["cycle"])

                cycle_mask.append(cond_cycle)
                all_mask.append(base_vis[i] and cond_cycle)

            ob["visibility"]["cycle"] = cycle_mask
            ob["visibility"]["all"] = all_mask

    def MaskStartEnd(self):
        for ob in tqdm(self.ob, desc="masking start and end time"):
            if "visibility" not in ob:
                continue

            t_all = numpy.array(ob["visibility"]["all"])
            nt = numpy.array(self.nightTime)

            # ---- Mask t_start ----
            if d_ver(ob, ["name", "t_start"]):
                ts = ob["t_start"]
                if ts:
                    ts = ts.replace("T", " ")
                    # jeśli ts podana jest tylko godzina
                    if "/" not in ts:
                        date_str = str(ephem.Date(self.start_time)).split()[0]
                        ts = f"{date_str} {ts}"

                    ts = ephem.Date(ts)

                    # jeśli ts przed startem nocy, dodaj 1 dzień
                    if ts < self.start_time:
                        ts += 1

                    mask_start = nt >= float(ts)
                    t_all = t_all & mask_start
                    ob["visibility"]["t_start"] = t_all.copy()

            # ---- Mask t_end ----
            if d_ver(ob, ["name", "t_end"]):
                te = ob["t_end"]
                if te:
                    te = te.replace("T", " ")
                    if "/" not in te:
                        date_str = str(ephem.Date(self.start_time)).split()[0]
                        te = f"{date_str} {te}"

                    te = ephem.Date(te)

                    # jeśli te przed startem nocy, dodaj 1 dzień
                    if te < self.start_time:
                        te += 1

                    mask_end = nt <= float(te)
                    t_all = t_all & mask_end
                    ob["visibility"]["t_end"] = t_all.copy()

            # Zaktualizuj "all" na końcu
            ob["visibility"]["all"] = t_all

    def MaskPhaseStartEnd(self):
        for ob in tqdm(self.ob, desc="masking phase start and end"):
            if not d_ver(ob, ["name", "ph_start", "ph_end", "visibility", "P", "hjd0"]):
                continue

            ps = float(ob["ph_start"])
            pe = float(ob["ph_end"])
            period = float(ob["P"])
            hjd0 = float(ob["hjd0"])

            # --- Obsługa jakości danych (opcjonalna) ---
            jd_good = None
            if d_ver(ob, ["obs_data"]):
                file = self.cfg[self.tel]["data_file"] + "/" + ob["obs_data"]
                try:
                    lc_tab = Table.read(file, format="ascii")
                    jd = numpy.array(lc_tab["jd_obs"])
                    jd = jd[~numpy.isnan(jd)]
                    if "quality" in lc_tab.colnames:
                        quality = numpy.array(lc_tab["quality"])
                        mask_quality = quality < 2
                        jd = jd[mask_quality]
                    if len(jd) > 0:
                        jd_good = jd
                except (FileNotFoundError, ValueError):
                    pass

            # Jeśli brak dobrych danych, pomijamy maskę
            if jd_good is None:
                continue

            t_all = numpy.array(ob["visibility"]["all"])
            nt = numpy.array(self.nightTime)
            jd_nt = numpy.array([ephem.julian_date(ephem.Date(tt)) for tt in nt])

            # Obliczamy fazę
            phase = ((jd_nt - hjd0) / period) % 1

            # Tworzymy maskę logiczną
            if pe < ps:
                mask_phase = ~((phase > pe) & (phase < ps))  # False = nieobserwowalne
            else:
                mask_phase = (phase >= ps) & (phase <= pe)

            # Nakładamy maskę
            t_all = t_all & mask_phase
            ob["visibility"]["ph_start"] = t_all.copy()
            ob["visibility"]["ph_end"] = t_all.copy()
            ob["visibility"]["all"] = t_all

    def MaskPhase(self):
        for ob in tqdm(self.ob, desc="masking phase coverage"):
            if not d_ver(ob, ["name", "ph_mk", "visibility", "P", "hjd0"]):
                continue

            if not ob["ph_mk"]:
                continue

            # rozpakowanie ph_mk: "min_count/filter/phase_width"
            try:
                ph_parts = ob["ph_mk"].split("/")
                min_count = int(ph_parts[0])
                filt = ph_parts[1]
                phase_width = float(ph_parts[2])
            except:
                continue  # zła specyfikacja ph_mk

            # --- wczytanie danych i filtr jakości ---
            jd = None
            if d_ver(ob, ["obs_data"]):
                file = self.cfg[self.tel]["data_file"] + "/" + ob["obs_data"]
            else:
                file = (self.cfg[self.tel]["data_file"] + ob["name"].lower() + "/" + filt + "/light-curve/" + ob["name"].lower() + "_" + filt + "_diff_light_curve.txt")

            try:
                lc_tab = Table.read(file, format="ascii")
                jd = numpy.array(lc_tab["jd_obs"])
                jd = jd[~numpy.isnan(jd)]
                if "quality" in lc_tab.colnames:
                    quality = numpy.array(lc_tab["quality"])
                    mask_quality = quality < 2
                    jd = jd[mask_quality]
                if len(jd) == 0:
                    continue

            except FileNotFoundError:
                txt = f'[{ob["name"]}] no data file {file}'
                self.log_msg(txt,always_loud=False)
                continue
            except ValueError:
                txt = f'[{ob["name"]}] ERROR loading {file}\n'
                self.log_msg(txt,always_loud=False)
                continue

            P = float(ob["P"])
            hjd0 = float(ob["hjd0"])

            # --- obliczenie faz ---
            ph = (jd - hjd0) / P % 1
            night_ph = (numpy.array(self.nightTime) - hjd0) / P % 1

            # --- przygotowanie maski dla każdego momentu nightTime ---
            t_all = numpy.array(ob["visibility"]["all"])
            for i, np_i in enumerate(night_ph):
                np_min = np_i - phase_width
                np_max = np_i + phase_width
                # liczymy ile punktów w danych obserwacji jest w tym oknie
                in_window = (ph >= np_min) & (ph <= np_max)
                if in_window.sum() >= min_count:
                    t_all[i] = False  # maskujemy
            # zapisujemy
            ob["visibility"]["ph_mk"] = t_all.copy()
            ob["visibility"]["all"] = t_all

    def Waga(self):
        any_observable = False

        for ob in self.ob:
            if not d_ver(ob, ["name", "visibility"]):
                continue

            mask_all = numpy.array(ob["visibility"]["all"])
            visible_count = mask_all.sum()
            total_count = len(self.nightTime)

            if visible_count == 0:
                waga = 0
            else:
                any_observable = True
                waga = 1.0 - float(visible_count) / float(total_count)
                if waga < 0.3:
                    waga = 0.3

            ob["waga"] = waga

        if not any_observable:
            raise ValueError("No observable objects!")

    def RandomizeList(self):
        # --- obiekty możliwe do losowania ---
        indices = numpy.arange(len(self.ob))
        wagi = numpy.array([ob.get("waga", 0) for ob in self.ob])

        mask = wagi > 0
        indices = indices[mask]
        wagi = wagi[mask]

        if len(indices) == 0:
            self.random_index_list = numpy.array([], dtype=int)
            return  # brak obiektów do losowania

        # --- normalizacja prawdopodobieństwa ---
        probability = wagi / wagi.sum()

        # --- ustawienie ziarna ---
        if not getattr(self, "seed", None):
            self.seed = random.randrange(0, 2 ** 32)

        txt = f"random seed: {self.seed}\n"
        self.log_msg(txt)

        numpy.random.seed(self.seed)

        # --- losowanie indeksów ---
        self.random_index_list = numpy.random.choice(indices, size=len(indices), replace=False, p=probability)

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

                    # usuwanie juz zaalocowanych
                    q = self.random_index_list[done]

                    if d_ver(self.ob[q], ["uobi"]):
                        uobi_tmp = self.ob[q]["uobi"]
                        if len(uobi_tmp.strip()) > 0:
                            indeksy = [i for i, item in enumerate(self.ob) if item.get("uobi") == uobi_tmp]
                            ni = [j for j, val in enumerate(self.random_index_list) if val in indeksy]
                    else:
                        ni = [done]

                    del_list = []
                    for ki in ni:
                        q = self.random_index_list[ki]
                        if self.ob[q].get("cycle",None):
                            dn = int(1440 * float(self.ob[q]["cycle"]))  # 1440 minut ma doba
                            if dn > len(self.ob[q]["visibility"]["all"][nn:]):
                                dn = len(self.ob[q]["visibility"]["all"][nn:])
                            self.ob[q]["visibility"]["all"][nn:nn+dn] = [False] * dn
                        else:
                            del_list.append(ki)

                    self.random_index_list = numpy.delete(self.random_index_list, del_list)



    def export(self):

        tmp = f'sunset: {self.sunset}'
        self.log_msg(tmp)
        tmp = f'sunrise: {self.sunrise}'
        self.log_msg(tmp)
        tmp = f'plan begining: {ephem.Date(self.nightTime[0])}'
        self.log_msg(tmp)

        ZNACZNIK = False
        if "start_at_sunset" in self.flags:
            txt = f'WAIT sunset={str(self.cfg[self.tel]["hsun"])}'
        else:
            txt = f'WAIT ut={str(ephem.Date(self.nightTime[0])).split()[1]}'
        self.plan.append(txt)
        self.log_msg(txt)
        for n, tmp in enumerate(self.nightTime):
            if self.nightFree[n]:
                ZNACZNIK = True
            if self.nightObChange[n]:  # tutaj sprawdzamy ze jest nowy obiekt alokowany
                if ZNACZNIK:   # tutaj wstawiamy wait, jak nie byla przerwa w alokacji
                    ZNACZNIK = False
                    self.oca.date = ephem.Date(self.nightTime[n])
                    txt = "WAIT ut=" + str(self.oca.date).split()[1]
                    self.plan.append(txt)
                    self.log_msg(txt)
                if True:
                    k = int(self.nightObject[n])
                    if d_ver(self.ob[k], ["name", "makro"]):
                            with open(self.ob[k]["makro"], "r", encoding="utf-8") as makro_file:
                                for line in makro_file:
                                    if line.strip():
                                        self.plan.append(line.strip())
                                        self.log_msg(txt)
                    else:

                        #par2remove = list(self.cfg["params_2_remove"])
                        schema = ObsValidator.load_schema("base_rules")["OBJECT"]
                        allowed = set(schema["allowed"])
                        filtered = {
                            k: v for k, v in self.ob[k].items()
                            if k in allowed
                        }

                        txt = ObsValidator.convert_from_obdict(filtered)

                        self.plan.append(txt)
                        self.log_msg(txt)


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
                txt = f"start_makro added"
                if self.loud: sys.stderr.write(txt)
                self.log_msg(txt)
            else:
                txt = f"NO start_makro FILE"
                if self.loud: sys.stderr.write(txt)
                self.log_msg(txt)



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
                txt = f"end_makro added"
                if self.loud: sys.stderr.write(txt)
                self.log_msg(txt)
            else:
                txt = f"NO end_makro FILE"
                if self.loud: sys.stderr.write(txt)
                self.log_msg(txt)


    def SavePlan(self):

            if self.save_plan:
                # tworzenie nazwy pliku z planem
                try:
                    if not os.path.exists(self.cfg["plan_catalog"]):
                        os.makedirs(self.cfg["plan_catalog"])  # tworzy wszystkie brakujące katalogi w ścieżce
                except PermissionError:
                    self.log_msg(f'No permission to create directory: {self.cfg["plan_catalog"]}')
                    #raise PermissionError(f"No permission to create directory: {self.cfg["plan_catalog"]}")

                sunset_date = str(ephem.Date(self.sunset)).split()[0].split("/")
                yyyymmdd = f"{sunset_date[0]}{sunset_date[1].zfill(2)}{sunset_date[2].zfill(2)}"

                # sprawdzamy czy sa juz takie pliki
                files = os.listdir(self.cfg["plan_catalog"])
                nums = []
                prefix = f"{self.tel}_{yyyymmdd}_"
                suffix = ".txt"
                for f in files:
                    if f.startswith(prefix) and f.endswith(suffix):
                        try:
                            n = int(f[len(prefix):-len(suffix)])
                            nums.append(n)
                        except ValueError:
                            pass
                nr = max(nums) + 1 if nums else 1
                filename = f"{self.tel}_{yyyymmdd}_{nr}.txt"
                logname = f"{self.tel}_{yyyymmdd}_{nr}_log.txt"

                self.plan_filename = os.path.join(self.cfg["plan_catalog"], filename)
                self.log_filename = os.path.join(self.cfg["plan_catalog"], logname)

                # zapisujemy
                try:
                    with open(self.plan_filename, "w") as plan_file:
                        for txt in tqdm(self.plan, desc=f'Saving plan to file {self.plan_filename}'):
                            plan_file.write(txt+"\n")
                        self.plan_saved = True

                except Exception as e:
                    txt = f"ERROR saving plan: {e}"
                    self.log_msg(txt)

                try:
                    with open(self.log_filename, "w") as log_file:
                        log_file.write(self.msg)

                except Exception as e:
                    txt = f"ERROR saving log: {e}"
                    self.log_msg(txt)


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

def merge_schemas(base: dict, extra: dict) -> dict:
    merged = base.copy()

    # properties
    merged.setdefault("properties", {})
    merged["properties"].update(extra.get("properties", {}))

    # required (jeśli masz)
    if "required" in base or "required" in extra:
        merged["required"] = list(set(base.get("required", []) + extra.get("required", [])))

    # additionalProperties – ostrożnie (tu przykład: AND)
    if "additionalProperties" in base or "additionalProperties" in extra:
        merged["additionalProperties"] = (
            base.get("additionalProperties", True)
            and extra.get("additionalProperties", True)
        )

    return merged