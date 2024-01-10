#!/usr/bin/env python3

#=============================
#        Marek Gorski
#=============================

# ./tpg.py zb08 ( or v16)  2019/08/20 [--nv --ph -o]  >| plan.txt
# --ph : reject object based on phase caverage"
# --loud : loud to sys.sterr
# -o : overload objects

import sys,os,numpy,ephem

class main:
  def __init__(self):

      self.cwd = os.getcwd()
      self.script_dir = os.path.dirname(os.path.abspath(__file__))

      #--------------- MAIN SCRIPT --------------
      self.config()
      self.getargs()                  # reads script arguments
      self.Initiate()                 # initiate script variables
      self.LoadObjects()              # loads objects to observe
      self.MakeTime()                 # prepares time variable
      self.CalcObject()               # calculates visibilities for all objects
      self.Waga()                     # calculate propability for randomization order of objects to observe
      self.RandomizeList()            # randomize listo of objects to observev
      self.MaskStartEnd()                  # mask visibilite base on t_start and t_end
      # self.MaskPhase()                # reject objects based on phase covarage
      # self.MakeCycle()                # changes prioryty based on cycle
      self.allocate()                 #   SCHEDULING OBJECTS
      self.export()                   # export plan

  def config(self):

      if not os.path.isfile(self.script_dir+'/tpg.cfg'):
          sys.stderr.write('WARNING - No config file tpg.cfg:\n')
          sys.exit()
      else:
         plik=open(self.script_dir+"/tpg.cfg",'r')
         for line in plik: 
             if line:  # no empty lines
                if line[0]!="#":
                   if "zb08_FILE=" in line:
                       self.zb08_file=line.split("=")[1].strip()

                   # if "IRISFILE=" in line:
                   #    self.irisfile=line.split("=")[1].strip()
                   # if "V16FILE=" in line:
                   #    self.v16file=line.split("=")[1].strip()
                   # if "BESTFILE=" in line:
                   #    self.bestfile=line.split("=")[1].strip()
                   # if "OBSDATA=" in line:
                   #    self.obsdata_dir=line.split("=")[1].strip()
                   # if "LOGDIR=" in line:
                   #    self.logdir=line.split("=")[1].strip()
                   # if "TABALL=" in line:
                   #    self.taball_path=line.split("=")[1].strip()
                   # if "OBSNUM=" in line:
                   #    self.obs_num=int(line.split("=")[1].strip())
                   # if "PHCOV=" in line:
                   #    self.phase_coverage=float(line.split("=")[1].strip())
                   # if "PHNUM=" in line:
                   #    self.phase_num=int(line.split("=")[1].strip())
                   if "HMAX=" in line:
                      self.h_max=line.split("=")[1].strip()
                   if "HMIN=" in line:
                      self.h_min=line.split("=")[1].strip()
                   if "HSUN=" in line:
                      self.h_sun=line.split("=")[1].strip()

       
  def getargs(self):
      if len(sys.argv)<3:
         print("wrong format: ./tpg.py telescope (zb08 or jk15) date (yyyy/mm/dd) [opt: master_file.txt] \n")
         print("example: ./tpg.py zb08 2023/08/12")
         sys.exit()
      else:
         self.flags=[]
         if sys.argv[1]=="iris": 
            self.Obfile=self.irisfile
            self.flags.append("iris")
            self.overhead=100.
         elif sys.argv[1]=="zb08":
            self.Obfile=self.zb08_file
            self.flags.append("zb08")
            self.overhead=20.

         self.date=sys.argv[2]+" 16:00:00"
         if "--ph" in sys.argv: self.flags.append("ph")
         if "--loud" in sys.argv: self.flags.append("loud")
         if "-o" in sys.argv: self.flags.append("overload")                   
                   

  def Initiate(self):

      #plan                      # 1 minute resolution
      self.nightTime=[]          # time in ephem.Date format
      self.nightObject=[]        # index of objects to observe
      self.nightFree=[]          # Bool marker of allocation 
      
      #objects to obserrve                
      self.name=[]        
      self.ra=[]
      self.dec=[]
      self.seq=[]
      self.comment=[]
      self.slotTime=[]           # minutes of allocation for object
      self.spec=[]               # specification for NV to observe
      self.tstart=[]             # 
      self.tend=[]               #
      self.obflag=[]             #
      


      #additional
      self.random_index_list=[]     # list with random roder of object indices to oberve
      self.pliki_hjd=[]             # files with phase coverage for objects
      

      self.probability=[]           # propability for object list randomization         
      self.visibility=[]            # visibility of each object (self.nightTime indicies)
      self.h=[]                     # altitude over horizon
      self.priority=[]              # priority for observations
      self.cycle=[]

  def LoadObjects(self):  # Load Object file

      if not os.path.isfile(self.Obfile):
          sys.stderr.write('No Objects Catalog file!\n')
          sys.stderr.write('Check tpg.cfg file\n')
          sys.exit()

      plik = open(self.Obfile, 'r')
      for line in plik:
          line = line.strip()
          if line:  # no empty lines
              if line.strip()[0] != "#":
                  flag = []
                  seq = "1/Ic/30"
                  priority = 0
                  slotTime = 0
                  comment=""

                  overhed = self.overhead
                  # cycle=1
                  # total=False
                  # tt=0
                  ts=False
                  te=False
                  if "comment=" in line:
                      tmp = line.split("comment=")[1]
                      if "\"" in line:
                          comment = tmp.split("\"")[1]
                      if "\'" in line:
                          comment = tmp.split("\'")[1]

                  if "seq=" in line:
                      seq=line.split("seq=")[1].split()[0]

                  if "priority=" in line:
                      priority=line.split("priority=")[1].split()[0]

                  if "t_start=" in line:
                      ts = line.split("t_start=")[1].split()[0]


                  if "t_end=" in line:
                      te = line.split("t_end=")[1].split()[0]

                  for i, tmp in enumerate(line.split()):

                      if i == 0:
                          name = line.split()[0]
                      elif i == 1:
                          ra = line.split()[1]
                      elif i == 2:
                          dec = line.split()[2]

                      #elif i > 2:

                          # if "flag=ignore_phs" in line.split()[i]:
                          #    flag.append("ignore_phs")
                          #
                          # if "flag=ignore_vis" in line.split()[i]:
                          #    flag.append("ignore_vis")
                          #
                          # if "nv_overhead=" in line.split()[i]:
                          #    overhed=float((line.split()[i].split("=")[1]))
                          # if "nv_total_time=" in line.split()[i]:
                          #    total=float((line.split()[i].split("=")[1]))
                          # if "cycle=" in line.split()[i]:
                          #    cycle=int(line.split()[i].split("=")[1])

                  slotTime = 0
                  k = 1
                  if "x(" in seq and ")" in seq:
                      k = int(seq.split("x")[0])
                      seq = seq.split("x(")[1].split(")")[0]
                  for x_seq in seq.split(","):
                      slotTime = slotTime + (float(x_seq.split("/")[0]) * (float(x_seq.split("/")[2]) + float(overhed)))/60.
                  slotTime=k*slotTime
                  # if len(time)>0 and len(time)==len(ndit):
                  #    for a,b in zip(time,ndit):
                  #        if "overload" in self.flags: a=float(a)/2.
                  #        tt=tt+(float(a)+overhed)*float(b)
                  # slotTime=int(float(tt)/60)
                  # if total: slotTime=int(float(total)/60)
                  # slotTime=slotTime+1

                  # if ts:
                  #    if "d" in ts:
                  #       tmp = ts
                  #       ts = tmp.split("d")[0]+" "+tmp.split("d")[1]
                  #       ts=float(ephem.Date(ts))
                  #    else:
                  #       ts=self.date.split()[0]+" "+ts
                  #       ts=ephem.Date(ts)
                  #       if ts<ephem.Date(self.date.split()[0]+" 16:00:00") :ts=ts+1
                  # if te:
                  #    if "d" in te:
                  #       tmp = te
                  #       te = tmp.split("d")[0]+" "+tmp.split("d")[1]
                  #       te=float(ephem.Date(te))
                  #    else:
                  #       te=self.date.split()[0]+" "+te
                  #       te=ephem.Date(te)
                  #       if te<ephem.Date(self.date.split()[0]+" 16:00:00") :te=te+1

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
                  # self.obflag.append(flag)
                  # self.cycle.append(cycle)

  def MakeTime(self):  # Make observing time
      self.oca = ephem.Observer()
      self.oca.date = self.date
      self.oca.lon = str(-70.2)
      self.oca.lat = str(-24.58)
      self.oca.elev = 2800
      self.oca.horizon = self.h_sun  # altitude of sun to start and stop observations

      self.sunset = self.oca.next_setting(ephem.Sun(), use_center=True)
      self.oca.date = self.sunset
      self.sunrise = self.oca.next_rising(ephem.Sun(), use_center=True)

      t = self.sunset
      while t < self.sunrise:
          t = t + 1 * ephem.minute  # time resolution : 1 minute
          self.nightTime.append(t)
          self.nightObject.append(-100)
          self.nightFree.append(True)

  def CalcObject(self):       # Colculate visibility for all objects
      for n,tmp in enumerate(self.name):
          s = ephem.FixedBody()
          s._ra = self.ra[n]
          s._dec = self.dec[n]
          self.oca.date = self.sunset
          self.oca.horizon=self.h_sun
          vis=[]
          hT=[]
          s.compute(self.oca)
          for i,tmp in enumerate(self.nightTime):
              self.oca.date=self.nightTime[i]
              s.compute(self.oca)
              h=float(str(s.alt).split(":")[0])+float(str(s.alt).split(":")[1])/60
              hT.append(h)
              if float(h)>float(self.h_min) and float(h)<float(self.h_max): T=True
              else: T=False
              #if "ignore_vis" in self.obflag[n]: T=True
              vis.append(T)
          self.visibility.append(vis)
          self.h.append(hT)

  def Waga(self):     #
      for n,tmp in enumerate(self.name):
          maska=numpy.array(self.visibility[n])
          maska=maska[maska]
          allminutes=len(self.nightTime)
          if len(maska)==0: waga=0
          else: waga=1.-float(len(maska))/float(allminutes)
          if waga<0.3 and len(maska)!=0 : waga=0.3
          self.probability.append(waga)

  def RandomizeList(self):          # randomize listo of objects to observev
      mozliwosci=numpy.array([i for i, _ in enumerate(self.name)])
      probability=numpy.array(self.probability)
      maska=(probability!=0)
      mozliwosci=mozliwosci[maska]
      probability=probability[maska]
      for n,tmp in enumerate(mozliwosci):
          suma=sum(probability)
          probability=[x/suma for x in probability]
      self.random_index_list=numpy.random.choice(mozliwosci, len(mozliwosci),replace=False, p=probability)

  def allocate(self):
      for nn, tmp in enumerate(self.nightTime):
          if self.nightFree[nn]:

              OK = False
              done = -1
              priormax = -1

              for nk, tmp in enumerate(self.random_index_list):
                  k = self.random_index_list[nk]

                  visibility = self.visibility[k]
                  slots = int(float((self.slotTime[k])) / 1.)
                  bounds = (nn, nn + slots)

                  if all(self.nightFree[bounds[0]:bounds[1] + 1]) and all(visibility[bounds[0]:bounds[1] + 1]) and \
                          int(self.priority[k]) > -50:
                      if int(self.priority[k]) > priormax:
                          priormax = int(self.priority[k])
                          bfin = bounds
                          kfin = k
                          OK = True
                          done = nk

              if OK:
                  self.nightFree[bfin[0]:bfin[1] + 1] = [False] * ((bfin[1] + 1) - bfin[0])
                  self.nightObject[bfin[0]:bfin[1] + 1] = [int(kfin)] * ((bfin[1] + 1) - bfin[0])
                  self.random_index_list = numpy.delete(self.random_index_list, done)

  def export(self):
      if True:
         prev=-100
         ZNACZNIK=False
         for n,tmp in enumerate(self.nightTime):
             if int(self.nightObject[n]) != prev:
                if ZNACZNIK:
                   ZNACZNIK=False
                   empty_time =  float(self.nightTime[n]-t0)*24.*60.*60.
                   self.oca.date = ephem.Date(self.nightTime[n])
                   txt = "WAIT ut="+str(self.oca.date).split()[1]
                   print(txt)
                prev = int(self.nightObject[n])
                if prev>-50:
                   k=int(self.nightObject[n])
                   #print(self.spec[k])
                   txt = f"OBJECT {self.name[k]} {self.ra[k]} {self.dec[k]} seq={self.seq[k]}"
                   if len(self.comment[k]) > 0: txt = txt + f" comment=\"{self.comment[k]}\""
                   print(txt)

                else:
                  t0=self.nightTime[n]
                  ZNACZNIK=True

  def MaskStartEnd(self):  # Masks visibility based on t_start t_end
      for n, tmp in enumerate(self.name):
          ZNACZNIK = False

          ts = self.tstart[n]
          if ts:
              if "T" in ts:
                  ts = ts.replace("T", " ")
              if "/" not in ts:
                  date = self.date.split()[0]
                  if int(ts.split(":")[0])<16:
                      date = ephem.Date(date)+1
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
                  if int(ts.split(":")[0])<16:
                      date = ephem.Date(date)+1
                      date = str(ephem.Date(date)).split()[0]
                  te = date + " " + te
              te = ephem.Date(te)
              time_e = te

              if time_s or time_e:
                  for i, tmp in enumerate(self.nightTime):
                      if float(self.nightTime[i]) > float(time_e) and self.visibility[n][i] and time_e:
                          self.visibility[n][i] = False
                          ZNACZNIK = True
              if "loud" in self.flags and ZNACZNIK:
                  txt = self.name[n] + " masked: t_start or t_end"
                  sys.stderr.write(txt + "\n")

  def MaskPhase(self):      # mask phase covarage
      if "ph" in self.flags:         

         lista_plikow=os.listdir(self.obsdata_dir)
         for zmienna in lista_plikow:
             pwd=os.path.abspath(self.obsdata_dir)
             pelna_sciezka=str(pwd)+"/"+str(zmienna)
             if "iris" in self.flags: case="iris"
             elif "v16" in self.flags: case="v16"
             if ((".hjd") in pelna_sciezka and case in pelna_sciezka):
                self.pliki_hjd.append(pelna_sciezka)  
       
         
         for n,tmp in enumerate(self.name):
             if "ignore_phs" not in self.obflag[n]:
                ZNACZNIK=False
                name = self.name[n].strip()
                P=0
                jd0=0
                ph_t=[]
                for i,curent_file in enumerate(self.pliki_hjd):
                    hjd_file=self.pliki_hjd[i]
                    if len(name)>3 and ( str(name) in str(hjd_file.split("/")[-1])) :
                       f=open(self.taball_path,'r')
                       for line in f:
                           if name.lower() in line.lower(): 
                              try:
                                 P=float(line.split()[3])
                              except: pass
                              try:
                                 jd0=float(line.split()[4])
                              except: pass   
                       if P>0.01:
                          jd_t=[]
                          f_ph = open(curent_file)
                          for linia in f_ph:
                              a=0  
                              try: a=float(linia.split()[0])
                              except: pass
                              jd_t.append(a)
                          ph_t=[(float(x)-float(jd0))/float(P)%1. for x in jd_t]
                          ph_t=numpy.array(ph_t)
                if len(ph_t)>self.obs_num:
                   for t,tmp3 in enumerate(self.nightTime):
                       if self.visibility[n][t]:
                          jd = ephem.julian_date(self.nightTime[t])
                          
                          faza=(float(jd)-float(jd0))/float(P)%1
                          roznica_faz=faza-ph_t
                          mk=roznica_faz<self.phase_coverage
                          roznica_faz=roznica_faz[mk]
                          mk=roznica_faz>-1.*self.phase_coverage
                          roznica_faz=roznica_faz[mk]
                          if len(roznica_faz)>self.phase_num:              
                             self.visibility[n][t]=False 
                             ZNACZNIK=True
                                                       
                if "loud" in self.flags and ZNACZNIK: 
                   txt=self.name[n]+" masked: phase"
                   sys.stderr.write(txt+"\n")                              
                             


  def MakeCycle(self):     
      for n,tmp in enumerate(self.name):
          cycle=int(self.cycle[n])
          if cycle>1:
             i=0
             ZNACZNIK_i=False
             ZNACZNIK_0=False
             while i<cycle+1:
                i=i+1
                dni=i
                name=self.spec[n].split()[0]
                filt=self.spec[n].split()[3]
                tmp = str(ephem.Date(ephem.Date(self.date)-24*dni*ephem.hour)).split()[0]
                filename=tmp.split("/")[0]+tmp.split("/")[1]+tmp.split("/")[2]+".txt"
                if "v16" in self.flags: filename="log_v16_"+filename
                if "iris" in self.flags: filename="log_iris_"+filename  
                if os.path.isfile(self.logdir+filename):
                   plik=open(self.logdir+filename,'r')
                   linie=plik.readlines()
                   log=[x.split()[1] for x in linie if name in x]
                   filt=filt.split(",")
                   ZNACZNIK=True
                   for f in filt:
                       if not f in log: ZNACZNIK=False 
                   if i==cycle and ZNACZNIK: ZNACZNIK_0=True  
                   if i<cycle and ZNACZNIK: ZNACZNIK_i=True
                   
             if not ZNACZNIK_0 and not ZNACZNIK_i: 
                self.priority[n]=self.priority[n]+1
                if "loud" in self.flags: 
                   txt=name+" over cycle: priority+1"
                   sys.stderr.write(txt+"\n")
             if ZNACZNIK_i: 
                self.priority[n]=1
                if "loud" in self.flags: 
                   txt=name+" within cycle: priority=1"
                   sys.stderr.write(txt+"\n")                




####### SCRIPT RUN
    
if __name__ == "__main__":
    m = main()
