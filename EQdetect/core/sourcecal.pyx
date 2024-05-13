import math
import inspect
import numpy as np
import traceback
from datetime import datetime
import multiprocessing
from .ipf import cal_ipf, arrp, delaz, resample, gauss_lklhood, RPF, resample1

def gridhypo(point, dx, dy, num_n=16):
    x_min, x_max = point[0] - dx, point[0] + dx
    y_min, y_max = point[1] - dy, point[1] + dy
    z_min, z_max = 0, 100
    point_grd = []
    z_values = np.random.uniform(z_min, z_max, num_n**3)
    x_values = np.arange(x_min, x_max, dx*2/num_n)
    y_values = np.arange(y_min, y_max, dy*2/num_n)
    z_values1 = np.arange(z_min, z_max, z_max/num_n)
    x_grid, y_grid, z_grid = np.meshgrid(x_values, y_values, z_values1, indexing='ij')
    for x, y, z in zip(x_grid.flatten(), y_grid.flatten(),z_values):
        point_grd.append([round(x,5), round(y,5), z, 0.0, 0.0, 0.0, 0.0, 0.0])
    return np.array(point_grd)
    
def pending_eq(eqdata, vorcel):
    eqdata.cal_hypo(vorcel)
    if eqdata.optsrc[3] > 10.0 or eqdata.optsrc[3] <= 0.0 or abs(eqdata.rmse)<0.001:
        return False, eqdata
    eqdata.endtime = datetime.utcnow()
    eqdata.update = eqdata.update + 1
    print(
        f"Pending Earthquake {eqdata.first.sta} No {eqdata.noeq} with paramater {eqdata.starttime}\n")
    file = open("output.txt", "a")
    file.write(
        f"Pending Earthquake No {eqdata.noeq} with paramater {eqdata.starttime}\n")
    file.write(
        f"Longitude : {round(eqdata.optsrc[0],3)}+-{round(eqdata.optsrc[5],3)}\nLatitude : {round(eqdata.optsrc[1],3)}+-{round(eqdata.optsrc[6],3)}\nDepth : {round(eqdata.optsrc[2],5)}\nOT : {round(eqdata.optsrc[4],2)}\nMag : {round(eqdata.optsrc[3],2)}\n")
    file.write(
        f"First Trigger : {eqdata.first.sta} {eqdata.triggroup}\n")
    file.write(
        f"First Trigger prior : {eqdata.first.weight} {eqdata.rmse}\n\n")
    file.close()
    return True, eqdata

class Phase:
    def __init__(self, pick):
        #pick = pick.split(" ")
        if len(pick) < 13:
            raise ValueError("Length of pick message is not appropiate")
        self.sta = pick[0]
        self.comp = pick[1]
        self.longitude = float(pick[4])
        self.latitude = float(pick[5])
        self.pa = float(pick[6])
        self.pv = float(pick[7])
        self.pd = float(pick[8])
        self.picktime = float(pick[10])
        self.weight = float(pick[11])
        self.telflag = float(pick[12])
        self.upd_sec = float(pick[13])
        self.reserr = 999
        self.status = 0

    def set_first(self, voronoi):
        #if inspect.isclass(voronoi) is False:
        #    raise ValueError("The voronoi input is not Class")
        #try:
        #    idx = voronoi.code.index(self.sta)
        #except:
        #    print(f"{self.sta} is not in list station database")
        #    return False
        #self.triggroup = voronoi.sta_trgg[idx]
        #self.estgroup = voronoi.sta_estg[idx]
        #self.init_eq = voronoi.init_eq[idx]
        #return True
        try:
            self.triggroup, self.estgroup, self.init_eq = voronoi.get_vor(self.sta)
            return True
        except:
            print(f"{self.sta} is not in list station database")
            return False
    
    def calerror(self, hypo):
        dist = delaz(self.latitude, self.longitude, hypo[1], hypo[0])
        ptime = arrp(dist, hypo)
        self.reserr = self.picktime - hypo[4] - ptime

class EQsrc:
    """
    status = 0 is pending earthquake
    status = -1 is cancel earthquake
    status = 1 is ongoing earthquke 
    status = -3 is change first pick
    """
    def __init__(self, ftrig, cfg, noeq):
        self.first = ftrig
        self.status = 0
        self.initsrc = ftrig.init_eq
        self.triggroup = ftrig.triggroup
        self.estgroup = ftrig.estgroup
        self.starttime = datetime.utcnow()
        self.endtime = datetime.utcnow()
        self.rmphase = []
        self.phase = [ftrig]
        self.code = [ftrig.sta]
        self.numtrg = 0
        self.optsrc = [self.initsrc[0], self.initsrc[1], 10.0, 0.0, 0.0, 0.0]
        self.update = 0
        self.cfg = cfg
        self.var = []
        self.report=False
        self.noeq = noeq
        self.need_resample = False
        self.rmse = 0

    def check_trgg(self, pick):
        pick.calerror(self.optsrc)
        if abs(pick.reserr) <= self.cfg.trig_err and (pick.sta in self.triggroup):
            if pick.sta in self.code:
                idx = self.code.index(pick.sta)
                picking = self.phase[idx]
                if picking.picktime == pick.picktime:
                    self.phase[idx] = pick
                else:
                    return False
            else:
                self.phase.append(pick)
                self.code.append(pick.sta)
                self.numtrg += 1
        elif abs(pick.reserr) <= self.cfg.est_trig_err and (pick.sta in self.estgroup):
            if pick.sta in self.code:
                idx = self.code.index(pick.sta)
                picking = self.phase[idx]
                if picking.picktime == pick.picktime:
                    self.phase[idx] = pick
                else:
                    return False
            else:
                self.phase.append(pick)
                self.code.append(pick.sta)            
        else:
            return False
        return True

    def check_estg(self, pick):
        if pick.sta in self.estgroup:
            if pick.sta in self.code:
                idx = self.code.index(pick.sta)
                picking = self.phase[idx]
                if picking.picktime == pick.picktime:
                    self.phase[idx] = pick
                else:
                    return False
            else:
                pick.calerror(self.optsrc)
                if abs(pick.reserr) <= self.cfg.trig_err and self.update < 1:
                    self.phase.append(pick)
                    self.code.append(pick.sta)
                elif abs(pick.reserr) <= self.cfg.est_trig_err and self.update > 1:
                    self.phase.append(pick)
                    self.code.append(pick.sta)
                else:
                    return False
        else:
            pick.calerror(self.optsrc)
            if abs(pick.reserr) <= self.cfg.est_trig_err:
                return True
            if (abs(pick.reserr) <= self.cfg.est_trig_err + 1.5) and (self.update>6):
                return True
            return False
        return True
    
    def eq_stat(self, vorcel):
        if self.numtrg >= 2 and self.status == 0:
            self.status = 1
        if self.numtrg >= 1 and self.status == 0 and len(self.phase)>=4:
            self.status = 1
        dttime = datetime.utcnow()-self.starttime
        if self.status == 0:
            if dttime.total_seconds() >= 80:
                self.status = -1
        if self.status >= 1:
            if (self.optsrc[3]>= 6 and dttime.total_seconds() > 8*60) or (self.optsrc[3]< 6 and dttime.total_seconds() > 5*60): 
                self.status = -1
        if self.status == 1:
            if self.update >= 50:
                self.status = 2
                return
            fp_rm = 0
            nsta = len(self.phase)
            #for i in range(len(self.phase)):
            i = 0
            while i < nsta:
                if i > len(self.phase) - 1:
                    break
                self.phase[i].calerror(self.optsrc)
                if abs(self.phase[i].reserr) > self.cfg.est_trig_err:
                    self.rmphase.append(self.phase[i])
                    if self.phase[i].sta == self.first.sta:
                        self.phase[i].status = -1
                        print(f"Remove first pick in estimation with code {self.phase[i].sta} {self.phase[i].reserr}")
                        fp_rm = 1
                    else:
                        print(f"Remove pick in estimation with code {self.phase[i].sta} {self.phase[i].reserr}")
                    if self.phase[i].sta in self.code:
                        idx = self.code.index(self.phase[i].sta)
                        self.code.pop(idx)
                    self.phase.pop(i)
                    if len(self.phase)==0 or i==nsta-1:
                        i = nsta + 1
                        continue
                else:
                    i += 1
                    if i==nsta-1:
                        i = nsta + 1
            if len(self.phase)<2:
                self.status=-1
                for pick in self.phase:
                    self.rmphase.append(pick)
                return
            else:
                if fp_rm == 1:
                    timepick = 99999999999999999
                    for pick in self.phase:
                        if pick.picktime < timepick:
                            self.first = pick
                            timepick = pick.picktime
                    self.first.set_first(vorcel)
                    self.initsrc = self.first.init_eq
                    self.triggroup = self.first.triggroup
                    self.estgroup = self.first.estgroup
                    self.numtrg = len(self.phase)
                    self.optsrc = [self.initsrc[0], self.initsrc[1], 10.0, 0.0, 0.0, 0.0]
                    self.status = 1
                return
            if self.update > 10 and len(self.phase)<3:
                self.status=-1
                for pick in self.phase:
                    if pick.sta != self.first.sta: 
                        self.rmphase.append(pick)

    def cal_hypo(self,vorcel):
        eq_cent = [self.optsrc[0], self.optsrc[1]]
        if not self.need_resample:
            particle = gridhypo(eq_cent, self.initsrc[2], self.initsrc[3],13)
        else:
            particle = self.particle
        allwi = []
        nottrigsta = []
        #find not triggered station
        #if self.update > 0:
        for sta in self.estgroup:
            if sta in self.code:
                continue
            df = vorcel.list_station[vorcel.list_station["Kode"] == sta]
            if df.empty:
                continue
            lats, longs = df["Lat"].values[0], df["Long"].values[0]
            nottrigsta.append([sta, longs, lats])

        for i in range(len(particle)):
            wi, pf = cal_ipf(particle[i], self.phase, nottrigsta, self.endtime.timestamp())
            particle[i] = pf
            allwi.append(wi)
        allwi /= np.sum(allwi)
        allwi = np.array(allwi) 
        x0 = np.sum(np.array(allwi)*np.array(particle)[:,0])
        y0 = np.sum(np.array(allwi)*np.array(particle)[:,1])
        z0 = np.sum(np.array(allwi)*np.array(particle)[:,2])
        varx0 = np.sqrt(np.sum((np.array(allwi)*np.array(particle)[:,0]**2))-x0**2)
        vary0 = np.sqrt(np.sum((np.array(allwi)*np.array(particle)[:,1]**2))-y0**2)
        varz0 = np.sqrt(np.sum((np.array(allwi)*np.array(particle)[:,2]**2))-z0**2)
        neff = 1. / np.sum(np.square(allwi))
        dmin = min(particle[:,2])
        dmax = max(particle[:,2])
        #Resampling
        if self.update > 0 and neff<len(allwi)/2:
            new_particle = resample(allwi, particle)
            self.need_resample = True
            if isinstance(new_particle, np.ndarray):
                new_particle = new_particle.tolist()
            #RPF
            try:
                allwi = np.ones(len(particle))*(1.0/len(particle))
                allwi /= np.sum(allwi)
                particle = RPF(np.array(new_particle), varx0, vary0, varz0, dmin, dmax)

                #recalculate the weight
                allwi = []
                for i in range(len(particle)):
                    wi, pf = cal_ipf(particle[i], self.phase, nottrigsta, self.endtime.timestamp())
                    particle[i] = pf
                    allwi.append(wi)
                allwi /= np.sum(allwi)
            except Exception as e:
                print(e, str(traceback.extract_stack()[-1][1]))
                print("RPF Error",neff)
                self.need_resample = False
        else:
            self.need_resample = False
            
               
        #mag0 = np.sum(np.array(allwi)*np.array(particle)[:,3])
        #ot0 = np.sum(np.array(allwi)*np.array(particle)[:,4])
        id_max = np.argmax(particle,axis=0)[5]
        if self.need_resample is False and self.update < 3:
            [x0, y0, z0, mag0, ot0] = particle[id_max][:5]
        else:
            #calculate optimal hypocenter 
            x0 = np.sum(np.array(allwi)*np.array(particle)[:,0])
            y0 = np.sum(np.array(allwi)*np.array(particle)[:,1])
        z0 = np.sum(np.array(allwi)*np.array(particle)[:,2])

        if self.update == 0:
            z0 = 10.0
        allot = []
        magall = []
        for pick in self.phase:
            dist = delaz(pick.latitude, pick.longitude, y0, x0)
            magall.append(5.067 + 1.281 * np.log10(pick.pd) + 1.760 * np.log10(dist))
            #magall.append((np.log(pick.pd*100)+1.2*np.log(dist)+0.0005*dist-0.005*z0 + 0.46)/0.72)
            ptime = arrp(dist, [x0, y0, z0])
            allot.append(pick.picktime-ptime)
        ot0 = round(np.median(allot),2)
        mag0 = np.median(magall)

        varx0 = np.sqrt(np.sum((np.array(allwi)*np.array(particle)[:,0]**2))-x0**2)
        vary0 = np.sqrt(np.sum((np.array(allwi)*np.array(particle)[:,1]**2))-y0**2)
        varz0 = np.sqrt(np.sum((np.array(allwi)*np.array(particle)[:,2]**2))-z0**2)
        varot0 = np.sqrt(np.sum((np.array(allwi)*np.array(particle)[:,4]**2))-ot0**2)
        #print(varx0, vary0, varz0, varot0, x0, y0, z0, ot0, np.sum((np.array(allwi)*np.array(particle)[:,0]**2))-x0**2)
        rmse = []
        for pick in self.phase:
            pick.calerror([x0, y0, z0, mag0, ot0, varx0, vary0, varz0, neff])
            rmse.append(pick.reserr)
        rmse = np.array(rmse)
        if len(rmse) == 1:
            self.rmse = rmse[0]
        else:
            self.rmse = np.sqrt(np.mean(rmse**2))
        
        if delaz(self.optsrc[1],self.optsrc[0],y0,x0) > 50 or max([varx0, vary0]) > 0.15 or len(self.phase)<3:
            self.need_resample = False

        self.optsrc = [x0, y0, z0, mag0, ot0, varx0, vary0, varz0, neff]
        wi, pf = cal_ipf(self.optsrc[:8], self.phase, nottrigsta, self.endtime.timestamp())
        self.first.weight = wi
        self.particle = particle
        self.var.append([x0, y0])  