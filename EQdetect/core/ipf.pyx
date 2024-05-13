import numpy as np
import math
from numpy.random import random
from scipy.linalg import cholesky

def RPF(particle1, varx, vary, varz, dmin, dmax):
    particle = particle1[:,:3]
    N = len(particle)
    #for i in range(N):
        #particle[i][0] = (particle[i][0] - eq_cent[0])*111.320
        #particle[i][1] = (particle[i][1] - eq_cent[1])*110.574
        #particle[i][2] = particle[i][2]/111
    #S = 1./(len(wi)) * (particle * wi[:, np.newaxis]).T.dot(particle)
    #L = cholesky(S,lower=False, overwrite_a=True)
    m = 1
    h=(8. * 3./4. * 5. * 8. * np.sqrt(np.pi) / N)**0.2
    if varx <0.3:
        varx = 0.3
    if vary <0.3:
        vary = 0.3
    if varz <0.3:
        varz = 0.3
    for i in range(N):
        #epsilon[i]=h * L @ np.random.beta(m)
        #particle[i] += h * (L @ epsilon[i]).T
        #particle[i] = []
        e = np.sqrt(np.random.beta(1.5,2.0,3))
        x = np.random.normal(0.0,1.0,3).reshape(1,-1)
        x = x / np.linalg.norm(x)
        epsilon = np.multiply(e,x)[0]
        while particle[i][2] < dmin or particle[i][2] > dmax:
            e = np.sqrt(np.random.beta(1.5,2.0,3))
            x = np.random.normal(0.0,1.0,3).reshape(1,-1)
            x = x / np.linalg.norm(x)
            epsilon = np.multiply(e,x)[0]
            particle[i][2] = particle[i][2]+h*epsilon[2]*varz
        particle1[i][0] = particle[i][0]+h*epsilon[1]*varx
        particle1[i][1] = particle[i][1]+h*epsilon[0]*vary
        particle1[i][2] = particle[i][2]
    return particle1
    
def resample1(weights, particles):
    N = len(weights)
    indexes = np.zeros(N, 'i')

    # take int(N*w) copies of each weight, which ensures particles with the
    # same weight are drawn uniformly
    num_copies = (np.floor(N*np.asarray(weights))).astype(int)
    k = 0
    for i in range(N):
        for _ in range(num_copies[i]): # make n copies
            indexes[k] = i
            k += 1

    # use multinormal resample on the residual to fill up the rest. This
    # maximizes the variance of the samples
    residual = weights - num_copies     # get fractional part
    residual /= sum(residual)           # normalize
    cumulative_sum = np.cumsum(residual)
    cumulative_sum[-1] = 1. # avoid round-off errors: ensures sum is exactly one
    indexes[k:N] = np.searchsorted(cumulative_sum, random(N-k))
    new_particles = particles[indexes]

    return new_particles

def resample(weights, particles):
    #normalized_weights = [w / sum(weights) for w in weights]
    Q = np.cumsum(weights)
    new_particles = []
    N = len(weights)
    m = 0
    for i in range(N):
        u0 = np.random.uniform(1.0 / (N*150), 1.0 / N, 1)[0]
        u = u0 + float(i) / N
        while Q[m] < u:
            m += 1  # no need to reset m, u always increases
        # Add state sample (weight, state)
        new_particles.append(particles[m])
    return new_particles

def arrp(dis, hypo):
    if hypo[2] < 24.0:
        v0 = 4.55220 # initial velocity in shallow layer 
        vg = 0.12040 # gradient velocity in shallow layer
    else:
        v0 = 7.12720 # initial velocity in deep layer
        vg = 0.00610 # gradient velocity in deep layer
    xc = (dis * dis - 2. * v0 / vg * hypo[2] - hypo[2] * hypo[2]) / (2. * dis)
    zc = -1. * v0 / vg
    ang1 = math.atan((hypo[2] - zc) / xc)
    if ang1 < 0.0:
        ang1 = ang1 + np.pi
    ang1 = np.pi - ang1
    ang2 = math.atan(-1. * zc / (dis - xc))
    ptime = (-1. / vg) * math.log(abs(math.tan(ang2 / 2.) / math.tan(ang1 / 2.))) # P-wave travel time
    return ptime

def delaz(lat1, lon1, lat2, lon2):
    # convert degree to kilometer
    dlat = (lat2 - lat1) * 111.12
    dlon = (lon2 - lon1) * 111.12 * math.cos(lat1 * np.pi / 180.0)
    return math.sqrt(dlat * dlat + dlon * dlon)

def phys(x):
    res = np.exp(-0.5*(x**2))/np.sqrt(2*np.pi)
    return res

def gauss_lklhood(sigmp, sigma, pd, pickerr, dist, hypo):
    fp = phys(pickerr/sigmp)/sigmp
    mua = np.exp(0.72*hypo[3]-1.2*np.log(dist)-0.0005*dist+0.005*hypo[2] - 0.46)
    fa = phys((pd*100-mua)/sigma)/sigma
    return fp, fa

def cal_ipf(hypo, station, nottrigsta, timenow):
    m = []
    allpt=[]
    allpd = []
    allot = []
    alldis = []
    allpo = []
    for pick in station:
        dist = delaz(pick.latitude, pick.longitude, hypo[1], hypo[0])
        ptime = arrp(dist, hypo)
        allpt.append(ptime)
        allpd.append(pick.pd)
        allot.append(pick.picktime-ptime)
        allpo.append(pick.picktime)
        alldis.append(dist)
        m.append((np.log(pick.pd*100)+1.2*np.log(dist)+0.0005*dist-0.005*hypo[2] + 0.46)/0.72)
    if len(m)==1:
        mag =  m[0]
        ot = allot[0]
    else:
        mag = round(np.median(m),2)
        ot = round(np.median(allot),2)
    hypo[3] = mag
    hypo[4] = ot
    weight = 1.0
    fp_all = 1.0
    fa_all = 1.0
    for i in range(len(allpt)):
        perr = allpo[i]-(ot+allpt[i])
        fp, fa = gauss_lklhood(1.5, 0.07, allpd[i], perr, alldis[i], hypo)
        if fp<0.0004:
            fp = 0.0004
        if fa<0.0004:
            fa = 0.0004
        fp_all *= fp
        fa_all *= fa
    for sta in nottrigsta:
        dist = delaz(sta[2], sta[1], hypo[1], hypo[0])
        if round(dist,3) == 0.0:
            ptime = 0.0 
        else:
            ptime = arrp(dist, hypo)
        perr = timenow - (ot+ptime)
        if perr >= 0:
            fp = 1
        else:
            fp, fa = gauss_lklhood(1.5, 0.07, 1.0, perr, dist, hypo)
        if fp<0.004:
            fp = 0.004
        fp_all *= fp

    weight = fp_all*fa_all
    hypo[5] = weight
    hypo[6] = fp_all
    hypo[7] = fa_all
    return weight, hypo
    




