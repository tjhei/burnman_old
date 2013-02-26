# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import numpy as np
import bisect
import matplotlib.pyplot as pyplot
from tools import *
import seismic
import scipy.integrate as integrate

# polynomial fit from Watson, Baxter, EPSL, 2007
# pressure: in GPa
# return: temperature in K
def watson_baxter(pressure):
    if (pressure <= 15e9):
        return 1900-1420*pow(0.8,pressure/1e9)
    else:
        return 1680+11.1*pressure/1e9





# geotherm from Brown and Shankland 81
def brown_shankland(pressure):
    depth = seismic.prem_model.depth(pressure)
    return lookup_and_interpolate(table_brown_depth, table_brown_temperature, depth)	

#This integrates dT/dP = gr * T / K_s
#TODO: what is this?
#def self_consistent(pressure, T0, params):
#    minP = 0
#    #upper mantle potential temperature, in K
#    # integrate the adiabatic gradient equation
#    lnT = integrate.quad( lambda x: (params['ref_grueneisen']/(1.e9*bm.bulk_modulus(x, params))), minP, pressure)
#    T = T0*np.exp(lnT[0])
#    return T
    
    


table_brown = read_table("data/brown_81.txt")
table_brown_depth = np.array(table_brown)[:,0]
table_brown_temperature = np.array(table_brown)[:,1]

# test geotherm
if __name__ == "__main__":
    p = np.arange(1.0e9,128.0e9,3e9)
    t1 = [watson_baxter(y) for y in p]
    t2 = [brown_shankland(y) for y in p]
    p1,=pyplot.plot(p,t1,'x--r')
    p2,=pyplot.plot(p,t2,'*-g')
    #pyplot.xlim(25,135)
    #pyplot.ylim(1600e9,3100e9)
    pyplot.legend([p1,p2],[ "watson", "brown"], loc=4)

    pyplot.show()
