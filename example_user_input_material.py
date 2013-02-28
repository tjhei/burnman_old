# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

"""
Shows user how to input a mineral of his/her choice and which physical values
need to be input for BurnMan to calculate Vs, Vp, Vphi and density at depth.

requires:
- creating minerals
- compute seismic velocities
- geotherms
- seismic models
- seismic comparison

teaches:
- how to create your own minerals

"""

import os, sys, numpy as np, matplotlib.pyplot as plt
#hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1,os.path.abspath('..')) 

import burnman
from burnman.minerals import material

if __name__ == "__main__":
    
    ### input variables ###
    #######################
    
    #INPUT for method
    method = 'slb' # choose 'slb' (finite-strain, stixrude and lithgow-bertelloni, 2005) or 'mgd' (mie-gruneisen-debeye, matas et al. 2007)
    
    class own_material (material):
            def __init__(self):
                    material.__init__(self)
                    self.params = {
                            'ref_V': 10.844e-6, #Unit cell volume [m^3/mole) at room pressure/temperature
                            'ref_K': 135.19e9, #Reference bulk modulus [Pa] at room pressure/temperature
                            'K_prime': 6.04, #pressure derivative of bulk modulus
                            'ref_mu': 175.0e9, #reference shear modulus at room pressure/temperature
                            'mu_prime': 1.7, #pressure derivative of shear modulus
                            'molar_mass': .055845, #molar mass in units of [kg/mol]
                            'n': 1, #number of atoms per molecule
                            'ref_Debye': 998.85, #Debye temperature for material. See Stixrude & Lithgow-Bertelloni, 2005 for values 
                            'ref_grueneisen': 1.368, #Gruneisen parameter for material. See Stixrude & Lithgow-Bertelloni, 2005 for values
                            'q0': 0.917, #q value used in caluclations. See Stixrude & Lithgow-Bertelloni, 2005 for values
                'eta_0s': 3.0} #eta value used in calculations. See Stixrude & Lithgow-Bertelloni, 2005 for values
    
    
    phases = [ own_material() ]
    molar_abundances = [ 1.0 ]
    
    
    #seismic model for comparison:
    seismic_model = burnman.seismic.prem() # pick from .prem() .slow() .fast() (see burnman/seismic.py)
    number_of_points = 20 #set on how many depth slices the computations should be done
    depths = np.linspace(700,2800, number_of_points)
    #depths = seismic_model.internal_depth_list()
    seis_p, seis_rho, seis_vp, seis_vs, seis_vphi = seismic_model.evaluate_all_at(depths)
    
            
    geotherm = burnman.geotherm.brown_shankland
    temperature = [geotherm(p) for p in seis_p]
    
    for ph in phases:
        ph.set_method(method)
    
    print "Calculations are done for:"
    for i in range(len(phases)):
        print molar_abundances[i], " of phase", phases[i].to_string()
    
    mat_rho, mat_vp, mat_vs, mat_vphi, mat_K, mat_mu = burnman.calculate_velocities(seis_p, temperature, phases, molar_abundances)    
    
    [rho_err,vphi_err,vs_err]=burnman.compare_with_seismic_model(mat_vs,mat_vphi,mat_rho,seis_vs,seis_vphi,seis_rho)
