import numpy as np
import seismic

# Voigt-Reuss-Hill average
# inputs are pressure, temperature, a list of phases, and
# a list of the same length of molar abundances of the phases.
# That is to say, given a mole of molecules, what fraction of them
# are in which phase.  They should add up to 1.

def voigt_reuss_hill(pressure, temperature, phases, molar_abundances):
    # from matas 2007, Appendix C

    n_phase = len(phases)
    assert n_phase == len(molar_abundances)
    assert sum(molar_abundances)  > 0.999
    assert sum(molar_abundances)  < 1.001
    it = range(n_phase)
    
    #determine the partial molar volumes of the phases
    for i in it:
      phases[i].set_state(pressure,temperature)

    V_i = [(molar_abundances[i]*phases[i].molar_volume()) for i in it]
    V_tot = sum(V_i)
    
    #calculate the density of the phase assemblage
    rho = (1./V_tot) * sum( molar_abundances[i]*phases[i].molar_mass() for i in it)

    #calculate the voigt and reuss averages of the phase assemblage for K and mu
    K_voigt = sum( V_i[i]/V_tot * phases[i].adiabatic_bulk_modulus() for i in it)
    K_reuss = 1./sum( V_i[i]/V_tot / phases[i].adiabatic_bulk_modulus() for i in it)
   
    mu_voigt = sum( V_i[i]/V_tot * phases[i].shear_modulus() for i in it)
    mu_reuss = 1./sum( V_i[i]/V_tot / phases[i].shear_modulus() for i in it)

    #average voigt and reuss for vrh
    K_vrh = (K_voigt+K_reuss)/2.
    mu_vrh = (mu_voigt+mu_reuss)/2.

    #compute seismic velocities
    v_s = np.sqrt( mu_vrh*1.e9 / rho)
    v_p = np.sqrt( (K_vrh*1.e9 + 4./3.*mu_vrh*1e9) / rho)
    v_phi = np.sqrt( (K_vrh*1e9) / rho)

    return rho, v_p, v_s, v_phi, K_vrh, mu_vrh

#calculate the voigt and reuss averages of the phase assemblage for a certain property X
#from: matas 2007, appendix D
def vhr_average(phase_volume, X):

    it = range(len(phase_volume))
    
    V_i = phase_volume
    V_tot = sum(V_i)

    X_voigt = sum( V_i[i]/V_tot * X[i] for i in it)
    X_reuss = 1./sum( V_i[i]/V_tot / X[i] for i in it)
    X_vrh = (X_voigt+X_reuss)/2.
    return X_vrh

def attenuation_correction(seismic_model, pressure,v_p,v_s,v_phi):
    Qphi, Qs = seismic_model.attenuation(pressure)
  
    beta = 0.3 # Matas et al. (2007) page 4        
    Qp  = 3./4.*pow((v_p/v_s),2.)*Qs    # Matas et al. (2007) page 4

    cot=1./np.tan(beta*np.pi/2.)
    v_p  = v_p*(1.-1./2.*cot*1./Qp)    # Matas et al. (2007) page 1
    v_s  = v_s*(1.-1./2.*cot*1./Qs)
    v_phi= v_phi*(1.-1./2.*cot*1./Qphi)

    # Kanamori and Anserson 1971, correcting from GHz to Hz
    #v_p  = v_p*(1.+1./(Qp*np.pi)*np.log(1./1.e9))    
    #v_s  = v_s*(1.+1./(Qs*np.pi)*np.log(1./1.e9))
    #v_phi= v_phi*(1.+1./(Qphi*np.pi)*np.log(1./1.e9))

    return v_p, v_s, v_phi

