# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

#system libs:
import numpy
import scipy.optimize as opt
import scipy.integrate as integrate
import math
import matplotlib.pyplot as plt

#own libs:
from tools import *
import mie_grueneisen_debye as mgd
import birch_murnaghan as bm
# TODO: add up weight percent and check <100 and tell them how much
## Based on Stixrude & Lithgow-Bertelloni (2005), all equation numbers refer to this paper. 

gas_constant=8.314462175

def volume(p,T,params):
    P_th = lambda x: mgd.mgd_thermal_pressure(T,x, params) #From Stixrude 2005: delta U (T)*gamma*ro (1/V)
    P_th_ref = lambda x: mgd.mgd_thermal_pressure(300.,x, params)#From Stixrude 2005: delta U (300 K)*gamma*ro (1/V)
      
    b_iikk= 9.*params['ref_K'] # EQ 28
    b_iikkmm= 27.*params['ref_K']*(params['K_prime']-4.) # EQ 29
    f = lambda x: 0.5*(pow(params['ref_V']/x,2./3.)-1.) # EQ 24

    func = lambda x: (1./3.)*(pow(1.+2.*f(x),5./2.))*((b_iikk*f(x)) \
            +(0.5*b_iikkmm*pow(f(x),2.))) + P_th(x) - P_th_ref(x) - p #EQ 21 
    
    V = opt.brentq(func, 0.09*params['ref_V'], 3.*params['ref_V']) 
    return V

def shear_modulus(T,V,params):

        # low T:
        # C_v = 234. * n * Av * boltzmann_constant (T/theta_D) ^ 3.  (3.43 from Poirier/EarthInterior)
        # high T:
        # C_v = 3. * n * gas_constant
        # n = number of moles
    #C_v = 3. * n * gas_constant # in J * mol / K
    P_th = lambda x: mgd.mgd_thermal_pressure(T,x, params)
    P_th_ref = lambda x: mgd.mgd_thermal_pressure(300.,x, params)
    f =.5*(pow(params['ref_V']/V,2./3.)-1.) # EQ 24
    
    a2_s = -2.*params['ref_grueneisen'] - 2.*params['eta_0s'] # EQ 47 
    a1_ii = 6. * params['ref_grueneisen'] # EQ 47
    a2_iikk = -12.*params['ref_grueneisen']+36.*pow(params['ref_grueneisen'],2.) - 18.*params['q0']*params['ref_grueneisen'] # EQ 47
    nu_o_nu0_sq = 1.+ a1_ii*f + (1./2.)*a2_iikk * pow(f,2.) # EQ 41
    #gamma = params['ref_grueneisen'] * (pow(params['ref_V']/V,-params['q0'])) # EQ 49
    gamma = 1./6.*pow(nu_o_nu0_sq,-1.)*(2*f+1.)*(a1_ii+(a2_iikk*f))
    eta_s = - gamma - (1./2. * pow(nu_o_nu0_sq,-1.) * pow((2.*f)+1.,2.)*a2_s) # EQ 46 NOTE the typo from Stixrude 2005
    
    delta_U = (P_th(V) - P_th_ref(V))*(V/gamma)
    
  
    G_threeterms = pow(1.+2.*f, 5./2.) * (params['ref_mu'] + (3.*(params['ref_K']*params['mu_prime']) - 5.*params['ref_mu'])*f \
                                   + ((6.*params['ref_K']*params['mu_prime']) - 24.*params['ref_K'] - 14.*params['ref_mu']+ 9./2.*params['ref_K']*params['K_prime']) * pow(f,2)) \
                                   - (eta_s*(1./V)*delta_U) # EQ 33 up to the second order

    return G_threeterms
    

def bulk_modulus(T,V,params):
    func_int = lambda t: pow(t,3.)/(math.exp(t)-1.)
    
    P_th = lambda x: mgd.mgd_thermal_pressure(T,x, params)
    P_th_ref = lambda x: mgd.mgd_thermal_pressure(300.,x, params)
    
    func_int_cv = lambda t: math.exp(t)*pow(t,4.)/pow(math.exp(t)-1.,2.)
    Dcv = integrate.quad(func_int_cv,0.,params['ref_Debye']/T) 
    Dcv_300 = integrate.quad(func_int_cv,0.,params['ref_Debye']/300.)
    cv_300 = 9.*params['n']*gas_constant*(pow(300./params['ref_Debye'],3.))*Dcv_300[0]
    Cv = (9.*params['n']*gas_constant*(pow(T/params['ref_Debye'],3.))*Dcv[0])-cv_300 #why substract this? Sanne #units of R

    #Kth_0 = params['ref_K'] - (pow(params['ref_grueneisen'],2.) * density_0 * delta_U_300/params['molar_mass']) 
    f =.5*(pow(params['ref_V']/V,2./3.)-1) # EQ 24
    gamma = params['ref_grueneisen'] * (pow(params['ref_V']/V,-params['q0'])) # EQ 49
    
    delta_U = (P_th(V) - P_th_ref(V))*(V/gamma) #convert into Stixrudian units
    K = pow(1.+2.*f, 5./2.) * ( params['ref_K'] + (3*params['ref_K']*params['K_prime'] -5.*params['ref_K'])*f) \
             + ((gamma+1.-params['q0'])*(gamma/V)  *delta_U ) \
            - ((pow(gamma,2.) / V )*(Cv*T - cv_300*300.) )   # EQ 32 up to the second order (the third order is basically zero when K'~4
    
    #+27./2.*(params['ref_K']*params['K_prime']-4.*params['ref_K'])*pow(f,2.)) \
    return K

#heat capacity at constant volume
def heat_capacity_v(T, V, params):
    return mgd.heat_capacity_v(T,V,params)

def thermal_expansivity(T,V,params):
    C_v = heat_capacity_v(T,V,params)
    gr = mgd.grueneisen_parameter(params['ref_V']/V, params)
    K = bulk_modulus(T,V,params)
    alpha = gr * C_v / K / V
    return alpha

#heat capacity at constant pressure
def heat_capacity_p(T,V,params):
    alpha = thermal_expansivity(T,V,params)
    gr = mgd.grueneisen_parameter(params['ref_V']/V, params)
    C_v = heat_capacity_v(T,V,params)
    C_p = C_v*(1. + gr * alpha * T)
    return C_p

#calculate the adiabatic bulk modulus (K_S) as a function of P, T, and V
# alpha is basically 1e-5
def bulk_modulus_adiabatic(T,V,params):
    K_T=bulk_modulus(T,V,params)
    alpha = thermal_expansivity(T,V,params)
    gr = mgd.grueneisen_parameter(params['ref_V']/V, params)
    K_S = K_T*(1. + gr * alpha * T)
    return K_S


