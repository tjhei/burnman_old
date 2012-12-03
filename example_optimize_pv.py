import os, sys, numpy as np
import matplotlib.pyplot as plt

lib_path = os.path.abspath('code/')
sys.path.append(lib_path)

from code import minerals 
from code import main as main

### input variables ###
#######################

#INPUT for method
method = 'slb' # choose 'slb' (finite-strain, stixrude and lithgow-bertelloni, 2005) or 'mgd' (mie-gruneisen-debeye, matas et al. 2007)

#INPUT for seismic models
name = 'prem'  					# choose from 'prem'(at 1 second, Dziewonski & Anderson, 1981), 'ref_fast' (Lekic et al. 2012), 'ref_slow' (Lekic et al. 2012)
attenuation_correction= 'off'   		# correction for attuation (not required for PREM at 1s, see Matas et al. 2007, page 4) choose from 'on', 'off'
depth_min = 0.0   				# minimum depth considered in km, choose 0.0 to use the limits given by the seismic model
depth_max = 0.0   				# maximum depth considered in km, choose 0.0 to use the limits given by the seismic model
depth_step = 100.0 				# steps at which the seismic velocity and calculated velocity are compared (in between is interpolated)

seis_p, seis_r, seis_vp, seis_vphi, seis_vs, seis_rho = main.seismo_load(name, attenuation_correction, depth_min, depth_max,depth_step)

geotherm = main.get_geotherm("brown_shankland")
temperature = [geotherm(p) for p in seis_p]

def material_error(amount_perovskite):
	phases = (minerals.Murakami_perovskite(), minerals.Murakami_fp_LS())
	molar_abundances = ( amount_perovskite, 1.0-amount_perovskite)

	for ph in phases:
		ph.set_method(method)

	print "Calculations are done for:"
	for i in range(len(phases)):
		print molar_abundances[i], " of phase", phases[i].to_string()

	mat_rho, mat_vs, mat_vp, mat_vphi, mat_K, mat_mu = main.calculate_velocities(seis_p, temperature, phases, molar_abundances)	

	[rho_err,vphi_err,vs_err]=main.compare_with_seismic_model(mat_vs,mat_vphi,mat_rho,seis_vs,seis_vphi,seis_rho/1e3) #check units on rho
	return vs_err, vphi_err

xx=np.linspace(0.0, 1.0, 20)
yy_vs=[material_error(x)[0] for x in xx]
yy_vphi=[material_error(x)[1] for x in xx]
plt.plot (xx,yy_vs,label=("vs error"))
plt.plot (xx,yy_vphi,label=("vphi error"))
plt.ylim(0,100)
plt.xlabel('% Perovskite')
plt.ylabel('Error')
plt.legend()
plt.show()
