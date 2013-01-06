"""
    BurnMan- a lower mantle toolkit
    Copyright (C) 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
"""


from minerals import *
from composition import calculate_partition_coefficient


class test_mineral (material):
	def __init__(self):
		self.params = {	'name':'test',
			'ref_V': 6.844e-6,
			'ref_K': 135.19,
			'K_prime': 6.04,
			'ref_mu': 175.,
			'mu_prime': 1.7,
			'molar_mass': .055845,
			'n': 1,
			'ref_Debye': 998.85,
			'ref_grueneisen': 1.368,
			'q0': 0.917}
 
class stishovite (material):
	def __init__(self):
		self.params = {	'name':'stishovite',
			'ref_V': 14.02e-6,
			'ref_K': 314.,
			'K_prime': 4.4,
			'ref_mu': 220.,
			'mu_prime': 1.6,
			'molar_mass': .0601,
			'n': 3,
			'ref_Debye': 1044.,
			'ref_grueneisen': 1.6,
			'q0': 2.4,
			'eta_0s': 3.0 }

class periclase (material):
	def __init__(self):
		self.params = {	'name':'periclase',
			'ref_V': 11.25e-6,
			'ref_K': 160.,
			'K_prime': 3.83,
			'ref_mu': 130.,
			'mu_prime': 2.2,
			'molar_mass': .0403,
			'n': 2,
			'ref_Debye': 773.,
			'ref_grueneisen': 1.5,
			'q0': 1.5,
			'eta_0s': 3.0 }
class wustite (material):
	def __init__(self):
		self.params = {	'name':'periclase',
			'ref_V': 12.26e-6,
			'ref_K': 160.,
			'K_prime': 3.83,
			'ref_mu': 46.,
			'mu_prime': 0.6,
			'molar_mass': .0718,
			'n': 2,
			'ref_Debye': 455.,
			'ref_grueneisen': 1.28,
			'q0': 1.5,
			'eta_0s': 3.0 }



class ferropericlase(material):
	def __init__(self, mg_num):
		self.mg = mg_num
		self.fe = 1.0-mg_num
		self.pe = periclase()
		self.wu = wustite()
		self.params = {'name':'ferropericlase',
			'ref_V': self.pe.params['ref_V']*self.mg + self.wu.params['ref_V']*self.fe,
			'ref_K': self.pe.params['ref_K']*self.mg + self.wu.params['ref_K']*self.fe,
			'K_prime': self.pe.params['K_prime']*self.mg + self.wu.params['K_prime']*self.fe,
			'ref_mu': self.pe.params['ref_mu']*self.mg + self.wu.params['ref_mu']*self.fe,
			'mu_prime': self.pe.params['mu_prime']*self.mg + self.wu.params['mu_prime']*self.fe,
			'molar_mass': self.pe.params['molar_mass']*self.mg + self.wu.params['molar_mass']*self.fe,
			'n': 5,
			'ref_Debye': self.pe.params['ref_Debye']*self.mg + self.wu.params['ref_Debye']*self.fe,
			'ref_grueneisen': self.pe.params['ref_grueneisen']*self.mg + self.wu.params['ref_grueneisen']*self.fe,
			'q0': self.pe.params['q0']*self.mg + self.wu.params['q0']*self.fe ,
			'eta_0s': self.pe.params['eta_0s']*self.mg + self.wu.params['eta_0s']*self.fe }


class mg_perovskite(material):
	def __init__(self):
		self.params = {	'name':'Mg perovskite',
			'ref_V': 24.43e-6,
			'ref_K': 250.,
			'K_prime': 4.0,
			#'ref_mu': 166.,
			#'mu_prime': 1.57,
                        'ref_mu': 175., #from S & L.-B. 2005
                        'mu_prime': 1.8,
			'molar_mass': .1020,
			'n': 5,
			'ref_Debye': 1070.,
			'ref_grueneisen': 1.48,
			'q0': 1.4, 
			'eta_0s': 2.4 }

class fe_perovskite(material):
	def __init__(self):
		self.params = {	'name':'Fe perovskite',
			'ref_V': 25.48e-6,
			'ref_K': 250.,
			'K_prime': 4.0,
			'ref_mu': 135.,
			'mu_prime': 1.3,
			'molar_mass': .1319,
			'n': 5,
			'ref_Debye': 1021.,
			'ref_grueneisen': 1.48,
			'q0': 1.4, 
			'eta_0s': 2.4 }

class mg_fe_perovskite(material):
	def __init__(self, mg_num):
		self.mg = mg_num
		self.fe = 1.0-mg_num
		self.mg_pv = mg_perovskite()
		self.fe_pv = fe_perovskite()
		self.params = {'name':'Mg-Fe perovskite',
			'ref_V': self.mg_pv.params['ref_V']*self.mg + self.fe_pv.params['ref_V']*self.fe,
			'ref_K': self.mg_pv.params['ref_K']*self.mg + self.fe_pv.params['ref_K']*self.fe,
			'K_prime': self.mg_pv.params['K_prime']*self.mg + self.fe_pv.params['K_prime']*self.fe,
			'ref_mu': self.mg_pv.params['ref_mu']*self.mg + self.fe_pv.params['ref_mu']*self.fe,
			'mu_prime': self.mg_pv.params['mu_prime']*self.mg + self.fe_pv.params['mu_prime']*self.fe,
			'molar_mass': self.mg_pv.params['molar_mass']*self.mg + self.fe_pv.params['molar_mass']*self.fe,
			'n': 5,
			'ref_Debye': self.mg_pv.params['ref_Debye']*self.mg + self.fe_pv.params['ref_Debye']*self.fe,
			'ref_grueneisen': self.mg_pv.params['ref_grueneisen']*self.mg + self.fe_pv.params['ref_grueneisen']*self.fe,
			'q0': self.mg_pv.params['q0']*self.mg + self.fe_pv.params['q0']*self.fe ,
			'eta_0s': self.mg_pv.params['eta_0s']*self.mg + self.fe_pv.params['eta_0s']*self.fe}

class Murakami_perovskite(material): #From Murakami's emails, see Cayman for details, represents 4 wt% Al X_mg = .94
	def __init__(self):
		self.params = {	'name':'Murakami Mg-Fe perovskite',
			'ref_V': 24.607e-6,
			'ref_K': 251.9,
			'K_prime': 4.01,
			#'ref_mu': 164.7,
			'ref_mu': 157.39,  #refitted to second order  
    		#'mu_prime': 1.58,
			'mu_prime': 2.08, #refitted to second order
			'molar_mass': .102165,
			'n': 5,
			'ref_Debye': 1054.,
			'ref_grueneisen': 1.48,
			'q0': 1.4, 
			'eta_0s': 2.4 }
			
class Murakami_fp_HS(material): #From Murakami's emails, see Cayman for details, represents Mg# = .79
	def __init__(self):
		self.params = {	'name':'Murakami Mg-Fe mw HS',
			'ref_V': 11.412e-6,
			'ref_K': 159.,
			'K_prime': 4.11,
			'ref_mu': 105.43,
			'mu_prime': 1.773,
			'molar_mass': .0494,
			'n': 2,
			'ref_Debye': 706.,
			'ref_grueneisen': 1.5,
			'q0': 1.5, 
			'eta_0s': 3.0 }
			
class Murakami_fp_LS(material): #From Murakami's emails, see Cayman for details, represents Mg# = .79
	def __init__(self):
		self.params = {	'name':'Murakami Mg-Fe mw LS',
			'ref_V': 11.171e-6,
			'ref_K': 170.,
			'K_prime': 4.00,
			'ref_mu': 116.34,
			'mu_prime': 1.668,
			'molar_mass': .0494,
			'n': 2,
			'ref_Debye': 706., 
			'ref_grueneisen': 1.5,
			'q0': 1.5, 
			'eta_0s': 3.0}

class mg_fe_perovskite_pt_dependent(material):
	def __init__(self, iron_number_with_pt):
		self.iron_number_with_pt = iron_number_with_pt
		self.params = {'name':'Depth dependent Mg Fe perovskite',
				'n':5 }
	def iron_number(self,pressure,temperature):
		return self.iron_number_with_pt(pressure,temperature)[1]
	def molar_mass(self,pressure,temperature):
		pv = mg_fe_perovskite(1.-self.iron_number_with_pt(pressure,temperature)[1])
		self.params = pv.params
		return pv.molar_mass(pressure,temperature)
	def density(self, pressure, temperature):
		pv = mg_fe_perovskite(1.-self.iron_number_with_pt(pressure,temperature)[1])
		self.params = pv.params
		return pv.density(pressure, temperature)
	def molar_volume(self,pressure, temperature):
		pv = mg_fe_perovskite(1.-self.iron_number_with_pt(pressure,temperature)[1])
		self.params = pv.params
		return pv.molar_volume(pressure, temperature)
	def bulk_modulus(self,pressure, temperature):
		pv = mg_fe_perovskite(1.-self.iron_number_with_pt(pressure,temperature)[1])
		self.params = pv.params
		return pv.bulk_modulus(pressure, temperature)
	def v_s(self,pressure, temperature):
		pv = mg_fe_perovskite(1.-self.iron_number_with_pt(pressure,temperature)[1])
		self.params = pv.params
		return pv.v_s(pressure, temperature)
	def v_p(self,pressure, temperature):
		pv = mg_fe_perovskite(1.-self.iron_number_with_pt(pressure,temperature)[1])
		self.params = pv.params
		return pv.v_p(pressure, temperature)
	def geotherm(self,pressure, temperature):
		pv = mg_fe_perovskite(1.-self.iron_number_with_pt(pressure,temperature)[1])
		self.params = pv.params
		return pv.v_s(pressure, temperature)

class ferropericlase_pt_dependent(material):
	def __init__(self, iron_number_with_pt):
		self.iron_number_with_pt = iron_number_with_pt
	def iron_number(self,pressure,temperature):
		return self.iron_number_with_pt(pressure,temperature)[0]
	def molar_mass(self,pressure,temperature):
		fp = ferropericlase(1.-self.iron_number_with_pt(pressure,temperature)[0])
		self.params = fp.params
		return fp.molar_mass(pressure,temperature)
	def density(self,pressure, temperature):
		fp = ferropericlase(1.-self.iron_number_with_pt(pressure,temperature)[0])
		self.params = fp.params
		return fp.density(pressure, temperature)
	def molar_volume(self,pressure, temperature):
		fp = ferropericlase(1.-self.iron_number_with_pt(pressure,temperature)[0])
		self.params = fp.params
		return fp.molar_volume(pressure, temperature)
	def bulk_modulus(self,pressure, temperature):
		fp = ferropericlase(1.-self.iron_number_with_pt(pressure,temperature)[0])
		self.params = fp.params
		return fp.bulk_modulus(pressure, temperature)
	def v_s(self,pressure, temperature):
		fp = ferropericlase(1.-self.iron_number_with_pt(pressure,temperature)[0])
		self.params = fp.params
		return fp.v_s(pressure, temperature)
	def v_p(self,pressure, temperature):
		fp = ferropericlase(1.-self.iron_number_with_pt(pressure,temperature)[0])
		self.params = fp.params
		return fp.v_p(pressure, temperature)
	def geotherm(self,pressure, temperature):
		fp = ferropericlase(1.-self.iron_number_with_pt(pressure,temperature)[0])
		self.params = fp.params
		return fp.v_s(pressure, temperature)

