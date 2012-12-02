from material import *
from composition import calculate_partition_coefficient

################ User input minerals

class user_mineral1 (material):
        def __init__(self):
                material.__init__(self)
                self.params = {
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
class user_mineral2 (material):
        def __init__(self):
                self.params = {
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
class user_mineral3 (material):
        def __init__(self):
                self.params = {
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

################ Minerals
class test_mineral (material):
	def __init__(self):
		self.params = {
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
		self.params = {
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
		self.params = {
			'ref_V': 11.24e-6,
			'ref_K': 161.,
			'K_prime': 3.9,
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
		self.params = {
			'ref_V': 12.06e-6,
			'ref_K': 152.,
			'K_prime': 4.9,
			'ref_mu': 47.,
			'mu_prime': 0.7,
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
		self.params = {
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
		self.params = {
			'ref_V': 24.45e-6,
			'ref_K': 251.,
			'K_prime': 4.1,
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
		self.params = {
			'ref_V': 25.48e-6,
			'ref_K': 281.,
			'K_prime': 4.1,
			'ref_mu': 161.,
			'mu_prime': 1.57,
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
		self.params = {
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
		self.params = {
			'ref_V': 24.607e-6,
			'ref_K': 251.9,
			'K_prime': 4.01,
			'ref_mu': 164.7,
			#'ref_mu': 157.39,  #refitted to second order  
    		        'mu_prime': 1.58,
			#'mu_prime': 2.08, #refitted to second order
			'molar_mass': .102165,
			'n': 5,
			'ref_Debye': 1054.,
			'ref_grueneisen': 1.48,
			'q0': 1.4, 
			'eta_0s': 2.4 }
			
class Murakami_fp_HS(material): #From Murakami's emails, see Cayman for details, represents Mg# = .79
	def __init__(self):
		self.params = {
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
		self.params = {
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

class fe_dependent_helper(material):
	def __init__(self, iron_number_with_pt, idx):
		self.iron_number_with_pt = iron_number_with_pt
		self.which_index = idx # take input 0 or 1 from iron_number_with_pt()

	def create_inner_material(self, iron_number):
		return [] # needs to be overwritten in class deriving from this one

        def set_state(self, pressure, temperature):
		self.pressure = pressure
		self.temperature = temperature
		self.fp = self.create_inner_material(self.iron_number())
		self.fp.method = self.method
		self.fp.set_state(pressure, temperature)
		self.params = self.fp.params
		material.set_state(self, pressure, temperature)

	def iron_number(self):
		return self.iron_number_with_pt(self.pressure,self.temperature)[self.which_index]
	def molar_mass(self):
		return self.fp.molar_mass()
	def density(self):
		return self.fp.density()
	def molar_volume(self):
		return self.fp.molar_volume()
	def bulk_modulus(self):
		return self.fp.bulk_modulus()
	def v_s(self):
		return self.fp.v_s()
	def v_p(self):
		return self.fp.v_p()
	def geotherm(self):
		return self.fp.v_s()


class mg_fe_perovskite_pt_dependent(fe_dependent_helper):
	def __init__(self, iron_number_with_pt, idx):
		fe_dependent_helper.__init__(self, iron_number_with_pt, idx)

	def create_inner_material(self, iron_number):
		return mg_fe_perovskite(iron_number)

class ferropericlase_pt_dependent(fe_dependent_helper):
	def __init__(self, iron_number_with_pt, idx):
		fe_dependent_helper.__init__(self, iron_number_with_pt, idx)

	def create_inner_material(self, iron_number):
		return ferropericlase(iron_number)
