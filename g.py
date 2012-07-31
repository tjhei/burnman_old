import numpy
import pylab
import scipy.linalg
import bisect
import math

# TODO: add up weight percent and check <100 and tell them how much


#print X
#    A = numpy.matrix([[ 1, 0, 1],
#    [ 0,1,0],
#    [ 0,1,1]])
#    b = numpy.matrix([[1],[2],[3]])
    
#    x = scipy.linalg.solve(A, b)
    #print x

def float_eq(a,b):
    return abs(a-b)<1e-10*max(1e-5,abs(a),abs(b))

def weight_pct_to_mol(element, amount):
    lower_mantle_mass = 4.043e27 # in g
    Av = 6.02214129e23 # in 1/mol

    molar_mass = {'Fe':55.845, 'Mg':24.305, 'O':15.999, 'Al':26.982, 'Ca':40.078, 'Si':28.085}
    return amount * lower_mantle_mass / molar_mass[element] * Av


def test_mol_conv():
    assert weight_pct_to_mol('Fe', 1.0) == 2*weight_pct_to_mol('Fe', 0.5)
    print float_eq(weight_pct_to_mol('Fe', 1.0), 4.35983834461e+49)


def conv_inputs(inp):
   

    names = {'Mg':'MgO','Fe':'FeO','Si':'SiO2', 'Ca':'Ca', 'Al':'Al'}
    out = {}
    for a in inp:
        out[names[a]] = weight_pct_to_mol(a,inp[a])
    return out
    


#
# inp = {'MgO':beta, 'FeO': , 'SiO2': gamma, 'Ca':, 'Al':} in mol
# params = {'Fe in pv': , 'Ca in pv':, 'Al in pv', 'Fe in fp':}
# returns: 'mol pv' A, 'mol fp' B, 'mol st' C in mol
# 'Mg in pv':0, 'Fe in pv':0, 'Ca in pv':0,'Si in pv':0, 'Al in pv':0
# 'Mg in fp':0,'Fe in fp':0
def determine_phases(inp, params):

    ret = {'mol pv':0, 'mol fp':0, 'mol st':0}
    ret['Mg in pv'] = 1-params['Fe in pv']-params['Ca in pv'] 
    ret['Fe in pv'] = params['Fe in pv']
    ret['Ca in pv'] = params['Ca in pv']
    ret['Si in pv'] = 1-params['Al in pv']
    ret['Al in pv'] = params['Al in pv']
    ret['Mg in fp'] = 1 - params['Fe in fp']
    ret['Fe in fp'] = params['Fe in fp']
 
    beta = inp['MgO']
    gamma = inp['SiO2']

    if (beta > gamma):
        ret['mol pv'] = beta - gamma
        ret['mol fp'] = beta - ret['mol pv']
        ret['mol st'] = 0
    elif (beta < gamma):
        ret['mol pv'] = beta
        ret['mol fp'] = 0
        ret['mol st'] = gamma - ret['mol pv']
    else:
        ret['mol pv'] = beta
        ret['mol fp'] = 0
        ret['mol st'] = 0
    
    return ret



inp = {'MgO':0.5, 'FeO': 0, 'SiO2':0.5, 'CaO':0, 'Al2O3':0.1}

params = {'Fe in pv': 0.8, 'Ca in pv':0.1, 'Al in pv':0, 'Fe in fp':1.0}

print determine_phases(inp, params)



def test_phases():
    # test everything into pv
    inp = {'MgO':20, 'FeO': 0, 'SiO2':20, 'CaO':0, 'Al2O3':0.0}
    params = {'Fe in pv': 0.0, 'Ca in pv':0.0, 'Al in pv':0.0, 'Fe in fp':0.0}
    t = determine_phases(inp, params)
    assert t['mol pv'] == 20
    assert t['mol fp'] == 0
    assert t['mol st'] == 0

    #
    inp = {'MgO':10, 'FeO': 0, 'SiO2':0, 'CaO':0, 'Al2O3':0.0}
    params = {'Fe in pv': 0.0, 'Ca in pv':0.0, 'Al in pv':0.0, 'Fe in fp':0.0}
    t = determine_phases(inp, params)
    assert t['mol pv'] == 10
    assert t['mol fp'] == 0
    assert t['mol st'] == 0

    #
    inp = {'MgO':10, 'FeO': 0, 'SiO2':3, 'CaO':0, 'Al2O3':0.0}
    params = {'Fe in pv': 0.0, 'Ca in pv':0.0, 'Al in pv':0.0, 'Fe in fp':0.0}
    t = determine_phases(inp, params)
    assert t['mol pv'] == 7
    assert t['mol fp'] == 3
    assert t['mol st'] == 0

    #
    inp = {'MgO':3, 'FeO': 0, 'SiO2':7, 'CaO':0, 'Al2O3':0.0}
    params = {'Fe in pv': 0.0, 'Ca in pv':0.0, 'Al in pv':0.0, 'Fe in fp':0.0}
    t = determine_phases(inp, params)
    assert t['mol pv'] == 3
    assert t['mol fp'] == 0
    assert t['mol st'] == 4



#input: pv, fp,
#return: bulk modulus, shear modulus, density
def eqn_of_state(inp):



    bla = 2.0

    out = {}
    out['density']= lambda pressure: 1+bla*pressure

    return out



geotherm_table=[]
for line in open("geotherm.txt").readlines():
    if (line[0]!='#'):
	numbers = map(float, line.split())
	geotherm_table.append(numbers)
geotherm_p=numpy.array(geotherm_table)[:,0]
geotherm_T=numpy.array(geotherm_table)[:,1]



# pressure: in GPa
# return: temperature
def geotherm(pressure):
    idx = bisect.bisect_left(geotherm_p, pressure)
    #    print geotherm_p[idx]

    return geotherm_T[idx]



# test geotherm
#p= numpy.arange(1.0,128.0,3)
#t= [geotherm(y) for y in p]
#pylab.plot(p,t,'+-')
#pylab.show()


# rho: density
# ref_rho: reference density
# ref_K: 
# K_prime: deriv. bul modulus
# return: pressure in GPa
def birch_murnaghan(rho, ref_rho, ref_K, K_prime):
    x = rho/ref_rho
    return 3.*ref_K/2. * (pow(x, 7./3.) - pow(x, 5./3.)) \
        * (1 + 0.75*(K_prime-4.)*(pow(x, 2./3.) - 1.))






# Reuss-Voigt-Hill average
# inp: bulk modulus, shear modulus, density
# return K_Si and gamma_i 
def voigt_reuss_hill(molar_abundance, molar_weight, modulus, thermal_exp, grueneisen, density, T):
    # from matas 2007, Appendix C

    n_phase = len(modulus)

    assert n_phase == len(molar_abundance)
    assert n_phase == len(molar_weight)
    assert n_phase == len(density)

    it = range(n_phase)

    n_i = molar_abundance  # molar abundance for phase i 
    M_i = molar_weight  # molar weight for phase i
    total_molar = sum(n_i)
    x_i = [(n_i[i] / total_molar) for i in it] # molar fraction of phase i
    
    #V_i = n_i * M_i / rho:
    V_i = [(n_i[i]*M_i[i]/density[i]) for i in it]
    #V = sum n_i V_i:
    V = sum(i*j for i, j in zip(n_i,V_i))
    #nu_i = x_i * V_i / V:
    nu_i = [(x_i[i] * V_i[i] / V) for i in it]

    #X_i = K_Si = K_Ti (1+alpha_i gamma_i T)
    #alpha_i: thermal expansion
    #gamma_i: grueneisen
    #K_Ti: bulk modulus    
    #not needed: X_i = [ (modulus[i] * (1.+thermal_exp[i] * grueneisen[i] * T)) for i in it]
    X_i = modulus

    #X_V = sum nu_i X_i:
    X_V = sum(i*j for i,j in zip(nu_i,X_i))
    #X_R = 1 / sum(nu_i/X_i):
    X_R = 1. / sum(i/j for i,j in zip(nu_i,X_i))

    return (X_V + X_R) / 2


def calc_velocities(molar_abundance, molar_weight, bulk_modulus, shear_modulus, thermal_exp, grueneisen,  density, T):

    it = range(len(molar_abundance))
    n_i = molar_abundance  # molar abundance for phase i 
    M_i = molar_weight  # molar weight for phase i
    V_i = [(n_i[i]*M_i[i]/density[i]) for i in it]
    #V = sum n_i V_i:
    V = sum(i*j for i, j in zip(n_i,V_i))    
    # avg_density = 1./ V sum(n_i M_i)
    avg_density = 1./ V * sum((n_i[i]*M_i[i]) for i in it)

    K_s = voigt_reuss_hill(molar_abundance, molar_weight, bulk_modulus, thermal_exp, grueneisen, density, T)
    mu = voigt_reuss_hill(molar_abundance, molar_weight, shear_modulus, thermal_exp, grueneisen, density, T)
    V_p = math.sqrt(K_s + 4./3. * mu / avg_density )
    V_s = math.sqrt(mu / avg_density)
    V_phi = math.sqrt(K_s / avg_density)

    return V_p,V_s,V_phi


molar_abundance=[1., 1., 1.]
molar_weight=[1., 1., 1.]
bulk_mod =[1., 1., 1.]
shear_mod =[1., 1., 1.]
thermal_exp = [1., 1., 1.]
grueneisen = [1., 1., 1.]
density = [1., 1., 1.]
T=1

print calc_velocities(molar_abundance, molar_weight, bulk_mod, shear_mod, thermal_exp, grueneisen,  density, T)



test_phases()
test_mol_conv()












inp1 = {'Mg':0.5, 'Fe': 0, 'Si':0.5, 'Ca':0.0, 'Al':0} # wt%
inp2 = conv_inputs(inp1)
print "in:", inp1
print "out:", inp2

params = {'Fe in pv': 0.0, 'Ca in pv':0.0, 'Al in pv':0.0, 'Fe in fp':0.0}
t = determine_phases(inp2, params)
print "phases:", t

ret = eqn_of_state(t)

print "eos:", ret

print ret['density'](42)








