#system libs:
import numpy
import pylab
import scipy.optimize as opt
import math

#own libs:
import geotherm
import prem
from tools import *
from eos_from_ian import bm_density, birch_murnaghan
import seismic

# TODO: add up weight percent and check <100 and tell them how much

molar_mass = {'Fe':55.845, 'Mg':24.305, 'O':15.999, 'Al':26.982, 'Ca':40.078, 'Si':28.085} # g/mol
Av = 6.02214129e23 # Avogadro constant in 1/mol 





# convert weight percentage (amount, 1.00 = 100%) of a given element to molar mass
def weight_pct_to_mol(element, amount):
    lower_mantle_mass = 4.043e27*.75 # in g

    return amount * lower_mantle_mass / molar_mass[element] * Av


def test_mol_conv():
    assert weight_pct_to_mol('Fe', 1.0) == 2*weight_pct_to_mol('Fe', 0.5)
    assert float_eq(weight_pct_to_mol('Fe', 1.0), 3.26987875846e+49)


def conv_inputs(inp):
    names = {'Mg':'MgO','Fe':'FeO','Si':'SiO2', 'Ca':'Ca', 'Al':'Al'}
    out = {}
    for a in inp:
        out[names[a]] = weight_pct_to_mol(a,inp[a])
    return out
    

# compute phases of pv, fp, st
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



#input: pv, fp, st in mol
#return: bulk modulus, shear modulus, density
def eqn_of_state(inp):
    # placeholder for now
    bla = 2.0

    out = {}
    out['density']= lambda pressure: 1+bla*pressure

    return out








#murakami test:

def test(p,T,V0,K_0,K_prime,dKdT,a_0,a_1,gamma_0):
    K_T = K_0 + dKdT*(T-300)
    alpha = a_0 + a_1*T
    P_th = alpha * K_T*(T-300)
    func = lambda x: birch_murnaghan (V0/x, 1., K_0, K_prime) + P_th - p
    #xx = numpy.arange(0.1, V0, V0/100.)
    #yy = [func(x) for x in xx]
    #pylab.plot(xx,yy,'o-r')
    #pylab.show()

    V = opt.brentq(func, 0.1, V0)
    bulk_mod = K_0 + K_prime * p + dKdT*(T-300.)
    bulk_mod = bulk_mod * (1. + alpha * gamma_0 * T)        # formula D6
    return V, bulk_mod


def murakami(molar_abundance):
    al_p = 0.075
    pv_X_Mg = 0.94
    fp_X_Mg = 0.79
    molar_weight=[pv_X_Mg*molar_mass['Mg']+(1.-pv_X_Mg)*molar_mass['Fe']+(1.-al_p)*molar_mass['Si']+al_p*molar_mass['Al']+3.*molar_mass['O'], \
                  fp_X_Mg*molar_mass['Mg']+(1.-fp_X_Mg)*molar_mass['Fe']+molar_mass['O']]

    list_p = []
    list_Vs = []
    list_Vp = []
    pv_density = []
    fp_density = []
    pv_shearmod = []
    fp_shearmod = []
    prem_shearmod = []

    for p in range(30,141,1):
        #T=geotherm.geotherm_formula(p)
        #T=geotherm.geotherm(p)
        T=geotherm.geotherm_brown(p) #by far the best fit

        density = [0., 0.]
        bulk_mod =[0., 0.]
        shear_mod = [0., 0.]

        # pv:
        # values from ricolleau table 1
        V0 = 164.
        K_0 = 245.
        K_prime = 4.
        dKdT = -.036
        a_0 = 3.19e-5
        a_1 = 0.88e-8
        gamma_0 = 1.48

        V, bulk_mod[0] = test(p,T,V0,K_0,K_prime,dKdT,a_0,a_1,gamma_0)
        atoms_per_unit_cell = 4. # correct for pv
        density[0] = molar_weight[0]*atoms_per_unit_cell / (Av*V*1e-24)    #correct according to jc and cayman

        pv_density.append(density[0])

        # fp:
        # values from ricolleau table 1
        if (p<=50.):
            V0 = 76.44
            K_0 = 158.
            K_prime = 4.
            dKdT = -.034
            a_0 = 2.20e-5
            a_1 = 3.61e-8
        else:
            V0 = 74.04
            K_0 = 170.
            K_prime = 4.
            dKdT = -.034
            a_0 = 2.20e-5
            a_1 = 3.61e-8
            #missing values, should we use 0?
            #dKdT = 0.
            #a_0 = 0.
            #a_1 = 0.

        gamma_0 = 1.50

        V, bulk_mod[1] = test(p,T,V0,K_0,K_prime,dKdT,a_0,a_1,gamma_0)
        atoms_per_unit_cell = 4. # correct for fp
        density[1] = molar_weight[1]*atoms_per_unit_cell/(Av*V*1e-24)    #correct according to jc and cayman
        fp_density.append(density[1])
        
        shear_mod[0] = 166. + 1.57*p -.02*(T-300)

        #if (p>=50.): # from the text
        #    shear_mod[1] = 130. + 2.04*p -.02*(T-300) # low spin
        #else:
        #    shear_mod[1] = 113. + 2.15*p -.02*(T-300)

        if (p>=50.): # reading from fig 3 in Murakami
            shear_mod[1] = 116. + 1.65*p -.02*(T-300.)
        else:
            shear_mod[1] = 103. + 1.78*p -.02*(T-300.)


        pv_shearmod.append(shear_mod[0])
        fp_shearmod.append(shear_mod[1])

        result = seismic.get_velocities(molar_abundance, molar_weight, bulk_mod, shear_mod, density, T)
        list_p.append(p)
        list_Vp.append(result[0])
        list_Vs.append(result[1])

        # shearmodulus = density * V_s^2: (Sanne)
        prem_shearmod.append(prem.prem_density(p)*pow(prem.prem_V(p)[1],2.0))

    return list_p, list_Vs, list_Vp, pv_density, fp_density, pv_shearmod, fp_shearmod, prem_shearmod


#compute prem
prem_p = numpy.arange(28.3,360.0,5)
prem_vp = [prem.prem_V(y)[0] for y in prem_p]
prem_vs = [prem.prem_V(y)[1] for y in prem_p]
prem_density = [prem.prem_density(y) for y in prem_p]

#compute murakami for 100% fp
molar_abundance=[0., 1.0]
list_p, fp_Vs, fp_Vp, pv_density, fp_density, pv_shearmod, fp_shearmod, prem_shearmod \
    = murakami(molar_abundance)

molar_abundance=[1.0, .0]
_, pv_Vs, pv_Vp, _,_,_,_,_ \
    = murakami(molar_abundance)

#molar_abundance=[0.95, .05]
molar_abundance=[0.93, .07]
_, mix_Vs, mix_Vp, _,_,_,_,_ \
    = murakami(molar_abundance)

mix_density = [molar_abundance[0] * pv_density[i] + molar_abundance[1] * fp_density[i] for i in range(len(pv_density))]

# plot Vs
pylab.subplot(2,2,1)
p1=pylab.plot(list_p,fp_Vs,'-k')
p2=pylab.plot(list_p,pv_Vs,'-b')
p3=pylab.plot(list_p,mix_Vs,'-r')
p4=pylab.plot(prem_p,prem_vs,'ok',markerfacecolor='white')
pylab.legend([p1,p2,p3,p4],["fp", "pv", "mix", "PREM"], loc=4)
pylab.title("Vs")
pylab.xlim(25,135)
pylab.ylim(5.,7.6)

# plot Vp
pylab.subplot(2,2,2)
p1=pylab.plot(list_p,fp_Vp,'-k')
p2=pylab.plot(list_p,mix_Vp,'-r')
p3=pylab.plot(prem_p,prem_vp,'ok',markerfacecolor='white')
pylab.legend([p1,p2,p3],["fp", "mix", "PREM"], loc=4)
pylab.title("Vp")
pylab.xlim(30,135)
pylab.ylim(9.25,14.)

# plot shear mod
pylab.subplot(2,2,3)
pylab.title("Shearmodulus comparison")
p1=pylab.plot(list_p,fp_shearmod,'-g')
p2=pylab.plot(list_p,pv_shearmod,'-b')
p3=pylab.plot(list_p,prem_shearmod,'ok',markerfacecolor='white',markevery=5)
pylab.legend([p1,p2,p3],["fp", "pv", "PREM"], loc=4)
pylab.xlim(30,135)

# plot density
pylab.subplot(2,2,4)
p1=pylab.plot(list_p,fp_density,'-k')
p2=pylab.plot(list_p,pv_density,'-b')
p3=pylab.plot(prem_p,prem_density,'ok',markerfacecolor='white')
p4=pylab.plot(list_p,mix_density,'-r')
pylab.legend([p1,p2,p3,p4],["fp", "pv", "PREM", "mix"], loc=4)
pylab.title("density")
pylab.xlim(30,135)
pylab.ylim(4.,6.5)



pylab.show()




test_phases()
test_mol_conv()












#print "full example:"

#inp1 = {'Mg':0.5, 'Fe': 0, 'Si':0.5, 'Ca':0.0, 'Al':0} # wt%
#inp2 = conv_inputs(inp1)
#print "in:", inp1
#print "out:", inp2

#params = {'Fe in pv': 0.0, 'Ca in pv':0.0, 'Al in pv':0.0, 'Fe in fp':0.0}
#t = determine_phases(inp2, params)
#print "phases:", t

#ret = eqn_of_state(t)
#
#print "eos:", ret
#
#print ret['density'](42)








