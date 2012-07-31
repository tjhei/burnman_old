import numpy
import scipy.linalg
import bisect

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








