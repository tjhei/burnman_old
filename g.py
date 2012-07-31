
import numpy
import scipy.linalg



#print X
#    A = numpy.matrix([[ 1, 0, 1],
#    [ 0,1,0],
#    [ 0,1,1]])
#    b = numpy.matrix([[1],[2],[3]])
    
#    x = scipy.linalg.solve(A, b)
    #print x




#
# inp = {'MgO':beta, 'FeO': 0, 'SiO2': gamma, 'CaO':0, 'Al2O3':0.1}
# params = {'Fe in pv': 0.9, 'Ca in pv':0.1, 'Al in pv', 'Fe in fp':1.0}
# ret = 'mol pv' A, 'mol fp' B, 'mol st' C
#, 'Mg in pv':0, 'Fe in pv':0, 'Ca in pv':0,'Si in pv':0, 'Al in pv':0, 'Mg in fp':0,'Fe in fp':0,}
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



# test everything into pv
inp = {'MgO':20, 'FeO': 0, 'SiO2':20, 'CaO':0, 'Al2O3':0.0}
params = {'Fe in pv': 0.0, 'Ca in pv':0.0, 'Al in pv':0.0, 'Fe in fp':0.0}
t = determine_phases(inp, params)
assert t['mol pv'] == 20
assert t['mol fp'] == 0
assert t['mol st'] == 0






