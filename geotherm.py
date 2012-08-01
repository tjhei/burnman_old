import numpy
import bisect
import pylab
from tools import *

# loads a simple geotherm from geotherm.txt


# pressure: in GPa
# return: temperature in K
def geotherm(pressure):
    idx = bisect.bisect_left(table_p, pressure) - 1
    if (idx < 0):
        return table_T[0]
    elif (idx < len(table_p)-1):
        return linear_interpol(pressure, table_p[idx], table_p[idx+1], table_T[idx], table_T[idx+1])
    else:
        return table_T[idx]


# polynomial fit from Watson, Baxter, EPSL, 2007
# pressure: in GPa
# return: temperature in K
def geotherm_formula(pressure):
    if (pressure <= 15):
        return 1900-1420*pow(0.8,pressure)
    else:
        return 1680+11.1*pressure




geotherm_table=[]
for line in open("geotherm.txt").readlines():
    if (line[0]!='#'):
	numbers = map(float, line.split())
	geotherm_table.append(numbers)

geotherm_table = sort_table(geotherm_table, 0)
table_p=numpy.array(geotherm_table)[:,0]
table_T=numpy.array(geotherm_table)[:,1]




# test geotherm
if __name__ == "__main__":
    print table_p
    p = numpy.arange(1.0,128.0,3)
    t = [geotherm(y) for y in p]
    t2 = [geotherm_formula(y) for y in p]
    pylab.plot(p,t,'+-')
    pylab.plot(p,t2,'x--r')
    pylab.show()
