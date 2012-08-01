import numpy
import bisect
import pylab
from tools import *

# this loads the PREM seismic velocities from prem_table.txt

# pressure: in GPa
# return: V_p, V_s in km/s
def prem_V(pressure):
    idx = bisect.bisect_left(table_p, pressure) - 1

    if (idx < 0):
        return table[0][3], table[0][4]
    elif (idx < len(table_p)-1):
        return linear_interpol(pressure, table_p[idx], table_p[idx+1], table[idx][3], table[idx+1][3]), \
            linear_interpol(pressure, table_p[idx], table_p[idx+1], table[idx][4], table[idx+1][4])
    else:
        return table[idx][3], table[idx][4]


def prem_radius(pressure):
    idx = bisect.bisect_left(table_p, pressure) - 1

    if (idx < 0):
        return table[0][0]
    elif (idx < len(table_p)-1):
        return linear_interpol(pressure, table_p[idx], table_p[idx+1], table[idx][0], table[idx+1][0])
    else:
        return table[idx][0]


#radius pressure density V_p V_s
table=[] 

for line in open("prem_table.txt").readlines():
    if (line[0]!='#'):
	numbers = map(float, line.split())
        numbers[1]=numbers[1]*0.1 # convert kbar to GPa
	table.append(numbers)

table = sort_table(table, 1)

table_p=numpy.array(table)[:,1]





# test
if __name__ == "__main__":
    p = numpy.arange(1.0,360.0,3)
    vs = [prem_V(y)[0] for y in p]
    pylab.plot(p,vs,'+-')
    pylab.show()
