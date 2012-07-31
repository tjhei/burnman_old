import numpy
import bisect
import pylab

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
if __name__ == "__main__":
    p = numpy.arange(1.0,128.0,3)
    t = [geotherm(y) for y in p]
    pylab.plot(p,t,'+-')
    pylab.show()
