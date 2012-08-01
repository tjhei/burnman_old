import operator

def sort_table(table, col=0):
    return sorted(table, key=operator.itemgetter(col))


# check if two floats are close to each other
def float_eq(a,b):
    return abs(a-b)<1e-10*max(1e-5,abs(a),abs(b))


def linear_interpol(x,x1,x2,y1,y2):
    print x,x1,x2
    assert(x1<=x)
    assert(x2>=x)
    assert(x1<=x2)

    alpha = (x - x1) / (x2-x1)
    return (1.-alpha)*y1 + alpha*y2
