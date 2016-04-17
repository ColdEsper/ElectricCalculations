import math

def discCase (q, r, z):
    sigma = q/(math.pi*r*r)
    epsilon=8.854187817*math.pow(10,-12)
    return (sigma/(2*epsilon))*(1-1/math.sqrt((r/z)*(r/z)+1))
