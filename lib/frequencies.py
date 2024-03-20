
import lib.my_constants as const

def PlasmaFreq(n,Z,mass):
    return (n*((Z*const.qe)**2)/(mass*const.e0))**.5

def CyclotronFreq(B,Z,mass):
    return B*Z*const.qe/mass
