import sys 
import numpy as np
from xml.dom import minidom
from spineigenfunctions import *
import itertools
from collections import Counter

###################################################################

# finp is the xml file
finp = sys.argv[1]
doc = minidom.parse(finp)

coll_input = doc.getElementsByTagName("CollInput")
b = coll_input[0].getElementsByTagName('ImpactParam')
btype = b[0].getAttribute("type")
bmin = b[0].getAttribute("bmin")
bmax = b[0].getAttribute("bmax")
nb = b[0].getAttribute("nb")
print btype, bmin, bmax, nb

vp = coll_input[0].getElementsByTagName('ImpactVel')
#vx = vp[0].getAttribute("vx")
#vy = vp[0].getAttribute("vy")
vz = vp[0].getAttribute("vz")
print vz

zgrid = coll_input[0].getElementsByTagName('Zgrid')
ztype = zgrid[0].getAttribute("type")
zmin = zgrid[0].getAttribute("zmin")
zmax = zgrid[0].getAttribute("zmax")
nzgrid = zgrid[0].getAttribute("nzgrid")
print ztype, zmin, zmax, nzgrid

states = coll_input[0].getElementsByTagName('Bound')
stamin = states[0].getAttribute("stamin")
stamax = states[0].getAttribute("stamax")
print stamin, stamax

states = coll_input[0].getElementsByTagName('SingIon')
stamin = states[0].getAttribute("stamin")
stamax = states[0].getAttribute("stamax")
print stamin, stamax

states = coll_input[0].getElementsByTagName('DoubIon')
stamin = states[0].getAttribute("stamin")
stamax = states[0].getAttribute("stamax")
print stamin, stamax

state = coll_input[0].getElementsByTagName('InitState')
istate = state[0].getAttribute("state")
print istate

potentials = coll_input[0].getElementsByTagName('Potential')
print len(potentials)
for pot in potentials:
    print pot.getAttribute("charge")


sys.exit()

