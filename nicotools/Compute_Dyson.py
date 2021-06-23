import sys
import numpy as np
from determinants import *

nstate_dyson = int(sys.argv[3])

# READS THE DETERMINANTS OF CISTATES 1 AND 2 AND LIST THOSE DIFFERING BY ONLY ONE ORBITAL (THE ONE "MISSING")
fcistate1 = open(sys.argv[1],"r")
mo_num, ndets1, nstate1 = ( int(x) for x in fcistate1.readline().split())
if(nstate_dyson>nstate1):
 print "nstate_dyson should be smaller than nstate1"
 sys.exit()
print 'Number of MO', mo_num
print 'Number of Determinants', ndets1

dets_alp1 = []
dets_bet1 = []
for i in range(ndets1):
  fcistate1.readline()
  det = [' ']
  l = fcistate1.readline().split()
  det += [ x for x in l[0] if (x=='+' or x=='-')]
  dets_alp1.append(det)
  det = [' ']
  l = fcistate1.readline().split()
  det += [ x for x in l[0] if (x=='+' or x=='-')]
  dets_bet1.append(det)

fcistate2 = open(sys.argv[2],"r")
mo_num, ndets2, nstate2 = ( int(x) for x in fcistate2.readline().split())
print 'Number of MO', mo_num
print 'Number of Determinants', ndets2
print 

dets_alp2 = []
dets_bet2 = []
for i in range(ndets2):
  fcistate2.readline()
  det = [' ']
  l = fcistate2.readline().split()
  det += [ x for x in l[0] if (x=='+' or x=='-')]
  dets_alp2.append(det)
  det = [' ']
  l = fcistate2.readline().split()
  det += [ x for x in l[0] if (x=='+' or x=='-')]
  dets_bet2.append(det)

ldets = []
for i in range(ndets1):
  for j in range(ndets2):
   # print dets_alp[i],dets_alp[j]
    ndiff, ndiff_alp, ndiff_bet, orbdiff_alp, orbdiff_bet = compare2dets(dets_alp1[i],dets_bet1[i],dets_alp2[j],dets_bet2[j])
   # print j,i, ndiff
    if(ndiff==1):
       ldets.append((j,i,orbdiff_alp,orbdiff_bet))

#print ldets

# READS THE CI COEFFS
esta1 = [] 
two_e_sta1 = [] 
dip_sta = [] 
cista1 = [] 
nucl_rep = 0.0
for i in range(nstate1):
  d = fcistate1.readline().split()
  nucl_rep = float(fcistate1.readline().split()[0])
  esta1.append(float(d[0])-nucl_rep)
  two_e_sta1.append(float(d[0])-float(d[1])-nucl_rep)
  dip_sta.append(float(d[2]))
  cicoeff = [] 
  for j in range(ndets1):
    cicoeff.append(float(fcistate1.readline().split()[0]))
  cista1.append(cicoeff)

#print esta1
#print cista1

esta2 = []
cista2 = []
for i in range(nstate2):
  d = fcistate2.readline().split()
  nucl_rep = float(fcistate2.readline().split()[0])
  esta2.append(float(d[0])-nucl_rep)
  cicoeff = []
  for j in range(ndets2):
    cicoeff.append(float(fcistate2.readline().split()[0]))
  cista2.append(cicoeff)

#print esta2
#print cista2

# COMPUTES NORM OF THE DYSON ORBITALS
fdyson = open('Dyson_norms.txt','w')
for i1 in range(nstate_dyson):
  normtot = 0.0
  print >> fdyson, i1, esta1[i1],two_e_sta1[i1], dip_sta[i1]
  for i2 in range(nstate2):
#    print "(N-1)e & (N)e states = ",i2,i1
    mocoeffs = np.zeros(mo_num+1)
    for l in ldets:
      j2, j1, orbalp, orbbet = l
#      print j2, j1, orbalp, orbbet
#!      print cista1[i1][j1], cista2[i2][j2]
      if(len(orbalp)==1):
        orb=orbalp[0]
      else:
        orb=orbbet[0]
      mocoeffs[orb]+=cista1[i1][j1]*cista2[i2][j2]
#      c1 = cista1[i1]
    normtot+=np.sum(np.square(mocoeffs))
    print >> fdyson, np.sum(np.square(mocoeffs)),
#    print 'Dyson orb. norm',np.sum(np.square(mocoeffs))
#    print 
  print >> fdyson, normtot
#  print normtot
#  print 
print "Dyson norms in Dyson_norms.txt"




