import sys
import numpy as np
from determinants import *


# arg 1 file with N states
# arg 2 file with N-1 states
# arg 3 nstate_dyson
# arg 4 Psit_collision.out file

fpsit =open(sys.argv[4],"r")

nstate_dyson = int(sys.argv[3])
nbb, nsta, ntime = ( int(x) for x in fpsit.readline().split())

esta = []
psibt = []
for i in range(nsta):
 esta.append(float(fpsit.readline().split()[0]))
for i in range(nbb):
  psib = []
  for j in range(ntime):
    psi = [float(x) for x in fpsit.readline().split()]
    psib.append(psi)
  psibt.append(psib)
 

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
dets_typ1 = []
for i in range(ndets1):
  fcistate1.readline()
  deta = [' ']
  l = fcistate1.readline().split()
  deta += [ x for x in l[0] if (x=='+' or x=='-')]
  nmo_occ = deta.count('+')
  dets_alp1.append(deta)
  detb = [' ']
  l = fcistate1.readline().split()
  detb += [ x for x in l[0] if (x=='+' or x=='-')]
  dets_bet1.append(detb)
  nholes=0
  for i in range(nmo_occ):
   if(not(deta[i+1])=='+'):
     nholes+=1
   if(not(detb[i+1])=='+'):
     nholes+=1
  dets_typ1.append(nholes)

fcistate2 = open(sys.argv[2],"r")
mo_num, ndets2, nstate2 = ( int(x) for x in fcistate2.readline().split())
print 'Number of MO', mo_num
print 'Number of Determinants', ndets2
print 

dets_alp2 = []
dets_bet2 = []
dets_typ2 = []
for i in range(ndets2):
  fcistate2.readline()
  deta = [' ']
  l = fcistate2.readline().split()
  deta += [ x for x in l[0] if (x=='+' or x=='-')]
  dets_alp2.append(deta)
  detb = [' ']
  l = fcistate2.readline().split()
  detb += [ x for x in l[0] if (x=='+' or x=='-')]
  dets_bet2.append(detb)
  nholes=0
  for i in range(nmo_occ):
   if(not(deta[i+1])=='+'):
     nholes+=1
   if(not(detb[i+1])=='+'):
     nholes+=1
#   if(nholes==2 and deta[0:nmo_occ+2]==detb[0:nmo_occ+2]):
#     nholes+=1
  dets_typ2.append(nholes)

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
w_singles = []
w_doubles = []
for i in range(nstate1):
  d = fcistate1.readline().split()
  nucl_rep = float(fcistate1.readline().split()[0])
  esta1.append(float(d[0])-nucl_rep)
  two_e_sta1.append(float(d[0])-float(d[1])-nucl_rep)
  dip_sta.append(float(d[2]))
  ws = 0.0
  wd = 0.0
  cicoeff = [] 
  for j in range(ndets1):
    d = fcistate1.readline().split()
    cicoeff.append(float(d[0]))
    if(dets_typ1[j]==1):
      ws += float(d[0])**2
    elif(dets_typ1[j]==2):
      wd += float(d[0])**2
  w_singles.append(ws)
  w_doubles.append(wd)
  cista1.append(cicoeff)

if(not(nsta==nstate1)):
  print 'nsta should be equal to nstate1'
  sys.exit()

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

modys_perstate = []
for i1 in range(nstate1):
  modys = []
  for i2 in range(nstate2):
    mocoeffs = np.zeros(mo_num+1)
    for l in ldets:
      j2, j1, orbalp, orbbet = l
      if(len(orbalp)==1):
        orb=orbalp[0]
      else:
        orb=orbbet[0]
      mocoeffs[orb]+=cista1[i1][j1]*cista2[i2][j2]
    modys.append(mocoeffs)
  modys_perstate.append(modys)

print len(modys_perstate),len(modys_perstate[0]),len(modys_perstate[0][0])

#k2 = 0
#k = 2
#ksta = k/2-1
#l = 1
#print (psibt[0][0][k]+psibt[0][0][k+1]*1j)
#print modys_perstate[ksta][k2][l]*(psibt[0][0][k]+psibt[0][0][k+1]*1j)
#sys.exit()

# COMPUTES TIME DEP Dyson Norms

for i in range(nbb):
 foutb = open('tdnorm'+str(i)+'.dat','w')
 for j in range(ntime):
   bimp = psibt[i][j][0]
   time = psibt[i][j][1]
   tdnorm = []
   for k2 in range(nstate2):
     tdmocoeffs = np.zeros(mo_num+1,dtype=complex)
     for k in range(2,2+2*nsta,2):
       ksta = k/2-1
       for l in range(1,mo_num+1):
         tdmocoeffs[l] += (psibt[i][j][k]+psibt[i][j][k+1]*1j)*modys_perstate[ksta][k2][l]
     tdnorm.append(np.sum( np.abs(tdmocoeffs)**2 ))

#     print np.abs(tdmocoeffs)**2,np.sum( np.abs(tdmocoeffs)**2 )
   print >> foutb,  time,' '.join(map(str, tdnorm)), np.sum(tdnorm)
 foutb.close()
# sys.exit()
      
 
#    print bimp,time, psibt[i][j][k], psibt[i][j][k+1]



sys.exit()


#      c1 = cista1[i1]
#    normtot+=np.sum(np.square(mocoeffs))
#    print >> fdyson, np.sum(np.square(mocoeffs)),
#    print 'Dyson orb. norm',np.sum(np.square(mocoeffs))
#    print 
#  print >> fdyson, normtot
#  print normtot
#  print 
#print "Dyson norms in Dyson_norms.txt"




