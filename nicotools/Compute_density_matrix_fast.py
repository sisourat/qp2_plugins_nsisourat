import sys
import numpy as np
from determinants import *
from numpy import linalg as LA

# READS THE DETERMINANTS OF CISTATES 1 AND 2 AND Computes the density matrix
fcistate1 = open(sys.argv[1],"r")
mo_num, ndets1, nstate1 = ( int(x) for x in fcistate1.readline().split())
print 'Number of MO', mo_num
print 'Number of Determinants', ndets1

dets_alp1 = []
dets_bet1 = []
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

# READS THE CI COEFFS
esta1 = [] 
two_e_sta1 = [] 
dip_sta = [] 
cista1 = [] 
nucl_rep = 0.0
for i in range(nstate1):
  d = fcistate1.readline().split()
  esta1.append(float(d[0])-nucl_rep)
  two_e_sta1.append(float(d[0])-float(d[1])-nucl_rep)
  dip_sta.append(float(d[2]))
  nucl_rep = float(fcistate1.readline().split()[0])
  cicoeff = [] 
  for j in range(ndets1):
    d = fcistate1.readline().split()
    cicoeff.append(float(d[0]))
  cista1.append(cicoeff)

fout = open('Density_matrices.txt','w')
feig = open('Natural_energies.txt','w')
print >> fout, mo_num, nstate1

ldets = []
for id1a in range(ndets1):
  for id1b in range(id1a+1,ndets1):
       ndiff, ndiff_alp, ndiff_bet, orbalp, orbbet = compare2dets(dets_alp1[id1a],dets_bet1[id1a],dets_alp1[id1b],dets_bet1[id1b])
       if(ndiff==2):
         if(len(orbalp)==2):
           ldets.append((id1b,id1a,orbalp))
           ldets.append((id1a,id1b,orbalp))

for i1a in range(nstate1):
 for i1b in range(i1a,nstate1):
   densmatmo = np.zeros((mo_num,mo_num))

   print >> fout, i1a, i1b
   print >> feig, i1a, i1b
   for id1 in range(ndets1):
     for imo in range(mo_num):
        if(dets_alp1[id1][imo]=='+'):
          densmatmo[imo-1][imo-1] += cista1[i1a][id1]*cista1[i1b][id1]*0.5
        if(dets_bet1[id1][imo]=='+'):
          densmatmo[imo-1][imo-1] += cista1[i1a][id1]*cista1[i1b][id1]*0.5

   for l in ldets:
     id1b, id1a, orbalp = l
     iorb1=orbalp[0]-1
     iorb2=orbalp[1]-1
     densmatmo[iorb1][iorb2] += cista1[i1a][id1a]*cista1[i1b][id1b]*0.5
     densmatmo[iorb2][iorb1] += cista1[i1a][id1a]*cista1[i1b][id1b]*0.5
   
   for imo1a in range(mo_num):
     for imo1b in range(mo_num):
       print >> fout, densmatmo[imo1a][imo1b]

   w, v = LA.eig(densmatmo)
   print >> feig, w
   print >> feig 

#   print densmatmo
   if(not(i1a==i1b)):

     print >> fout, i1b, i1a
     print >> feig, i1b, i1a

     for imo1a in range(mo_num):
       for imo1b in range(mo_num):
         print >> fout, densmatmo[imo1a][imo1b]

     print >> feig, w
     print >> feig 
