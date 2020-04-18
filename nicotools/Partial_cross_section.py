import numpy as np
import sys

sig = open(sys.argv[1],"r")
dys = open(sys.argv[2],"r")
kmean = open(sys.argv[4],"r")
nprojector = int(sys.argv[3])

dat = np.loadtxt(kmean)
ikmean = dat[:,2]
normDyson = dat[:,1]

esta = []
sigsta = []
for l in sig:
 d = l.split()
 esta.append(float(d[1]))
 sigsta.append(float(d[2]))

nsta = len(esta)

i = 0
ista = 0
sigexc = 0
sigsitot = 0
sigditot = 0
sigdi = 0
normtot = 0 
nbound = 0
nsi = 0
nsitot = 0

for l in dys:
 i+=1
 if(not(np.mod(i,2)==0)):
   d = l.split()
#   two_e_sta = float(d[2])
 else:
   if(ista>nsta-1):
    break  
   d = l.split()

   if(ista>0 and esta[ista]<-2.0):
    sigexc += sigsta[ista]
    nbound+=1

   if(esta[ista]>-2.0 and esta[ista]<0):
    sigsitot += sigsta[ista]
    nsi+=1
    nsitot+=1

   if(esta[ista]>0):
    norm = 0
    for j in range(nprojector):
     norm += float(d[j])
    normtot += norm
#    if ikmean[ista-(nbound+nsi+1)]==1:
    if normDyson[ista-(nbound+nsi+1)]>0.3:
     sigsitot += sigsta[ista] 
     nsitot+=1
    else :
     sigdi += sigsta[ista] 
   ista += 1

print sigexc, sigsitot, sigdi
#print normtot/nsta
