import numpy as np
import sys

sig = open(sys.argv[1],"r")
dys = open(sys.argv[2],"r")

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
for l in dys:
 i+=1
 if(np.mod(i,2)==0):
   if(ista>nsta):
    break  
   d = l.split()

   if(ista>0 and esta[ista]<-2.0):
    sigexc += sigsta[ista]

   if(esta[ista]>-2.0 and esta[ista]<0):
    sigsitot += sigsta[ista]

   if(esta[ista]>0):
    sigsi = 0
    for j in range(19):
     sigsi += sigsta[ista]*float(d[j])
    sigsitot += sigsi
    sigditot += sigsta[ista]
    sigdi += sigsta[ista]-sigsi
#    print esta[ista], sigsta[ista], sigsi, sigsta[ista]-sigsi
   ista += 1

print
print
print sigexc, sigsitot, sigdi
