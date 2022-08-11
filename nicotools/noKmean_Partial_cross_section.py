import numpy as np
import sys

sig = open(sys.argv[1],"r")
dys = open(sys.argv[2],"r")
kmean = open(sys.argv[4],"r")
nprojector = int(sys.argv[3])
sip = float(sys.argv[5])
dip = float(sys.argv[6])

#dat = np.loadtxt(kmean)
#ikmean = dat[:,2]
#normDyson = dat[:,1]

sigsta = []
for l in sig:
 d = l.split()
# esta.append(float(d[1]))
 sigsta.append(float(d[2]))

esta = []
i=0
for l in dys:
 i+=1
 if(not(np.mod(i,2)==0)):
   d = l.split()
 else:
   esta.append(float(d[1]))
dys.close()
dys = open(sys.argv[2],"r")
nsta = len(esta)

i = 0
ista = 0
sigexc = 0
sigsi = 0
sigdiupper = 0
sigdi = 0
nbound = 0
nsi = 0
nsitot = 0

sigsi_sta_tot = np.zeros(nprojector)
norm_sta_tot = np.zeros(nprojector)

sigsi_sta_di = np.zeros(nprojector)
sigsi_di = 0.0

for l in dys:
 i+=1
 if(not(np.mod(i,2)==0)):
   d = l.split()
#   two_e_sta = float(d[2])
 else:
#   print(ista,esta[ista],sip)
   if(ista>nsta-1):
    break  
   d = l.split()

   if(ista>0 and esta[ista]<sip):
    sigexc += sigsta[ista]
    nbound+=1

   if(esta[ista]>sip and esta[ista]<dip):
    sigsi += sigsta[ista]
    nsi+=1
    nsitot+=1
    for j in range(nprojector):
     sigsi_sta_tot[j]+=sigsta[ista]*float(d[j])
     norm_sta_tot[j]+=float(d[j])

   if(esta[ista]>dip):
    sigdiupper+=sigsta[ista]
    nsitot+=1

    for j in range(nprojector):
     sigsi_sta_tot[j]+=sigsta[ista]*float(d[j])
     sigsi_sta_di[j]+=sigsta[ista]*float(d[j])
     norm_sta_tot[j]+=float(d[j])

    
   ista += 1

#print  sigsi_sta_tot
sigsi += np.sum(sigsi_sta_di)
sigsi_sta_tot = sigsi_sta_tot/np.sum(sigsi_sta_tot)*sigsi
#print sigexc, sigsi,  np.sum(sigsi_sta_tot), sigsi_sta_tot[0], np.sum(sigsi_sta_tot)-sigsi_sta_tot[0], sigdiupper-np.sum(sigsi_sta_di), sigdiupper
#print sigdiupper-np.sum(sigsi_sta_di), sigdiupper
print  sigsi, sigsi_sta_tot[0]



