import numpy as np
import sys

sig = open(sys.argv[1],"r")
dys = open(sys.argv[2],"r")
kmean = open(sys.argv[4],"r")
nprojector = int(sys.argv[3])
sip = float(sys.argv[5])
dip = float(sys.argv[6])

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
sigdiupper = 0
sigdi = 0
normtot = 0 
nbound = 0
nsi = 0
nsitot = 0

sigsi_sta_tot = np.zeros(nprojector)
norm_sta_tot = np.zeros(nprojector)

for l in dys:
 i+=1
 if(not(np.mod(i,2)==0)):
   d = l.split()
#   two_e_sta = float(d[2])
 else:
   if(ista>nsta-1):
    break  
   d = l.split()

   if(ista>0 and esta[ista]<sip):
    sigexc += sigsta[ista]
    nbound+=1

   if(esta[ista]>sip and esta[ista]<dip):
    sigsitot += sigsta[ista]
    nsi+=1
    nsitot+=1
    normtot = 0.0
    norm_sta = np.zeros(nprojector)
    sigsi_sta = np.zeros(nprojector)
    for j in range(nprojector):
     sigsi_sta[j]+=sigsta[ista]*float(d[j])
     norm_sta[j]=float(d[j])
     normtot += float(d[j])
    norm_sta_tot += norm_sta
    sigsi_sta_tot += sigsi_sta/np.sum(norm_sta)

   if(esta[ista]>dip):
    sigdiupper+=sigsta[ista]
    if ikmean[ista-(nbound+nsi+1)]==0:
#    if normDyson[ista-(nbound+nsi+1)]>0.01:
     sigsitot += sigsta[ista] 
     nsitot+=1

     normtot = 0.0
     norm_sta = np.zeros(nprojector)
     sigsi_sta = np.zeros(nprojector)
     for j in range(nprojector):
      sigsi_sta[j]+=sigsta[ista]*float(d[j])
      norm_sta[j]=float(d[j])
#      normtot += float(d[j])
#     normtot += float(d[-1])
     sigsi_sta_tot += sigsi_sta/np.sum(norm_sta)
    
     norm_sta_tot += norm_sta
    else :
      sigdi += sigsta[ista] 
   ista += 1

#for i in range(nprojector):
#  print i, sigsi_sta_tot[i]#/np.sum(sigsi_sta_tot)*sigsitot, norm_sta_tot[i]/np.sum(norm_sta_tot), norm_sta_tot[i]/np.sum(norm_sta_tot)*sigsitot

#sigsi_sta_tot= sigsi_sta_tot/np.sum(sigsi_sta_tot)*sigsitot
#norm_sta_tot=norm_sta_tot/np.sum(norm_sta_tot)
#print np.sum(norm_sta_tot)

#print 
#print np.sum(sigsi_sta_tot), sigsitot
#print np.sum(sigsi_sta_tot[0:3])
#psys.exit()
#print "1s", np.sum(sigsi_sta_tot[0]),
#print "2s", np.sum(sigsi_sta_tot[1]),
#print "2p", np.sum(sigsi_sta_tot[2:5]),
#print "3s", np.sum(sigsi_sta_tot[5]),
#print "3p", np.sum(sigsi_sta_tot[6:9]),
#print "n>3", np.sum(sigsi_sta_tot[10:nprojector+1])

#print norm_sta_tot

#print np.sum(sigsi_sta_tot[0]),
#print np.sum(sigsi_sta_tot[1]),
#print np.sum(sigsi_sta_tot[2:5]),
#print np.sum(sigsi_sta_tot[5]),
#print np.sum(sigsi_sta_tot[6:9]),
#print np.sum(sigsi_sta_tot)

print 
print sigexc, sigsitot, sigdi, sigdiupper
#print normtot/nsta
