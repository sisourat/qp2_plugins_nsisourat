import numpy as np
import sys

fprob = open(sys.argv[1],"r")
fdys = open(sys.argv[2],"r")
nprojector = int(sys.argv[3])
sip = float(sys.argv[4])
dip = float(sys.argv[5])


w = fprob.readline().split()
nbb = int(w[0])
nsta = int(w[1])

esta = np.zeros(nsta)
for i in range(nsta):
 esta[i] = float(fprob.readline().split()[0])

i=0
cdys = np.zeros(nsta)
for l in fdys:
 i+=1
 if((np.mod(i,2)==0)):
   d = l.split()
   norm = 0.0
   for j in range(nprojector):
    norm += float(d[j])
   cdys[(i-1)/2] = norm
fdys.close()



prob = np.zeros(nsta)
si = 0.0
di = 0.0
siupp = 0.0
diupp = 0.0
for i in range(nbb):
 p0 = 0.0
 p1 = 0.0
 p1upp = 0.0
 p2 = 0.0
 p2upp = 0.0
 w = fprob.readline().split()
 for j in range(1,nsta+1):
   k = j -1
   prob[k] = float(w[j])
   if(esta[k]>sip and esta[k]<dip):
     p1 += prob[k]
     p1upp += prob[k] 
   elif(esta[k]>dip):
     p1 += prob[k]*cdys[k]*(2.0-cdys[k]) 
     p2 += (1.0-cdys[k])**2*prob[k] 
#     print prob[k], prob[k]*cdys[k]*(2.0-cdys[k])+(1.0-cdys[k])**2*prob[k]
     p2upp += prob[k] 
   else:
     p0 += prob[k]

 siupp +=float(w[0])*p1upp*2.0*np.pi*0.28
 si +=float(w[0])*p1*2.0*np.pi*0.28
 diupp +=float(w[0])*p2upp*2.0*np.pi*0.28
 di +=float(w[0])*p2*2.0*np.pi*0.28
# print  w[0], float(w[0])*p1*2.0*np.pi*0.28, float(w[0])*p2*2.0*np.pi*0.28
 print w[0],p0,p1,p2,2*p0*p2-0.5*p1**2,"     ",p0,p1upp,p2upp,2*p0*p2upp-0.5*p1upp**2#,p0+p1+p2, p0+p1upp+p2upp
 
print si*0.2, di*0.2,siupp*0.2,diupp*0.2
fprob.close()

