import numpy as np
import sys

from sklearn.cluster import KMeans
import matplotlib.pyplot as plt

dys = open(sys.argv[1],"r")
nbound = int(sys.argv[2])

fout = open("features_kmean.txt","w")

i = 0
ista = 0
normtot = 0 
esta = []
for l in dys:
 i+=1
 if(not(np.mod(i,2)==0)):
   d = l.split()
   esta.append(float(d[1]))
   two_e_sta = float(d[2])
 else:
   d = l.split()

   norm = 0
   for j in range(nbound):
     norm += float(d[j])
   if(esta[ista]>0):
    print(two_e_sta, norm, esta[ista],file=fout)
   ista+=1

fout.close()

plt.figure(figsize=(12, 12))
plt.xlabel('e-e repulsion')
plt.ylabel('Sum of the norm of the Dyson orbitals')

dat=np.loadtxt("features_kmean.txt")
X = dat[:,0:2]

random_state = 170
y_pred = KMeans(n_clusters=2, random_state=random_state).fit_predict(X)

for i in range(len(y_pred)):
        print(X[i,0],X[i,1],y_pred[i])

plt.scatter(X[:, 0], X[:, 1], c=y_pred)

plt.show()

