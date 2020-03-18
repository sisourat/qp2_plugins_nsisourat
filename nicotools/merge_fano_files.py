import sys
import gzip

fdir = sys.argv[1]
nfiles = int(sys.argv[2])

ftot = open(fdir+'/cippres/fanotot.txt','w')

print >> ftot, '1'
for i in range(nfiles):
 fcinp = gzip.open(fdir+'/cippres/cfano_cippres'+str(i+1)+'.gz','r')

 feinp = gzip.open(fdir+'/cippres/efano_cippres'+str(i+1)+'.gz','r')
 ltmp = feinp.readline()
 ltmp = feinp.readline()

 ntmp = fcinp.readline()
 l = fcinp.readline()
 d = l.split()
 nsta = int(d[0])

 if(i==0):
  print >> ftot, nsta*nfiles

 for j in range(nsta):
   lc = fcinp.readline()
   dc = lc.split()
   le = feinp.readline()
   de = le.split()
   print >> ftot, de[0], dc[0]
