import sys
import gzip

fdir = sys.argv[1]
nfiles = int(sys.argv[2])


for i in range(nfiles):
 fcinp = gzip.open(fdir+'/cippres/cfano_cippres'+str(i+1)+'.gz','r')

 feinp = gzip.open(fdir+'/cippres/efano_cippres'+str(i+1)+'.gz','r')
 ltmp = feinp.readline()
 ltmp = feinp.readline()

 ntmp = fcinp.readline()
 l = fcinp.readline()
 d = l.split()
 nsta_d = int(d[0])
 nsta_f = int(d[1])

 for j in range(nsta_d):
  ftot = open(fdir+'/cippres/fanotot'+str(j)+'.txt','a')
  if(i==0):
    print >> ftot, '1'
    print >> ftot, nsta_d*nfiles
  for k in range(nsta_f):
    lc = fcinp.readline()
    dc = lc.split()
    le = feinp.readline()
    de = le.split()
    print >> ftot, de[0], dc[0]
  ftot.close()
