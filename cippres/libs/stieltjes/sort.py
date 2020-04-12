import sys


f = open(sys.argv[1],'r')
x = []
y = []
for l in f:
    w=l.split()
    x.append(float(w[0]))
    y.append(w[1:])

list1, list2 = zip(*sorted(zip(x, y)))

for i in range(len(list1)):
    print list1[i], ''.join(list2[i])
