source /home/common/quantum_package/mendeleiev_quantum_package.rc

rm -f fanotot.txt
for i in {2..71}
do
 rm Ne.ezfio/cippres/*fano*
 qp set cippres ici2 $i
 qp get cippres ici2
 qp run cippres_fano
 cat fano.out >> fanotot.txt
done
