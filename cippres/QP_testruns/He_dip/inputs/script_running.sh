source /home/common/quantum_package/mendeleiev_quantum_package.rc

qp create_ezfio -b "aug-cc-pvtz" He.xyz -o He.ezfio
qp run scf
qp set cippres finput_cippres he.xml

qp run cippres_gencsf
qp run cippres_runci

qp set cippres ici1 1
qp set cippres ici2 2
qp run cippres_dip

