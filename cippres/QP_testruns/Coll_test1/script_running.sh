source /home/common/quantum_package/mendeleiev_quantum_package.rc

qp create_ezfio -b "aug-cc-pvtz" He.xyz -o He.ezfio
qp run scf
qp set cippres finput_cippres he.xml

qp run cippres_gencsf
qp run cippres_runci

qp set cippres finput_coll coll_input.xml
qp set cippres ici1 1
qp run cippres_setup_collision
qp run cippres_collision
qp run cippres_prop_collision

