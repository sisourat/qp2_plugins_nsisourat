source /home/common/quantum_package/mendeleiev_quantum_package.rc

qp create_ezfio -b "he-test" LiHe.xyz -o LiHe.ezfio -m 2
qp set electrons elec_alpha_num 2
qp set electrons elec_beta_num 2

qp run scf
qp set cippres finput_cippres lihe.xml

qp set electrons elec_alpha_num 3
qp run cippres_gencsf

qsub job_qp.sh

OR


qp run cippres_runci

qp set cippres ici1 1
qp set cippres ici2 2
qp run cippres_fano

run stieltjes

