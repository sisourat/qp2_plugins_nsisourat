source /home/common/quantum_package/mendeleiev_quantum_package.rc

qp create_ezfio -b "he-test" LiHe.xyz -o LiHe.ezfio -m 2
qp set electrons elec_alpha_num 2
qp set electrons elec_beta_num 2

qp run scf
qp set cippres finput_cippres lihe.xml

qp set electrons elec_alpha_num 3
qp run cippres_gencsf

qsub job_qp.sh

## i_stj_job=1 Fano; i_stj_job=2 Dipole
qp set cippres i_stj_job 1  
## ifanosta=index of the resonance state (from ici1 run)
qp set cippres ifanosta 1
qp run cippres_stieltjes

