source /home/common/quantum_package/mendeleiev_quantum_package.rc

qp create_ezfio -b "cc-pvdz" He.xyz -o He.ezfio
qp run scf
qp set cippres finput_cippres he.xml

qp run cippres_gencsf
qp run cippres_runci
qp run cistate_analysis
# cp cistate in a file

qp create_ezfio -b "cc-pvdz" He.xyz -o Hep.ezfio
qp set electrons elec_alpha_num 1
qp set electrons elec_beta_num 0
qp run scf
qp set cippres finput_cippres hep.xml

qp run cippres_gencsf
qp run cippres_runci
qp run cistate_analysis
# cp cistate in a file

python /home/nico/Workspace/qp2/plugins/qp2_plugins_nsisourat/nicotools/Compute_Dyson.py det_N det_Nm1






