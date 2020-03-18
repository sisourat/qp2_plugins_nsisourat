source /home/common/quantum_package/mendeleiev_quantum_package.rc

qp create_ezfio -b "ne-eventmp" Ne.xyz -o Ne.ezfio
qp run scf 

qp run swap_mos << EOF
 1 
 5
EOF

qp set electrons elec_alpha_num 5
qp set electrons elec_beta_num 4 

qp set_mo_class -a [1-4,6-109] -c [5]

qp set scf_utils frozen_orb_scf True
qp run scf 


qp set cippres finput_cippres ne.xml
qp run cippres_gencsf

qsub job_qp.sh
qsub job_qpfano.sh
python /home/common/quantum_package/plugins/qp2_plugins_nsisourat/cippres/nicotools/merge_fano_files.py Ne.ezfio 50
/home/common/quantum_package/plugins/qp2_plugins_nsisourat/cippres/libs/stieltjes/stieltjes < Ne.ezfio/cippres/fanotot.txt

## then one can analyze stieltjes.order.* files


