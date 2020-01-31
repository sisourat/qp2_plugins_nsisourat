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
