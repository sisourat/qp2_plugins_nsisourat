 qp_create_ezfio -b "he-cs-nico" -p bfd HeCs.xyz -m 3
 qp run scf

 qp set cippres finput_cippres cshe.xml
 qp set electrons elec_alpha_num 10
 qp set electrons elec_beta_num 9
 qp run cippres_gencsf

 qsub job_qp.sh
 qsub job_qpfano.sh

 python /home/common/quantum_package/plugins/qp2_plugins_nsisourat/cippres/nicotools/merge_fano_files.py HeCs.ezfio 188
 /home/common/quantum_package/plugins/qp2_plugins_nsisourat/cippres/libs/stieltjes/stieltjes < HeCs.ezfio/cippres/fanotot1.txt

