=======
cippres
=======

Authors : Nico/Manu
Date : 12/09/2019

CIPPRES stands for Configuration Interaction Plugin for Photoionized and Resonant Electronic States.

ORMAS-CI calculations are performed to compute bound and approximate continuum states. Hamiltonian and Dipole coupling matrix elements can then be computed between these states. Finally, Stieltjes imaging technique is applied to recover the correct continuum coupling matrix elements.

REQUIRMENTS:
* QP (obviously)
* python
* python libs : xml.dom, itertools, collections

USERS GUIDE:

0) (create or set EZFIO database prior to cipres runs)
qp create_ezfio -b "aug-cc-pvtz" He.xyz -o He.ezfio

or 

~/Workspace/qp2/bin/qpsh
qp set_file He.ezfio

0') (run scf prior to cipres runs)
qp run scf

0'') (2e integrals AO to MO transformation, not mandatory but highly recommanded for better performance)
 qp run four_idx_transform

1) set the xml input file in EZFIO
 qp set cippres finput_cippres test.xml

2) generate CSFs in header.txt and list.txt
 qp run cippres_gencsf
 (maybe you will have to qp set electrons elec_beta_num X)

3)  run CI (if ifcsf 1)
 qp run cippres_runci

4)  run Fano (if ifcsf 2)
 qp set cippres ici1 X
 qp set cippres ici2 X
 qp run cippres_fano

4')  run dip (if ifcsf 2)
 qp set cippres ici1 X
 qp set cippres ici2 X
 qp run cippres_dip

5)  run Stieltjes (in libs)
 $QP_ROOT/plugins/qp_plugins_sisourat-/libs/stieltjes/stieltjes < fano.txt

TESTS:







