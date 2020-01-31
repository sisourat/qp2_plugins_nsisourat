program cippres_runci
  use general
  implicit none
  BEGIN_DOC
! CIPPRES stands for Configuration Interaction Plugin for Photoionized and Resonant Electronic States.
! ORMAS-CI calculations are performed to compute bound and approximate continuum states. Hamiltonian and Dipole coupling matrix elements can then be computed between these states. Finally, Stieltjes imaging technique is applied to recover the correct continuum coupling matrix elements.
! cippres_runci performs ORMAS-CI calculations using lists of CSFs generated previously (with cippres_gencsf)
  END_DOC

  character(len=lenmax) :: finput
  integer :: ilen, jlen
  logical :: file_e

  integer :: i, j

! TODO Read lists of CSFs from EZFIO (see cippres_gencsf first)


! GENERAL
! TODO Compute dipole matrix elements between two different CI runs
! TODO Include Stieltjes in qp

  if(ifcsf==1) then

    do i = 1, n_ciruns_cippres
!      print*,'ncsfs',i, n_csf_cippres(i) 
      eigvalues_cippres(1:n_csf_cippres(i),i) += nuclear_repulsion
      print*,'CI eigval =', eigvalues_cippres(1:n_csf_cippres(i),i)
    enddo

   call ezfio_set_cippres_eigvalues_cippres(eigvalues_cippres)
   call ezfio_set_cippres_eigvectors_cippres(eigvectors_cippres)
   call ezfio_set_cippres_ifcsf(2)

  else

   print*, "Please run cippres_gencsf first or type qp set cippres ifcsf 1"
   print*, "Note that if you rerun cippres_runci, you should first delete the eigval/vec from EZFIO"

  endif

end program cippres_runci
