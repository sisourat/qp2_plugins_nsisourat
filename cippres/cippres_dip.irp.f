program cippres_dip
  use general
 ! check the routine i_H_j_one_e with the variable mo_dipole_x
 ! create a routine for one-e general matrix with determinants

  implicit none
  BEGIN_DOC
! CIPPRES stands for Configuration Interaction Plugin for Photoionized and Resonant Electronic States.
! ORMAS-CI calculations are performed to compute bound and approximate continuum states. Hamiltonian and Dipole coupling matrix elements can then be computed between these states. Finally, Stieltjes imaging technique is applied to recover the correct continuum coupling matrix elements.
! cippres_dip computes the dipole matrice couplings between the CI eigenvectors of ici1 and ici2 runs
  END_DOC

  double precision :: nucl_contrib
  integer :: i, j

! TODO Read the info (ici1, ici2,...) from an input

! GENERAL
! MANU : how to get input filename from command line qp run??

! TODO Compute dipole matrix elements between two different CI runs
! TODO Include Stieltjes in qp

  print*,ifcsf
  if(ifcsf==2) then

    print*,ici1,ici2
    if(ici1==0) then
       print*, "Please set ici1 (initial states)"
       print*, "qp set cippres ici1 X "
       stop
    endif

    if(ici2==0) then
       print*, "Please set ici2 (final states)"
       print*, "qp set cippres ici2 X "
       stop
    endif

     
     if(ici1/=0 .and. ici2/=0) then
!      print*,dip_couplings_cippres
      call ezfio_set_cippres_cdip_cippres(dip_couplings_cippres)
      call ezfio_set_cippres_edip_cippres(edip_couplings_cippres)
      nucl_contrib = 0d0
      do i = 1, nucl_num
       do j = 1, 3
        nucl_contrib += nucl_charge(i)*nucl_coord(i,j)
       enddo
      enddo
      print*,'Nuclear Contribution to be added to "diagonal" elements',nucl_contrib 

      open(unit=10,file='dip.out')   
      do i = 1, n_csf_cippres(ici2) 
        write(10,*), edip_couplings_cippres(i,1), dip_couplings_cippres(i,1)
      enddo
      close(10)

!      call ezfio_set_cippres_ifcsf(3)
     endif


  else 

    print*, "ifcsf = ", ifcsf
    print*, "but it should be equal to 2 for Fano calculations"
    print*, "Please run cippres_runci first or type qp set cippres ifcsf 2"

  endif

end program cippres_dip


