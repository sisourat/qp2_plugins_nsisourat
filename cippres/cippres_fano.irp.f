program cippres_fano
  use general
  implicit none
  BEGIN_DOC
! CIPPRES stands for Configuration Interaction Plugin for Photoionized and Resonant Electronic States.
! ORMAS-CI calculations are performed to compute bound and approximate continuum states. Hamiltonian and Dipole coupling matrix elements can then be computed between these states. Finally, Stieltjes imaging technique is applied to recover the correct continuum coupling matrix elements.
! cippres_fano computes the H matrice couplings between the CI eigenvectors of ici1 and ici2 runs
  END_DOC

 integer :: i, j

! TODO Read the info (ici1, ici2,...) from an input

! GENERAL
! MANU : how to get input filename from command line qp run??

! TODO Compute dipole matrix elements between two different CI runs
! TODO Include Stieltjes in qp

  if(ifcsf==2) then

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

!    if(ifanosta==0) then
!       print*, "Please set ifanosta (the initial state)"
!       print*, "qp set cippres ifanosta X "
!       stop
!    endif

!      call ezfio_get_cippres_ici1(ici1)
!      call ezfio_get_cippres_ici2(ici2)
     if(ici1/=0 .and. ici2/=0) then
!      print*,ici1,ici2
!      print*, twoe_couplings_cippres(:,:)
!      print*, e_couplings_cippres(:,:)
 
!      open(unit=10,file='fano.out')   
!      do i = 1, n_csf_cippres(ici2) 
!        write(*,'(100(e24.16,1X))')  (e_couplings_cippres(i,j), twoe_couplings_cippres(i,j), j=ifanosta,ifanosta+10)
!        write(10,'(100(e24.16,1X))') (e_couplings_cippres(i,j), twoe_couplings_cippres(i,j), j=ifanosta,ifanosta+10)
!      enddo
!      close(10)

      call ezfio_set_cippres_cfano_cippres(twoe_couplings_cippres)
      call ezfio_set_cippres_efano_cippres(e_couplings_cippres)
!      call ezfio_set_cippres_ifcsf(3)
     endif

  else 

    print*, "ifcsf = ", ifcsf
    print*, "but it should be equal to 2 for Fano calculations"
    print*, "Please run cippres_runci first or type qp set cippres ifcsf 2"
    print*, "Note that if you rerun cippres_fano, you should first delete the cfano/efano from EZFIO"

  endif

end program cippres_fano
