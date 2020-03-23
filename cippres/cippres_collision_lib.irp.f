use bitmasks ! you need to include the bitmasks_module.f90 features
use general

 BEGIN_PROVIDER [double precision, v_coll]
  implicit none
     call ezfio_get_cippres_v_coll(v_coll)
 END_PROVIDER 

 BEGIN_PROVIDER [integer, i_state_coll]
  implicit none
     call ezfio_get_cippres_i_state_coll(i_state_coll)
 END_PROVIDER 

 BEGIN_PROVIDER [integer, stamin_coll]
  implicit none
     call ezfio_get_cippres_stamin_coll(stamin_coll)
 END_PROVIDER 

 BEGIN_PROVIDER [integer, stamax_coll]
  implicit none
     call ezfio_get_cippres_stamax_coll(stamax_coll)
 END_PROVIDER 

 BEGIN_PROVIDER [integer, n_time]
  implicit none
     call ezfio_get_cippres_n_time(n_time)
 END_PROVIDER 

 BEGIN_PROVIDER [integer, n_bimp]
  implicit none
     call ezfio_get_cippres_n_bimp(n_bimp)
 END_PROVIDER 

 BEGIN_PROVIDER [integer, n_pcenter]
  implicit none
     call ezfio_get_cippres_n_pcenter(n_pcenter)
 END_PROVIDER 

 BEGIN_PROVIDER [double precision, charge_pcenter, (n_pcenter)]
  implicit none
     call ezfio_get_cippres_charge_pcenter(charge_pcenter)
 END_PROVIDER 

 BEGIN_PROVIDER [double precision, bgrid, (n_bimp)]
  implicit none
     call ezfio_get_cippres_bgrid(bgrid)
 END_PROVIDER 

 BEGIN_PROVIDER [double precision, zgrid, (n_time)]
  implicit none
     call ezfio_get_cippres_zgrid(zgrid)
 END_PROVIDER 

 BEGIN_PROVIDER [double precision, tgrid, (n_time)]
  implicit none
     call ezfio_get_cippres_tgrid(tgrid)
 END_PROVIDER 

 BEGIN_PROVIDER [double precision, coll_couplings_cippres, (n_csf_max,n_csf_max,n_time,n_bimp)]
 use general
 implicit none
 integer :: i, j, k, l
 integer :: ib, ic, it
 double precision, allocatable :: eigval1(:),eigvec1(:,:),eigval2(:),eigvec2(:,:),coll_csf_mat(:,:),coll_mat(:,:), mattmp(:,:)
 double precision, dimension(mo_num,mo_num) :: w1e

 integer :: nsta, ncsf
 double precision :: hij

 logical :: exists

 PROVIDE ezfio_filename !HF_bitmask mo_coef
 call ezfio_has_cippres_coll_couplings_cippres(exists)

 if (exists) then

   call ezfio_get_cippres_coll_couplings_cippres(coll_couplings_cippres)

 else

   coll_couplings_cippres(:,:,:,:) = 0d0

   allocate(coll_csf_mat(n_csf_cippres(ici1),n_csf_cippres(ici1)))
 
   allocate(eigval1(n_csf_cippres(ici1)),eigval2(n_csf_cippres(ici1)))
   eigval1(:) = eigvalues_cippres(1:n_csf_cippres(ici1),ici1)
   eigval2(:) = eigvalues_cippres(1:n_csf_cippres(ici1),ici1)

   allocate(eigvec1(n_csf_cippres(ici1),n_csf_cippres(ici1)))
   allocate(eigvec2(n_csf_cippres(ici1),n_csf_cippres(ici1)))
   eigvec1(:,:) = eigvectors_cippres(1:n_csf_cippres(ici1),1:n_csf_cippres(ici1),ici1)
   eigvec2(:,:) = eigvectors_cippres(1:n_csf_cippres(ici1),1:n_csf_cippres(ici1),ici1)


   if (mpi_master) then
     call ezfio_has_cippres_n_pcenter(exists)
      if (exists) then
        call ezfio_has_cippres_charge_pcenter(exists)
      endif
   endif

   call ezfio_get_cippres_n_pcenter(n_pcenter)
   call ezfio_get_cippres_charge_pcenter(charge_pcenter)
   n_coulomb_center = n_pcenter
   touch n_coulomb_center
   charge_coulomb_center(1:n_pcenter) = charge_pcenter(1:n_pcenter)
   touch charge_coulomb_center

   call ezfio_get_cippres_stamin_coll(stamin_coll)
   call ezfio_get_cippres_stamax_coll(stamax_coll)

   nsta = stamax_coll-stamin_coll+1
   ncsf = n_csf_cippres(ici1)
   allocate(coll_mat(nsta,nsta))
   coll_mat(:,:) = 0d0
   allocate(mattmp(ncsf,nsta))

 do ib = 1, n_bimp
  do it = 1, n_time

    w1e(:,:) = 0d0
    do ic = 1, n_coulomb_center
      coulomb_center(1,ic) = bgrid(ib)
      coulomb_center(2,ic) = 0d0 
      coulomb_center(3,ic) = zgrid(it)
      touch coulomb_center
      w1e(:,:) += charge_coulomb_center(ic)*mo_integrals_coulomb_center(:,:,ic)
    enddo

    coll_csf_mat(:,:) = 0d0

!$OMP PARALLEL DO PRIVATE(i,j,k,l,hij)
! SCHEDULE(DYNAMIC) 
    do i = 1, n_csf_cippres(ici1) ! first loop on the csf of the space ispace 
     do j = i, n_csf_cippres(ici1)
      do k = 1, n_det_csf_cippres(i,ici1) ! then on the determinants belonging to the ith CSF of space ispace
       do l = 1, n_det_csf_cippres(j,ici1)
          call i_w1e_j(csf_basis(1,1,k,i,ici1),csf_basis(1,1,l,j,ici1),N_int,w1e,hij)
          coll_csf_mat(j,i) += hij * coef_det_csf_basis(k,i,ici1) * coef_det_csf_basis(l,j,ici1)
       enddo
      enddo
      coll_csf_mat(i,j) = coll_csf_mat(j,i)
     enddo
    enddo
!$OMP END PARALLEL DO

    coll_mat(:,:) = 0d0
!    do i = 1, n_csf_cippres(ici1) ! first loop on the first eigenvectors
!     do j = 1, n_csf_cippres(ici1) ! then on the second eigenvectors

!!$OMP PARALLEL DO PRIVATE(i,j,k,l)
!! SCHEDULE(DYNAMIC) 
!    do i = stamin_coll, stamax_coll
!     do j = stamin_coll, stamax_coll
!      do k = 1, n_csf_cippres(ici1) ! loop over the csfs of the ici1 run
!       do l = 1, n_csf_cippres(ici1) ! then over the csfs of the ici1 run
!          coll_mat(j,i) += coll_csf_mat(l,k) * eigvec1(k,i) * eigvec2(l,j)
!       enddo
!      enddo
!     enddo
!    enddo
!!$OMP END PARALLEL DO

    CALL DGEMM('N','N',ncsf,nsta,ncsf,1.d0,coll_csf_mat(1:ncsf,1:ncsf),ncsf,eigvec1(1:ncsf,stamin_coll:stamax_coll),ncsf,0.d0,mattmp,ncsf)
    CALL DGEMM('N','N',nsta,nsta,ncsf,1.d0,transpose(eigvec2(1:ncsf,stamin_coll:stamax_coll)),nsta,mattmp,ncsf,0.d0,coll_mat(1:nsta,1:nsta),nsta)

    coll_couplings_cippres(stamin_coll:stamax_coll,stamin_coll:stamax_coll,it,ib) = coll_mat(1:nsta,1:nsta)
   enddo
 enddo

 deallocate(coll_csf_mat,eigval1,eigval2,eigvec1,eigvec2,coll_mat,mattmp) 

 endif

 END_PROVIDER

