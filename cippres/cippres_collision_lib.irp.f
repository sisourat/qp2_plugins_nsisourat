use bitmasks ! you need to include the bitmasks_module.f90 features
use general

 BEGIN_PROVIDER [integer, n_sta_coll_max]
  implicit none
     if(n_csf_max<5000) then
        n_sta_coll_max = n_csf_max
     else
        n_sta_coll_max = 5000
     endif
 END_PROVIDER 

 BEGIN_PROVIDER [double precision, b_coll]
  implicit none
   b_coll = 0d0
 END_PROVIDER 

 BEGIN_PROVIDER [double precision, v_coll]
  implicit none
     call ezfio_get_cippres_v_coll(v_coll)
 END_PROVIDER 

 BEGIN_PROVIDER [integer, i_state_coll]
  implicit none
     call ezfio_get_cippres_i_state_coll(i_state_coll)
 END_PROVIDER 

 BEGIN_PROVIDER [integer, stamin_bound]
  implicit none
     call ezfio_get_cippres_stamin_bound(stamin_bound)
 END_PROVIDER 

 BEGIN_PROVIDER [integer, stamax_bound]
  implicit none
     call ezfio_get_cippres_stamax_bound(stamax_bound)
 END_PROVIDER 

 BEGIN_PROVIDER [integer, stamin_si]
  implicit none
     call ezfio_get_cippres_stamin_si(stamin_si)
 END_PROVIDER 

 BEGIN_PROVIDER [integer, stamax_si]
  implicit none
     call ezfio_get_cippres_stamax_si(stamax_si)
 END_PROVIDER 

 BEGIN_PROVIDER [integer, stamin_di]
  implicit none
     call ezfio_get_cippres_stamin_di(stamin_di)
 END_PROVIDER 

 BEGIN_PROVIDER [integer, stamax_di]
  implicit none
     call ezfio_get_cippres_stamax_di(stamax_di)
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

 BEGIN_PROVIDER [double precision, coll_couplings_mo, (mo_num,mo_num,n_time)]
 use general
 implicit none
 integer :: i, j, k, l
 integer :: ib, ic, it
 double precision, dimension(mo_num,mo_num) :: w1e

 double precision :: t1, t2

 logical :: exists

 PROVIDE ezfio_filename !HF_bitmask mo_coef

   print*,'Computing coll_couplings', b_coll
   call cpu_time(t1)
   coll_couplings_mo(:,:,:) = 0d0

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

  do it = 1, n_time

    w1e(:,:) = 0d0
    do ic = 1, n_coulomb_center
      coulomb_center(1,ic) = b_coll
      coulomb_center(2,ic) = 0d0
      coulomb_center(3,ic) = zgrid(it)
      touch coulomb_center
      w1e(:,:) += charge_coulomb_center(ic)*mo_integrals_coulomb_center(:,:,ic)
    enddo
    coll_couplings_mo(:,:,it) = w1e(:,:)
   enddo

   call cpu_time(t2)
   print*,t2-t1
   print*,' '

 END_PROVIDER


 BEGIN_PROVIDER [double precision, coll_couplings, (n_sta_coll_max,n_sta_coll_max,n_time)]
 use general
 implicit none
 integer :: i, j, k, l
 integer :: ib, ic, it
 double precision, allocatable :: eigval1(:),eigvec1(:,:),eigval2(:),eigvec2(:,:),coll_csf_mat(:,:),coll_mat(:,:), mattmp(:,:)
 double precision, dimension(mo_num,mo_num) :: w1e

 integer :: nsta, ncsf
 double precision :: hij

 double precision :: t1, t2

 integer :: ni, nf

 logical :: exists

 PROVIDE ezfio_filename !HF_bitmask mo_coef

   print*,'Computing coll_couplings', b_coll
   call cpu_time(t1)
   coll_couplings(:,:,:) = 0d0

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

   if(stamax_di == 0) then
    nsta = stamax_bound-stamin_bound+1
    ni = stamin_bound
    nf = stamax_bound
   else
    nsta = stamax_di-stamin_bound+1
    ni = stamin_bound
    nf = stamax_di
   endif

   if(nsta>n_sta_coll_max) then
    print*, "nsta > n_sta_coll_max, I stop"
    stop
   endif

   ncsf = n_csf_cippres(ici1)
   allocate(coll_mat(nsta,nsta))
   coll_mat(:,:) = 0d0
   allocate(mattmp(ncsf,nsta))

  do it = 1, n_time

    w1e(:,:) = 0d0
    do ic = 1, n_coulomb_center
      coulomb_center(1,ic) = b_coll
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
    CALL DGEMM('N','N',ncsf,nsta,ncsf,1.d0,coll_csf_mat(1:ncsf,1:ncsf),ncsf,eigvec1(1:ncsf,ni:nf),ncsf,0.d0,mattmp,ncsf)
    CALL DGEMM('N','N',nsta,nsta,ncsf,1.d0,transpose(eigvec2(1:ncsf,ni:nf)),nsta,mattmp,ncsf,0.d0,coll_mat(1:nsta,1:nsta),nsta)

    coll_couplings(ni:nf,ni:nf,it) = coll_mat(1:nsta,1:nsta)
   enddo

   call cpu_time(t2)
   print*,t2-t1
   print*,' '

 deallocate(coll_csf_mat,eigval1,eigval2,eigvec1,eigvec2,coll_mat,mattmp) 

 END_PROVIDER

