program cippres_prop_collision_mo
  use general
  use propdyn
  use map_module
 ! create a routine for one-e general matrix with determinants

  implicit none
  BEGIN_DOC
! CIPPRES stands for Configuration Interaction Plugin for Photoionized and Resonant Electronic States.
! cippres_prop_collision_mo solves the TDSE using the matrix elements in the MO basis
  END_DOC

  double precision :: t1, t2, tdyn, p1, p2, p3
  double precision :: hfen, r12en, en
  logical :: exists

  integer :: i, j, ib, k, l, it
  integer :: nocc, nact, ioccmin, ioccmax, iactmin, iactmax
  integer :: isave_time, nsave_time

!!nico  integer :: nmo
  integer(key_kind) :: i1
  double precision  :: mo_two_e_integral, temp

  double precision, dimension(:), allocatable :: prob_save, estatmp
  complex (kind=8), dimension(:,:), allocatable :: psit_save

! call to mo_two_e_integral function does not function if one call to big_array_coulomb_integrals is not done before (load ints in memory?)
  nmo = mo_num
  temp = big_array_coulomb_integrals(nmo,nmo,nmo)

  PROVIDE ezfio_filename !HF_bitmask mo_coef

  if (mpi_master) then
   call ezfio_has_cippres_n_time(exists)
    if (exists) then
      call ezfio_has_cippres_n_bimp(exists)
      if (exists) then
       call ezfio_has_cippres_tgrid(exists)
       if (exists) then
         call ezfio_has_cippres_zgrid(exists)
         if (exists) then
           call ezfio_has_cippres_bgrid(exists)
            if (exists) then
             call ezfio_has_cippres_v_coll(exists)
             if (exists) then
               call ezfio_has_cippres_i_state_coll(exists)
               if (exists) then
                 call ezfio_has_cippres_stamin_bound(exists)
                 if (exists) then
                   call ezfio_has_cippres_stamax_bound(exists)
                 endif
               endif
             endif
           endif
         endif
       endif
     endif
   endif
 endif

  if (exists) then
   call ezfio_get_cippres_n_time(n_time)
   call ezfio_get_cippres_n_bimp(n_bimp)
   call ezfio_get_cippres_zgrid(zgrid)
   call ezfio_get_cippres_tgrid(tgrid)
   call ezfio_get_cippres_bgrid(bgrid)
   call ezfio_get_cippres_stamin_bound(stamin_bound)
   call ezfio_get_cippres_stamax_bound(stamax_bound)
   call ezfio_get_cippres_stamin_si(stamin_si)
   call ezfio_get_cippres_stamax_si(stamax_si)
   call ezfio_get_cippres_stamin_di(stamin_di)
   call ezfio_get_cippres_stamax_di(stamax_di)
   call ezfio_get_cippres_i_state_coll(i_state_coll)
   call ezfio_get_cippres_v_coll(v_coll)

   ioccmin = stamin_bound
   ioccmax = stamax_bound
   nocc = ioccmax - ioccmin + 1
   iactmin = stamin_si
   iactmax = stamax_si
   nact = iactmax - iactmin + 1

   print*, "" 
   print*, "Collision info"
   print*, "" 
   print*,"n_bimp,n_time,i_state_coll,v_coll =", n_bimp, n_time, i_state_coll, v_coll
   print*,"actmin, actmax, nact =",iactmin,iactmax,nact 
   print*,"impact. param. =", bgrid
   print*,"tgrid = ", tgrid
   print*,"zgrid = ", zgrid
   print*, "" 

   ntime= n_time
   nsave_time= n_time

   nsta = nact + 1 ! +1 for the active occupied MO electron
   allocate(mcoup(ntime,nsta,nsta),timegrid(ntime),esta(nsta),estatmp(nsta),mat(ntime,nsta,nsta))
   allocate(psi(nsta),psit_save(nsta,nsave_time),prob_save(nsta))
   allocate(rmat2intrp(ntime,nsta,nsta),cmat2intrp(ntime,nsta,nsta))
   allocate(matintrp(nsta,nsta))

   timegrid(1:ntime) = tgrid(1:n_time)
   
!   hfen = 0d0
!   r12en = 0d0
!   do i = 1, nocc
!     j = ioccmin + i - 1
!     hfen = hfen + mo_one_e_integrals(j,j)
!     do k = 1, nocc
!        l = ioccmin + k - 1
!       r12en = r12en + 2d0*mo_two_e_integral(j,l,j,l)-mo_two_e_integral(j,j,l,l)
!     enddo
!   enddo

   j = 1 
   i = i_state_coll
   esta(j) =  mo_one_e_integrals(i,i) 
    do k = 1, nocc
      l = ioccmin + k - 1
      esta(j) = esta(j) + 2d0*mo_two_e_integral(i,l,i,l)-mo_two_e_integral(i,i,l,l)
    enddo

   do i = iactmin, iactmax
    j = j + 1
    esta(j) =  mo_one_e_integrals(i,i) 
    do k = 1, nocc
        l = ioccmin + k - 1
        esta(j) = esta(j) + mo_two_e_integral(i,l,i,l)+mo_two_e_integral(i,i,l,l)
    enddo
   enddo
   write(*,*)esta
   estatmp(:) = esta(:)

! loop over b

   open(unit=20,file='Prop_collision_mo.out')
   open(unit=21,file='Prop_collision_mo_testsi.out')
   open(unit=22,file='Prop_collision_mo_testdi.out')
   open(unit=30,file='Psit_collision_mo.out')
   write(20,*)n_bimp,nsta
   write(21,*)n_bimp,nsta
   write(22,*)n_bimp,nsta
   write(30,*)n_bimp,nsta,nsave_time
   do i = 1, nsta
     write(20,*)esta(i)-nuclear_repulsion
     write(21,*)esta(i)-nuclear_repulsion
     write(22,*)esta(i)-nuclear_repulsion
     write(30,*)esta(i)-nuclear_repulsion
     write(*,*)esta(i)
   enddo


   print*,'start dyn'
   do ib = 1, n_bimp
    esta(:) = estatmp(:)

    b_coll = bgrid(ib)
    touch b_coll
    do it = 1, ntime
      mat(it,1,1) = coll_couplings_mo(i_state_coll,i_state_coll,it)
      mat(it,1,2:nsta) = dsqrt(2d0)*coll_couplings_mo(i_state_coll,iactmin:iactmax,it)
      mat(it,2:nsta,1) = dsqrt(2d0)*coll_couplings_mo(iactmin:iactmax,i_state_coll,it)
      do i = 2, nsta
        do j = 2, nsta
          k = iactmin + i - 2
          l = iactmin + j - 2
          mat(it,j,i) = coll_couplings_mo(l,k,it)
        enddo
      enddo
      mat(it,nsta,nsta) = coll_couplings_mo(iactmax,iactmax,it)
      mat(it,nsta,2:nsta) = coll_couplings_mo(iactmax,iactmin:iactmax,it)
      mat(it,2:nsta,nsta) = coll_couplings_mo(iactmin:iactmax,iactmax,it)
    enddo

    mcoup(:,:,:) = 0d0
    mcoup(:,1:nsta,1:nsta) = mat(:,1:nsta,1:nsta)
!   do it = 1, ntime
!     print*,timegrid(it),coll_couplings_mo(1,1,it),coll_couplings_mo(1,2,it)
!     print*,timegrid(it),cdabs(mcoup(it,1,1)),cdabs(mcoup(it,1,2))
!   enddo
!   stop

    psi(:) = 0d0
    psi(1) = 1d0

! propagation 
    call cpu_time(t1)
    call dyn
    call cpu_time(t2)
    tdyn = t2-t1
    write(*,*)'DYN takes',tdyn

!    print*,bgrid(ib),nsta,psi
    prob_save(:) = cdabs(psi(:))**2
    write(*,'(5000(f12.6,1X))')bgrid(ib),(cdabs(psi(i))**2,i=1,nsta), sum(cdabs(psi(:))**2),v_coll
    write(20,'(5000(f12.6,1X))')bgrid(ib),(cdabs(psi(i))**2,i=1,nsta), sum(cdabs(psi(:))**2),v_coll

    do isave_time = 1, nsave_time
      write(30,'(5000(f22.16,1X))')bgrid(ib),zgrid(isave_time),psit_save(:,isave_time)
    enddo

   enddo 

   close(20)
   close(21)
   close(22)

   deallocate(mcoup,timegrid,esta,estatmp,mat)
   deallocate(matintrp,psi,rmat2intrp,cmat2intrp)

  else

    print*, "Z and/or b grids are not setup correctly"
    stop

  endif

end program cippres_prop_collision_mo


