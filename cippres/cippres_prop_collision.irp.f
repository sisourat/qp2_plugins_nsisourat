program cippres_prop_collision
  use general
  use propdyn
 ! create a routine for one-e general matrix with determinants

  implicit none
  BEGIN_DOC
! CIPPRES stands for Configuration Interaction Plugin for Photoionized and Resonant Electronic States.
! cippres_prop_collision solves the TDSE using the matrix elements from cippres_collision
  END_DOC

  double precision :: t1, t2, tdyn
  logical :: exists

  integer :: i, j, ib, k, l, it, ni, nf

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
   print*, "" 
   print*, "Collision info"
   print*, "" 
   print*,"n_bimp,n_time,i_state_coll,v_coll =", n_bimp, n_time, i_state_coll, v_coll
   print*,"stamin, stamax, nsta =", stamin_bound, stamax_bound, stamax_bound-stamin_bound+1
   print*,"impact. param. =", bgrid
   print*,"tgrid = ", tgrid
   print*,"zgrid = ", zgrid
   print*, "" 

   call ezfio_get_cippres_ici1(ici1)
   call ezfio_get_cippres_n_csf_cippres(n_csf_cippres)

   ntime= n_time

   nsta = stamax_bound-stamin_bound+1
   nsi  = stamax_si-stamin_si+1
   ndi  = stamax_di-stamin_di+1
   ni = stamin_bound
   nf = stamax_di
   if(stamax_di == 0) then  
      nsi = 0
      ndi = 0 
      ni = stamin_bound
      nf = stamax_bound
   endif
   ntotsta = nsta + nsi + ndi

   allocate(mcoup(ntime,nsta,nsta),timegrid(ntime),esta(ntotsta),mat(ntime,ntotsta,ntotsta,n_bimp))
   allocate(psi(nsta))
   allocate(rmat2intrp(ntime,nsta,nsta),cmat2intrp(ntime,nsta,nsta))
   allocate(matintrp(nsta,nsta))

   allocate(g_si(1:ntime,1:nsta,1:nsta),g_ddi(1:ntime,1:nsta,1:nsta),g_sdi(1:ntime,1:nsta))
   allocate(g_si_intrp(1:ntime,1:nsta,1:nsta),g_ddi_intrp(1:ntime,1:nsta,1:nsta),g_sdi_intrp(1:ntime,1:nsta))

   timegrid(1:ntime) = tgrid(1:n_time)
   esta(1:ntotsta) = eigvalues_cippres(ni:nf,ici1)

! loop over b

   open(unit=20,file='Prop_collision.out')
   write(20,*)n_bimp,nsta
   do i = 1, nsta
     write(20,*)esta(i)
     write(*,*)esta(i)
   enddo

   do it = 1, ntime
     mat(it,1:ntotsta,1:ntotsta,1:n_bimp) = coll_couplings(ni:nf,ni:nf,it,1:n_bimp)
   enddo

   print*,'start dyn'
   do ib = 1, n_bimp

    g_si(:,:,:) = 0d0; p_si = 0d0
    g_ddi(:,:,:) = 0d0; p_sdi = 0d0
    g_sdi(:,:) = 0d0; p_ddi = 0d0

    mcoup(:,:,:) = 0d0
    mcoup(:,1:nsta,1:nsta) = mat(:,1:nsta,1:nsta,ib)

    do i = 1, nsta
      do j = 1, nsta

! coupling between bound states and single ionization states
!        do k = nsta+1, nsta+nsi
!          do l = nsta+1, nsta+nsi
!            g_si(:,j,i) = g_si(:,j,i) + mat(:,l,j,ib)*mat(:,k,i,ib)
!          enddo
!        enddo

! coupling between bound states and double ionization states
!         do k = nsta+nsi+1, nsta+nsi+ndi
!           do l = nsta+nsi+1, nsta+nsi+ndi
!             g_ddi(:,j,i) = g_ddi(:,j,i) + mat(:,l,j,ib)*mat(:,k,i,ib)
!           enddo
!         enddo

      enddo ! loop over j

! coupling between bound states and single ionization states
        do k = nsta+1, nsta+nsi
            g_si(:,i,i) = g_si(:,i,i) + mat(:,k,i,ib)**2/abs(esta(i)-esta(k))
        enddo
!   do it = 1, ntime
!     print*,timegrid(it),g_si(it,1,1)
!   enddo

       
! coupling between bound states and double ionization states via single ionization states
      do k = nsta+nsi+1, nsta+nsi+ndi
        do l = nsta+1, nsta+nsi
          g_sdi(:,i) = g_sdi(:,i) + mat(:,l,i,ib)*mat(:,l,k,ib)  
        enddo
      enddo

    enddo ! loop over i

!   do it = 1, ntime
!    write(*,*)g_si(it,1,1)
!   enddo

    psi(:) = 0d0
    psi(i_state_coll) = 1d0

! propagation 
    call cpu_time(t1)
    call dyn
    call cpu_time(t2)
    tdyn = t2-t1
    write(*,*)'DYN takes',tdyn

!     print*,bgrid(ib),nsta,psi
    write(*,'(5000(f12.6,1X))')bgrid(ib),(cdabs(psi(i))**2,i=1,nsta), sum(cdabs(psi(:))**2),v_coll
    write(20,'(5000(f12.6,1X))')bgrid(ib),(cdabs(psi(i))**2,i=1,nsta), sum(cdabs(psi(:))**2),v_coll

   enddo 

   close(20)

   deallocate(mcoup,timegrid,esta,mat)
   deallocate(matintrp,psi,rmat2intrp,cmat2intrp)
   deallocate(g_si,g_sdi,g_ddi)
   deallocate(g_si_intrp,g_sdi_intrp,g_ddi_intrp)

  else

    print*, "Z and/or b grids are not setup correctly"
    stop

  endif

end program cippres_prop_collision


