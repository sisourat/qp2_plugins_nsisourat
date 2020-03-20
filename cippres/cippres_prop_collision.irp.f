program cippres_prop_collision
  use general
  use propdyn
 ! create a routine for one-e general matrix with determinants

  implicit none
  BEGIN_DOC
! CIPPRES stands for Configuration Interaction Plugin for Photoionized and Resonant Electronic States.
! cippres_prop_collision solves the TDSE using the matrix elements from cippres_collision
  END_DOC

  logical :: exists

  integer :: i, j, ib

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
              call ezfio_has_cippres_coll_couplings_cippres(exists)
              if (exists) then
                call ezfio_has_cippres_v_coll(exists)
                if (exists) then
                  call ezfio_has_cippres_i_state_coll(exists)
                  if (exists) then
                    call ezfio_has_cippres_stamin_coll(exists)
                    if (exists) then
                      call ezfio_has_cippres_stamax_coll(exists)
                    endif
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
   call ezfio_get_cippres_coll_couplings_cippres(coll_couplings_cippres)
   call ezfio_get_cippres_stamin_coll(stamin_coll)
   call ezfio_get_cippres_stamax_coll(stamax_coll)
   call ezfio_get_cippres_i_state_coll(i_state_coll)
   call ezfio_get_cippres_v_coll(v_coll)
   print*, "" 
   print*, "Collision info"
   print*, "" 
   print*,"n_bimp,n_time,i_state_coll,v_coll =", n_bimp, n_time, i_state_coll, v_coll
   print*,"stamin, stamax, ntotsta =", stamin_coll, stamax_coll, stamax_coll-stamin_coll
    print*,"impact. param. =", bgrid
   print*,"tgrid = ", tgrid
   print*,"zgrid = ", zgrid
   print*, "" 

   call ezfio_get_cippres_ici1(ici1)
   call ezfio_get_cippres_n_csf_cippres(n_csf_cippres)

   ntime= n_time
   ntotsta = stamax_coll-stamin_coll

   allocate(mcoup(ntime,ntotsta,ntotsta),timegrid(ntime),esta(ntotsta),mat(ntotsta,ntotsta))
   allocate(psi(ntotsta))
   allocate(rmat2intrp(ntime,ntotsta,ntotsta),cmat2intrp(ntime,ntotsta,ntotsta))
   allocate(matintrp(ntotsta,ntotsta))

   timegrid(1:ntime) = tgrid(1:n_time)
   esta(1:ntotsta) = eigvalues_cippres(stamin_coll:stamax_coll,ici1)

! loop over b

   open(unit=20,file='Prop_collision.out')
   write(20,*)n_bimp,ntotsta
   do i = 1, ntotsta
     write(20,*)esta(i)
     write(*,*)esta(i)
   enddo

   print*,'start dyn'
   do ib = 1, n_bimp

    do i = 1, ntime
     mcoup(i,1:ntotsta,1:ntotsta) = dcmplx(coll_couplings_cippres(stamin_coll:stamax_coll,stamin_coll:stamax_coll,i,ib),0d0)
!     print*,'mcoup',i,mcoup(i,1,1)
    enddo

    psi(:) = 0d0
    psi(i_state_coll) = 1d0
! propagation 
!  call cpu_time(t1)
    call dyn
!  call cpu_time(t2)
!  tdyn = t2-t1
!  write(*,*)'DYN takes',tdyn
!     print*,bgrid(ib),ntotsta,psi
    write(*,'(5000(f12.6,1X))')bgrid(ib),(cdabs(psi(i))**2,i=1,ntotsta), sum(cdabs(psi(:))**2),v_coll
    write(20,'(5000(f12.6,1X))')bgrid(ib),(cdabs(psi(i))**2,i=1,ntotsta), sum(cdabs(psi(:))**2),v_coll
   enddo 

   close(20)

   deallocate(mcoup,timegrid,esta,mat)
   deallocate(matintrp,psi,rmat2intrp,cmat2intrp)

  else

    print*, "Z and/or b grids are not setup correctly"
    stop

  endif

end program cippres_prop_collision


