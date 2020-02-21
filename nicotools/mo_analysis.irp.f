program mo_analysis
 implicit none
 integer :: nx,i
 double precision :: xmax,dx,r(3),weight,accu
 double precision, allocatable :: dm_a(:),dm_b(:),mos_array(:)
 allocate(dm_a(N_states),dm_b(N_states),mos_array(mo_num))

!xmax = 5.d0 
!nx = 1000
!dx = xmax/dble(nx)
!r = 0.d0
!do i = 1, nx
! call dm_dft_alpha_beta_at_r(r,dm_a,dm_b)  ! density 
! call give_all_mos_at_r(r,mos_array) ! mos 
! r(1) += dx
! write(33,'(100(F16.10,X))')r(1),dm_a(1:10)+dm_b(1:10)
!enddo

 call myprint_mulliken_per_mo

end

