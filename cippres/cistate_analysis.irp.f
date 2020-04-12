program cistate_analysis
 implicit none
 integer :: irun,nx,i,istate
 irun = 1
 print*,'irun=',irun
 call from_csf_to_det(irun)
 double precision :: xmax,dx,r(3),weight,accu
 double precision, allocatable :: dm_a(:),dm_b(:),mos_array(:)
 allocate(dm_a(N_states),dm_b(N_states),mos_array(mo_num))
!xmax = 5.d0 
!nx = 1000
!dx = xmax/dble(nx)
!r = 0.d0
!do i = 1, nx
! call dm_dft_alpha_beta_at_r(r,dm_a,dm_b)  ! density 
! r(1) += dx
! write(33,'(100(F16.10,X))')r(1),dm_a(1:10)+dm_b(1:10)
!enddo
!
!accu = 0.d0
!do i = 1, n_points_final_grid
! r(1) = final_grid_points(1,i)
! r(2) = final_grid_points(2,i)
! r(3) = final_grid_points(3,i)
! weight = final_weight_at_r_vector(i)
! call dm_dft_alpha_beta_at_r(r,dm_a,dm_b)  ! density 
! accu += (dm_a(1) + dm_b(1)) * weight
!enddo
!print*,'accu = ',accu
!

 call myprint_mulliken
end

