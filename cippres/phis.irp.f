program phis
  use map_module
  implicit none
  BEGIN_DOC
! plugin to obtain input data for Fano-CI calculations
! scf and four_idx_transformation must be run before
  END_DOC
  integer :: nmo
  integer ::  i, j, k, l
  integer(key_kind) :: i1
  double precision ::  mo_two_e_integral, temp

! call to mo_two_e_integral function does not function if one call to big_array_coulomb_integrals is not done before (load ints in memory?)
  temp = big_array_coulomb_integrals(nmo,nmo,nmo)

  nmo = mo_num
  write(*,*)'2e integrals',nmo
  do i = 1, nmo
    do j = 1, nmo
     do k = 1, nmo
       do l = 1, nmo
!          call mo_two_e_integrals_index(l,k,j,i,i1)
          write(*,*)l,k,j,i,mo_two_e_integral(l,k,j,i)
          !write(*,*)l,k,j,i,mo_two_e_integral(11,11,11,11)
       enddo
     enddo
   enddo
  enddo

  write(*,*)''
  write(*,*)'dipMO x integrals',nmo
  do i = 1, nmo
    do j = 1, nmo
      write(*,*)j,i,mo_dipole_x(j,i)
    enddo
  enddo

  write(*,*)''
  write(*,*)'dipMO y integrals',nmo
  do i = 1, nmo
    do j = 1, nmo
      write(*,*)j,i,mo_dipole_y(j,i)
    enddo
  enddo

  write(*,*)''
  write(*,*)'dipMO z integrals',nmo
  do i = 1, nmo
    do j = 1, nmo
      write(*,*)j,i,mo_dipole_z(j,i)
    enddo
  enddo

end
