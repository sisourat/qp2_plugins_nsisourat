program phis
  use map_module
  implicit none
  BEGIN_DOC
! plugin to obtain input data for Fano-CI calculations
! scf and four_idx_transformation must be run before
  END_DOC
!  integer :: nmo
  integer ::  i, j, k, l, ik, jl
  integer(key_kind) :: i1
  double precision ::  mo_two_e_integral, temp

! call to mo_two_e_integral function does not function if one call to big_array_coulomb_integrals is not done before (load ints in memory?)
  temp = big_array_coulomb_integrals(nmo,nmo,nmo)

! e1=(i,k) and e2=(j,l)
!        triangular canonical order i>=k, j>=l,  ik>=jl

  nmo = mo_num
  open(unit=10,file="mo2eint.txt")
  write(10,*)'2e integrals',nmo
  do i = 1, nmo
    do j = 1, nmo
     do k = 1, i !, nmo
       do l = 1, j !, nmo
          ik = i*(i-1)/2+k
          jl = j*(j-1)/2+l
          if(ik>=jl) write(10,*)l,k,j,i,mo_two_e_integral(l,k,j,i)
       enddo
     enddo
   enddo
  enddo
  close(10)

  open(unit=10,file="dipmoint_x.txt")
  write(10,*)'dipMO x integrals',nmo
  do i = 1, nmo
    do j = i, nmo
      write(10,*)j,i,mo_dipole_x(j,i), mo_spread_x(j,i)
    enddo
  enddo
  close(10)

  open(unit=10,file="dipmoint_y.txt")
  write(10,*)'dipMO y integrals',nmo
  do i = 1, nmo
    do j = i, nmo
      write(10,*)j,i,mo_dipole_y(j,i), mo_spread_y(j,i)
    enddo
  enddo
  close(10)

  open(unit=10,file="dipmoint_z.txt")
  write(10,*)'dipMO z integrals',nmo
  do i = 1, nmo
    do j = i, nmo
      write(10,*)j,i,mo_dipole_z(j,i), mo_spread_z(j,i)
    enddo
  enddo
  close(10)

  open(unit=10,file="h1e_moint.txt")
  write(10,*)'h1e MO integrals',nmo
  do i = 1, nmo
    do j = i, nmo
      write(10,*)j,i,mo_one_e_integrals(j,i)
    enddo
  enddo
  close(10)


end
