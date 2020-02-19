subroutine read_csf_basis(csf)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Reads the determinants from the |EZFIO| file
  END_DOC

  integer(bit_kind), intent(out) :: csf(N_int,2,n_det_max_csf, n_csf_max ,n_ciruns_cippres)
  integer                        :: N_int2
  integer                        :: i,k

  call ezfio_get_determinants_N_int(N_int2)
  ASSERT (N_int2 == N_int)
  call ezfio_get_determinants_bit_kind(k)
  ASSERT (k == bit_kind)
  if(bit_kind .ne. 8)then
   print*,'The bit_kind used for the csf is ',bit_kind
   print*,'But it should be 8'
   stop
  endif
  call ezfio_get_cippres_csf_basis(csf)

end

