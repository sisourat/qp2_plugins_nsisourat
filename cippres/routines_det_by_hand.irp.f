subroutine create_det(n_a,n_b,occ_a,occ_b,det)
 use bitmasks ! you need to include the bitmasks_module.f90 features
 integer, intent(in) :: n_a,n_b,occ_a(mo_num),occ_b(mo_num)
 integer(bit_kind), intent(out) :: det(N_int,2)
 integer :: iorb,i
 det = 0_bit_kind
 do i = 1, n_a
  iorb = occ_a(i)
  call set_bit_to_integer(iorb,det(1,1),N_int)
 enddo
 do i = 1, n_b
  iorb = occ_b(i)
  call set_bit_to_integer(iorb,det(1,2),N_int)
 enddo
end

subroutine i_w1e_j(key_i,key_j,Nint,w1e,hij)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns $\langle i| w1e |j \rangle$  where $i$ and $j$ are determinants and w1e the matrix elements of a one-e operator in the MO basis.
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: key_i(Nint,2), key_j(Nint,2)
  double precision, intent(in) :: w1e(mo_num,mo_num)
  double precision, intent(out)  :: hij

  integer :: degree,m,p
  double precision :: diag_w1e_mat_elem_one_e,phase
  integer                        :: exc(0:2,2,2)
  call get_excitation_degree(key_i,key_j,degree,Nint)
  hij = 0.d0
  if(degree>1)then
   return
  endif
  if(degree==0)then
   hij = diag_w1e_mat_elem_one_e(key_i,N_int,w1e)
  else
   call get_single_excitation(key_i,key_j,exc,phase,Nint)
   if (exc(0,1,1) == 1) then
     ! Mono alpha
     m = exc(1,1,1)
     p = exc(1,2,1)
   else
     ! Mono beta
     m = exc(1,1,2)
     p = exc(1,2,2)
   endif
   hij = phase * w1e(m,p)
  endif

end

double precision function diag_w1e_mat_elem_one_e(det_in,Nint,w1e)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Computes $\langle i|w1e|i \rangle$.
  END_DOC
  double precision, intent(in) :: w1e(mo_num,mo_num)
  integer,intent(in)             :: Nint
  integer(bit_kind),intent(in)   :: det_in(Nint,2)

  integer(bit_kind)              :: hole(Nint,2)
  integer(bit_kind)              :: particle(Nint,2)
  integer                        :: i, nexc(2), ispin
  integer                        :: occ_particle(Nint*bit_kind_size,2)
  integer                        :: occ_hole(Nint*bit_kind_size,2)
  integer(bit_kind)              :: det_tmp(Nint,2)
  integer                        :: na, nb

  ASSERT (Nint > 0)
  ASSERT (sum(popcnt(det_in(:,1))) == elec_alpha_num)
  ASSERT (sum(popcnt(det_in(:,2))) == elec_beta_num)

  diag_w1e_mat_elem_one_e = 0.d0

  !call debug_det(det_in,Nint)
  integer                        :: tmp(2)
  !DIR$ FORCEINLINE
  call bitstring_to_list_ab(det_in, occ_particle, tmp, Nint)
  do ispin = 1,2
   do i = 1, tmp(ispin)
!    print*,occ_particle(i,ispin),w1e(occ_particle(i,ispin),occ_particle(i,ispin))
    diag_w1e_mat_elem_one_e +=  w1e(occ_particle(i,ispin),occ_particle(i,ispin))
   enddo
  enddo

end

