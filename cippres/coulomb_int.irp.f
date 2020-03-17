BEGIN_PROVIDER [integer, n_coulomb_center]
 implicit none
 n_coulomb_center = 1
END_PROVIDER 

BEGIN_PROVIDER [double precision, coulomb_center, (3,n_coulomb_center)]
 implicit none
 coulomb_center = 0.d0
END_PROVIDER 

BEGIN_PROVIDER [double precision, charge_coulomb_center, (n_coulomb_center)]
 implicit none
 charge_coulomb_center = 0.d0
END_PROVIDER

BEGIN_PROVIDER [ double precision, ao_integrals_coulomb_center, (ao_num,ao_num,n_coulomb_center)]
  BEGIN_DOC
! Nucleus-electron interaction in the |AO| basis set, per atom Coulomb center 
!
! :math:`\langle \chi_i | -\frac{1}{|r-R_A|} | \chi_j \rangle`
  END_DOC
  implicit none
  double precision               :: alpha, beta, gama, delta
  integer                        :: i_c,num_A,num_B
  double precision               :: A_center(3),B_center(3),C_center(3)
  integer                        :: power_A(3),power_B(3)
  integer                        :: i,j,k,l,n_pt_in,m
  double precision               :: overlap_x,overlap_y,overlap_z,overlap,dx,NAI_pol_mult

  ao_integrals_coulomb_center = 0.d0

  !$OMP PARALLEL                                                    &
      !$OMP DEFAULT (NONE)                                          &
      !$OMP PRIVATE (i,j,k,l,m,alpha,beta,A_center,B_center,power_A,power_B,&
      !$OMP  num_A,num_B,c,n_pt_in,C_center)                        &
      !$OMP SHARED (ao_num,ao_prim_num,ao_expo_ordered_transp,ao_power,ao_nucl,nucl_coord,ao_coef_normalized_ordered_transp,&
      !$OMP  n_pt_max_integrals,ao_integrals_coulomb_center,n_coulomb_center,coulomb_center)
  n_pt_in = n_pt_max_integrals
  !$OMP DO SCHEDULE (dynamic)

  double precision               :: c
  do j = 1, ao_num
    power_A(1)= ao_power(j,1)
    power_A(2)= ao_power(j,2)
    power_A(3)= ao_power(j,3)
    num_A = ao_nucl(j)
    A_center(1) = nucl_coord(num_A,1)
    A_center(2) = nucl_coord(num_A,2)
    A_center(3) = nucl_coord(num_A,3)
    do  k = 1, n_coulomb_center
      C_center(1) = coulomb_center(1,k)
      C_center(2) = coulomb_center(2,k)
      C_center(3) = coulomb_center(3,k)
      do i = 1, ao_num
        power_B(1)= ao_power(i,1)
        power_B(2)= ao_power(i,2)
        power_B(3)= ao_power(i,3)
        num_B = ao_nucl(i)
        B_center(1) = nucl_coord(num_B,1)
        B_center(2) = nucl_coord(num_B,2)
        B_center(3) = nucl_coord(num_B,3)
        c = 0.d0
        do l=1,ao_prim_num(j)
          alpha = ao_expo_ordered_transp(l,j)
          do m=1,ao_prim_num(i)
            beta = ao_expo_ordered_transp(m,i)
            c = c + NAI_pol_mult(A_center,B_center,power_A,power_B,  &
                alpha,beta,C_center,n_pt_in)                         &
                * ao_coef_normalized_ordered_transp(l,j)             &
                * ao_coef_normalized_ordered_transp(m,i)
          enddo
        enddo
        ao_integrals_coulomb_center(i,j,k) = -c
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
END_PROVIDER

BEGIN_PROVIDER [double precision, mo_integrals_coulomb_center, (mo_num,mo_num,n_coulomb_center)]
 implicit none
 BEGIN_DOC
! mo_integrals_coulomb_center(i,j,k) =
! $\langle \phi_i| -\frac{1}{|r-R_k|} | \phi_j \rangle$.
! where R_k is the coordinate of the k-th nucleus.
 END_DOC

 integer :: k
 mo_integrals_coulomb_center = 0.d0
 do k = 1, n_coulomb_center
   call ao_to_mo(                                                 &
       ao_integrals_coulomb_center(1,1,k),                        &
       size(ao_integrals_coulomb_center,1),                       &
       mo_integrals_coulomb_center(1,1,k),                        &
       size(mo_integrals_coulomb_center,1)                        &
       )
 enddo

END_PROVIDER

