BEGIN_PROVIDER [ double precision, pop_ao_per_mo, (ao_num,ao_num,mo_num) ]
&BEGIN_PROVIDER [double precision, gop_per_mo, (ao_num,mo_num)]
&BEGIN_PROVIDER [double precision, mulliken_population_per_mo, (nucl_num,mo_num)]

   BEGIN_DOC
   ! 
   END_DOC
   implicit none
   integer                        :: i,j,k,l
   double precision               :: dm_mo

   pop_ao_per_mo = 0.d0
   do k = 1, ao_num
     do l = 1, ao_num
       do i = 1, mo_num
           pop_ao_per_mo(l,k,i) += mo_coef(k,i) * mo_coef(l,i) * ao_overlap(k,l)
       enddo
     enddo
   enddo

   gop_per_mo = 0.d0
   do i = 1, ao_num
    do j = 1, ao_num
     do k = 1, mo_num
       gop_per_mo(i,k) += pop_ao_per_mo(j,i,k)
     enddo
    enddo
   enddo

   mulliken_population_per_mo = 0.d0
   do i = 1, ao_num
    do j = 1, mo_num
     mulliken_population_per_mo(ao_nucl(i),j) += gop_per_mo(i,j)
    enddo
   enddo

END_PROVIDER


BEGIN_PROVIDER [ double precision, one_e_dm_ao, (ao_num,ao_num,N_states) ]
   BEGIN_DOC
   ! 
   END_DOC
   implicit none
   integer                        :: i,j,k,l,istate
   double precision               :: dm_mo

   one_e_dm_ao = 0.d0
   do k = 1, ao_num
     do l = 1, ao_num
       do i = 1, mo_num
         do j = 1, mo_num
          do istate = 1, N_states
            dm_mo = myone_e_dm_mo(j,i,istate)
            !    if(dabs(dm_mo).le.1.d-10)cycle
            one_e_dm_ao(l,k,istate) += mo_coef(k,i) * mo_coef(l,j) * dm_mo
          enddo
         enddo
       enddo
     enddo
   enddo

END_PROVIDER

BEGIN_PROVIDER [ double precision, myone_e_dm_mo, (mo_num,mo_num,N_states) ]
   implicit none
   BEGIN_DOC
   ! One-body density matrix
   END_DOC
   myone_e_dm_mo = one_e_dm_mo_alpha + one_e_dm_mo_beta
END_PROVIDER

BEGIN_PROVIDER [double precision, my_population, (ao_num,ao_num,N_states)]
 implicit none
 integer :: i,j,istate
 BEGIN_DOC
! 
! 
 END_DOC
 my_population = 0.d0
 do i = 1, ao_num
  do j = 1, ao_num
   do istate = 1, N_states
     my_population(j,i,istate) = one_e_dm_ao(i,j,istate) * ao_overlap(i,j)
   enddo
  enddo
 enddo
END_PROVIDER

BEGIN_PROVIDER [double precision, pop, (ao_num,N_states)]
 implicit none
 pop = 0.d0
 integer :: i,j,istate
 BEGIN_DOC
! 
 END_DOC
 do i = 1, ao_num
  do j = 1, ao_num
   do istate = 1, N_states
     pop(i,istate) += my_population(j,i,istate)
   enddo
  enddo
 enddo

END_PROVIDER

BEGIN_PROVIDER [double precision, mulliken_population, (nucl_num,N_states)]
 implicit none
 integer :: i,j,istate
 BEGIN_DOC
!ATOMIC POPULATION
 END_DOC
 mulliken_population = 0.d0
 do i = 1, ao_num
  do istate = 1, N_states
    mulliken_population(ao_nucl(i),istate) += pop(i,istate)
  enddo
 enddo
END_PROVIDER

subroutine myprint_mulliken
 implicit none
 double precision :: accu
 integer :: i
 integer :: j
 integer :: istate
 do istate = 1, N_states
 print*,"*********************************************************************************"
 print*,'Mulliken population for state', istate
   accu= 0.d0
   do i = 1, nucl_num
    print*,i,nucl_charge(i),mulliken_population(i,istate)
    accu += mulliken_population(i,istate)
   enddo
 print*,'Sum of Mulliken Pop. = ',accu
 print*,'AO POPULATIONS'
 accu = 0.d0
 do i = 1, ao_num
  accu += pop(i,istate)
  write(*,'(1X,I3,1X,A4,1X,I2,1X,A4,1X,F10.7)')i,trim(element_name(int(nucl_charge(ao_nucl(i))))),ao_nucl(i),trim(l_to_character(ao_l(i))),pop(i,istate)
 enddo
 print*,'sum = ',accu
 print*,"*********************************************************************************"
 enddo

end

subroutine myprint_mulliken_per_mo
 implicit none
 double precision :: accu
 integer :: i
 integer :: j
 integer :: imo
 do imo = 1, mo_num
 print*,"*********************************************************************************"
 print*,'Mulliken population for MO', imo
   accu= 0.d0
   do i = 1, nucl_num
    print*,i,trim(element_name(int(nucl_charge(i)))),nucl_charge(i),mulliken_population_per_mo(i,imo)
    accu += mulliken_population_per_mo(i,imo)
   enddo
 print*,'Sum of Mulliken Pop. = ',accu
 print*,'AO POPULATIONS'
 accu = 0.d0
 do i = 1, ao_num
  accu += pop(i,imo)
  write(*,'(1X,I3,1X,A4,1X,I2,1X,A4,1X,F10.7)')i,trim(element_name(int(nucl_charge(ao_nucl(i))))),ao_nucl(i),trim(l_to_character(ao_l(i))),gop_per_mo(i,imo)
 enddo
 print*,'sum = ',accu
 print*,"*********************************************************************************"
 enddo

end

