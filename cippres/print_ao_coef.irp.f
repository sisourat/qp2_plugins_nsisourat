program print_ao_coef
  implicit none
  integer :: i,j

   do i=1,ao_num
     do j=1,ao_prim_num(i)
      print*, ao_coef_normalized(i,j)
     enddo
   enddo

end


