program print_mo_overlap
  implicit none
  integer :: i,j,k,l
  double precision :: overlap(mo_num,mo_num), maxovl
  integer :: list_mo_del(mo_num), n_mo_del, skip

  integer :: m,p,s
  integer :: i1,i2
  double precision :: x, cutoff

  cutoff = 1d-12

  m = size(mo_coef,1)
  p = size(mo_overlap,1)
  call ortho_qr(mo_overlap,p,mo_num,mo_coef,m,ao_num)

!  call orthonormalize_mos
  
  do i=1,mo_num
    do j=1,mo_num
      do k=1,ao_num
        do l=1,ao_num

          overlap(j,i) += mo_coef(l,j) * mo_coef(k,i) * ao_overlap(k,l)

        enddo
      enddo
    enddo
  enddo

  
  n_mo_del = 0
  list_mo_del(:) = 0
  do i=1,mo_num-1
    skip = 0
    do k = 1, n_mo_del
     if(i==list_mo_del(k)) then
      skip = 1
      exit
     endif
    enddo
    if(skip==1) cycle
    do j=i+1,mo_num
     skip = 0
     do k = 1, n_mo_del
      if(j==list_mo_del(k)) then
         skip = 1
         exit
      endif
     enddo
     if(skip==1) cycle


     if(overlap(j,i)>cutoff) then
      n_mo_del = n_mo_del + 1  
      list_mo_del(n_mo_del) = j
!      print*, j, i
     endif
    ! print*, overlap(j,i) / sqrt( overlap(j,j)*overlap(i,i) ) , j, i
    !  print*, overlap(j,i) 
    enddo
  enddo

  print*,"Nb of AOs", ao_num
  print*,"Nb of dropped MOs", n_mo_del 
  print*,"Nb of MOs to be used", mo_num-n_mo_del 
  print*," "
  print*,"The following MOs will be swapped at the end of the MOs list"
  do k = 1, n_mo_del
   i1 = mo_num + 1 - k
   i2 = list_mo_del(k)
   do i=1,ao_num
     x = mo_coef(i,i1)
     mo_coef(i,i1) = mo_coef(i,i2)
     mo_coef(i,i2) = x
   enddo
   print*, list_mo_del(k), "=>", i1
  enddo

  call save_mos
end


