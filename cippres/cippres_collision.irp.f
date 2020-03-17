program cippres_collision
  use general
 ! create a routine for one-e general matrix with determinants

  implicit none
  BEGIN_DOC
! CIPPRES stands for Configuration Interaction Plugin for Photoionized and Resonant Electronic States.
! cippres_collision computes the V_projectile matrice couplings between the CI eigenvectors of ici1=1 run
  END_DOC

  character(len=lenmax) :: finput
  integer :: ilen, jlen, ijob
  logical :: file_e, exists


  print*, 'Reads ', finput_coll

  finput=finput_coll
  jlen=index(finput,' ')
  inquire( file="./"//finput(1:jlen-1), exist=file_e )
  if ( file_e .eqv. .false. ) then
   write(*,*) finput(1:jlen-1), " does not exist"
   stop
  endif

  print*,"COMPUTING VP MATRIX ELEMENTS"
  call ezfio_set_cippres_coll_couplings_cippres(coll_couplings_cippres)
  print*,"JOB DONE"

end program cippres_collision

