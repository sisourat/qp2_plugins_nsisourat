program cippres_collision
  use general
 ! create a routine for one-e general matrix with determinants

  implicit none
  BEGIN_DOC
! CIPPRES stands for Configuration Interaction Plugin for Photoionized and Resonant Electronic States.
! cippres_collision computes the dipole matrice couplings between the CI eigenvectors of ici1 and ici2 runs
  END_DOC

  print*, "coll = ", coll_couplings_cippres

end program cippres_collision


