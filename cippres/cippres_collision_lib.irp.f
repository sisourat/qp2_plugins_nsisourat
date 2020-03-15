use bitmasks ! you need to include the bitmasks_module.f90 features
use general

 BEGIN_PROVIDER [integer, nb_b]
  implicit none
    call ezfio_get_cippres_nb_b(nb_b)
 END_PROVIDER 

 BEGIN_PROVIDER [double precision, coll_couplings_cippres, (n_csf_max,n_csf_max,nb_b)]
 use general
 implicit none
 integer :: i, j, k, l
 integer :: ib, ic
 double precision, allocatable :: eigval1(:),eigvec1(:,:),eigval2(:),eigvec2(:,:),coll_csf_mat(:,:),coll_mat(:,:,:)
 double precision, dimension(mo_num,mo_num) :: w1e

 double precision :: hij

 logical :: exists

   coll_couplings_cippres(:,:,:) = 0d0

   allocate(coll_csf_mat(n_csf_cippres(ici2),n_csf_cippres(ici1)))
 
   allocate(eigval1(n_csf_cippres(ici1)),eigval2(n_csf_cippres(ici2)))
   eigval1(:) = eigvalues_cippres(1:n_csf_cippres(ici1),ici1)
   eigval2(:) = eigvalues_cippres(1:n_csf_cippres(ici2),ici2)

   allocate(eigvec1(n_csf_cippres(ici1),n_csf_cippres(ici1)))
   allocate(eigvec2(n_csf_cippres(ici2),n_csf_cippres(ici2)))
   eigvec1(:,:) = eigvectors_cippres(1:n_csf_cippres(ici1),1:n_csf_cippres(ici1),ici1)
   eigvec2(:,:) = eigvectors_cippres(1:n_csf_cippres(ici2),1:n_csf_cippres(ici2),ici2)

   allocate(coll_mat(n_csf_cippres(ici2),n_csf_cippres(ici1),nb_b))

   coll_mat(:,:,:) = 0d0

   do ib = 1, nb_b
    w1e(:,:) = 0d0
    do ic = 1, n_coulomb_center
      coulomb_center(1,ic) = 0d0 !b(ib)
      coulomb_center(2,ic) = 0d0 
      coulomb_center(3,ic) = 0d0 ! z(time)
      w1e(:,:) += mo_integrals_coulomb_center(:,:,ic)
    enddo

    coll_csf_mat(:,:) = 0d0
    do i = 1, n_csf_cippres(ici1) ! first loop on the csf of the space ispace 
     do j = 1, n_csf_cippres(ici2)
      do k = 1, n_det_csf_cippres(i,ici1) ! then on the determinants belonging to the ith CSF of space ispace
       do l = 1, n_det_csf_cippres(j,ici2)
          call i_w1e_j(csf_basis(1,1,k,i,ici1),csf_basis(1,1,l,j,ici2),N_int,w1e,hij)
          coll_csf_mat(j,i) += hij * coef_det_csf_basis(k,i,ici1) * coef_det_csf_basis(l,j,ici2)
       enddo
      enddo
     enddo
    enddo

    do i = 1, n_csf_cippres(ici1) ! first loop on the first eigenvectors
     do j = 1, n_csf_cippres(ici2) ! then on the second eigenvectors
      do k = 1, n_csf_cippres(ici1) ! loop over the csfs of the ici1 run
       do l = 1, n_csf_cippres(ici2) ! then over the csfs of the ici2 run
          coll_mat(j,i,ib) += coll_csf_mat(l,k) * eigvec1(k,i) * eigvec2(l,j)
       enddo
      enddo
     enddo
    enddo

   enddo

  call ezfio_set_cippres_coll_couplings_cippres(coll_mat)

 deallocate(coll_csf_mat,eigval1,eigval2,eigvec1,eigvec2,coll_mat) 

 END_PROVIDER

