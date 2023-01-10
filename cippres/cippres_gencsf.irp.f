
program cippres_gencsf
  use general
  implicit none
  BEGIN_DOC
! CIPPRES stands for Configuration Interaction Plugin for Photoionized and Resonant Electronic States.
! ORMAS-CI calculations are performed to compute bound and approximate continuum states. Hamiltonian and Dipole coupling matrix elements can then be computed between these states. Finally, Stieltjes imaging technique is applied to recover the correct continuum coupling matrix elements.
!  cippres_gencsf generates the lists of CSFs used in ORMAS-CI calculations 
  END_DOC

  character(len=lenmax) :: finput
  integer :: ilen, jlen
  logical :: file_e

  integer :: i, j

! TODO Save the lists of CSFs in EZFIO (ideally from the python script, i.e. without writing any txt files)

! READ XML file containing the ORMAS info
   print*, 'Reads ', finput_cippres
   print*, 'And generate the list of CSFs '

   finput=finput_cippres
   jlen=index(finput,' ')
   inquire( file="./"//finput(1:jlen-1), exist=file_e )
   if ( file_e .eqv. .false. ) then
    write(*,*) finput(1:jlen-1), " does not exist"
    stop
   endif
   call generate_csfs(finput)
   call ezfio_set_cippres_ifcsf(1)

end program cippres_gencsf
 
subroutine generate_csfs(finput)
 use general
 use bitmasks
 implicit none
 character(len=lenmax), intent(in) :: finput
 
 integer :: irun
 integer :: i,j,k
 integer :: nalpha, nbeta
 character(len=lenmax)             :: fileigvec

 integer, allocatable :: occ_a(:,:)
 integer, allocatable :: occ_b(:,:)
 character(1) :: ca, cb

!!nico integer :: nciruns, ncsfmax, ndetmax
 integer :: nciruns, ndetmax
 integer, dimension(:), allocatable :: nsta, ncsf
 integer, dimension(:,:), allocatable :: ndet
 double precision, dimension(:), allocatable :: prttol
 double precision, dimension(:,:,:), allocatable :: coefdet
 integer(bit_kind), allocatable :: csf_basis_tmp(:,:,:,:,:)

 integer :: n1, n2, n3, n4
 character(len=lenmax) :: fh, fl

! the python script will write the info and list of CSFs, for as many CIruns as given in finput, into header$irun.txt and list$irun.txt ($irun=1,2,3,....)
  call system('python $QP_ROOT/plugins/qp2_plugins_nsisourat/cippres/libs/generate_csfs.py '//trim(finput))

! save the info and list of CSFs into EZFIO
  open(unit=10,file='parser.txt') 
   read(10,*) j
   call ezfio_set_cippres_n_ciruns_cippres(j)
   nciruns = j
   read(10,*) j
   call ezfio_set_cippres_n_csf_max(j)
   ncsfmax = j
   read(10,*) j
   call ezfio_set_cippres_n_det_max_csf(j)
   ndetmax = j
  close(10)

 allocate(csf_basis_tmp(N_int,2, ndetmax, ncsfmax ,nciruns))
 allocate(occ_a(mo_num,ndetmax), occ_b(mo_num,ndetmax)) 
 allocate(nsta(nciruns),prttol(nciruns),ncsf(nciruns))
 allocate(ndet(ncsfmax,nciruns),coefdet(ndetmax,ncsfmax,nciruns))

 csf_basis_tmp(:,:,:,:,:) = 0_bit_kind

! read the info for each CI run
  do irun = 1, nciruns

   n1 = modulo(floor(irun/1000d0),10)
   n2 = modulo(floor(irun/100d0),10)
   n3 = modulo(floor(irun/10d0),10)
   n4 = modulo(irun,10)

   if(irun<10) then
     fh = 'header'//achar(48+irun)//'.txt'
     fl = 'list'//achar(48+irun)//'.txt'
   elseif(irun<100) then
     fh = 'header'//achar(48+n3)//achar(48+n4)//'.txt'
     fl = 'list'//achar(48+n3)//achar(48+n4)//'.txt'
   elseif(irun<1000) then
     fh = 'header'//achar(48+n2)//achar(48+n3)//achar(48+n4)//'.txt'
     fl = 'list'//achar(48+n2)//achar(48+n3)//achar(48+n4)//'.txt'
   elseif(irun<10000) then
     fh = 'header'//achar(48+n1)//achar(48+n2)//achar(48+n3)//achar(48+n4)//'.txt'
     fl = 'list'//achar(48+n1)//achar(48+n2)//achar(48+n3)//achar(48+n4)//'.txt'
   endif

   open(unit=21,file=fh)
   open(unit=22,file=fl)

    read(21,*)nsta(irun),prttol(irun)
    read(21,*)nalpha,nbeta
    if(nalpha/=elec_alpha_num) then
      print*,"Input inconsistent, nalpha in cirun",irun
    endif
    if(nbeta/=elec_beta_num) then
      print*,"Input inconsistent, nbeta in cirun",irun
    endif
    read(21,*)ncsf(irun)

   do i = 1, ncsf(irun)
    read(22,*)ndet(i,irun) 
    do j = 1, ndet(i,irun)
     read(22,*)coefdet(j,i,irun),ca,(occ_a(k,j),k=1,nalpha),cb,(occ_b(k,j),k=1,nbeta)
     call create_det(nalpha,nbeta,occ_a(1,j),occ_b(1,j),csf_basis_tmp(1,1,j,i,irun))
    enddo
   enddo 

   close(21)
   close(22)
  enddo

  call ezfio_set_cippres_coef_det_csf_basis(coefdet)
  call ezfio_set_cippres_n_csf_cippres(ncsf)
  call ezfio_set_cippres_n_sta_cippres(nsta)
  call ezfio_set_cippres_prttol_cippres(prttol)
  call ezfio_set_cippres_n_det_csf_cippres(ndet)
  call ezfio_set_cippres_csf_basis(csf_basis_tmp)

  deallocate(nsta,prttol,ncsf,ndet,coefdet,csf_basis_tmp)
  deallocate(occ_a,occ_b)

  PROVIDE ezfio_filename
  call system('mv parser.txt '//trim(ezfio_filename)//'/cippres/')
  call system('mv header*.txt '//trim(ezfio_filename)//'/cippres/')
  call system('mv list*.txt '//trim(ezfio_filename)//'/cippres/')

end subroutine generate_csfs

