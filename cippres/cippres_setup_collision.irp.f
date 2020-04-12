program cippres_setup_collision
  use general
 ! create a routine for one-e general matrix with determinants

  implicit none
  BEGIN_DOC
! CIPPRES stands for Configuration Interaction Plugin for Photoionized and Resonant Electronic States.
! cippres_collision computes the V_projectile matrice couplings between the CI eigenvectors of ici1=1 run
  END_DOC

  character(len=lenmax) :: finput
  integer :: ilen, jlen, ijob
  logical :: file_e

  print*, 'Reads ', finput_coll

  finput=finput_coll
  jlen=index(finput,' ')
  inquire( file="./"//finput(1:jlen-1), exist=file_e )
  if ( file_e .eqv. .false. ) then
   write(*,*) finput(1:jlen-1), " does not exist"
   stop
  endif

  print*,"SET UP"
  call setup_coll(finput)
  print*,"SET UP DONE"

end program cippres_setup_collision

subroutine setup_coll(finput)
 use general
 use bitmasks
 implicit none

 character(len=lenmax) :: finput
 integer :: i, j
 double precision :: vproj ! along z
 double precision :: bmin, bmax
 integer :: nb
 character(len=lenmax) :: btype
 double precision, allocatable :: btmp(:)

 double precision :: zmin, zmax, spac
 integer :: nz
 character(len=lenmax) :: ztype
 double precision, allocatable :: ztmp(:), ttmp(:)

 integer :: npcenter
 double precision, allocatable :: ch_pcenter(:)

 integer :: istate, istamin, istamax, sistamin, sistamax, distamin, distamax
 logical :: exists

 call system('python $QP_ROOT/plugins/qp2_plugins_nsisourat/cippres/libs/coll_input.py '//trim(finput)//' > collinp')

 open(unit=10,file='collinp')
  read(10,*)btype,bmin,bmax,nb
  read(10,*)vproj
  read(10,*)ztype,zmin,zmax,nz
  read(10,*)istamin,istamax
  read(10,*)sistamin,sistamax
  read(10,*)distamin,distamax
  read(10,*)istate
  read(10,*)npcenter
  allocate(ch_pcenter(npcenter))
  do i = 1, npcenter
    read(10,*)ch_pcenter(i)
  enddo
 close(10)
 call ezfio_set_cippres_n_pcenter(npcenter)
 call ezfio_set_cippres_charge_pcenter(ch_pcenter)

 print*,nb,nz,npcenter
 call ezfio_set_cippres_v_coll(vproj)
 call ezfio_set_cippres_stamin_bound(istamin)
 call ezfio_set_cippres_stamax_bound(istamax)
 call ezfio_set_cippres_stamin_si(sistamin)
 call ezfio_set_cippres_stamax_si(sistamax)
 call ezfio_set_cippres_stamin_di(distamin)
 call ezfio_set_cippres_stamax_di(distamax)
 call ezfio_set_cippres_i_state_coll(istate)

 call ezfio_set_cippres_n_time(nz)
 allocate(ztmp(nz),ttmp(nz))

 if(ztype=='linear' .or. ztype=='Linear' .or. ztype=='LINEAR') then

   do i = 1, nz
     ztmp(i) = zmin + (i-1)*(zmax-zmin)/nz
   enddo

 elseif(ztype=='exp' .or. ztype=='Exp' .or. ztype=='EXP') then

  spac=0d0
  do i = 1, nz/2
   spac=spac+1.1**i
  enddo
   spac=(zmax-zmin)/(2d0*spac)

  ztmp(nz/2)=0d0
   j = 0
  do i = nz/2+1, nz-1
    j = j + 1
    ztmp(i) = ztmp(i-1)+1.1**j*spac
    ztmp(nz-i) = -ztmp(i)
  enddo
    ztmp(nz) = zmax

 else
   write(*,*)"No other z grids implemented"
   stop
 endif

 do i = 1, nz
   ttmp(i) = ztmp(i)/vproj
 enddo

 call ezfio_set_cippres_n_bimp(nb)
 allocate(btmp(nb))

 if(btype=='linear' .or. btype=='Linear' .or. btype=='LINEAR') then
 
  do i = 1, nb
    btmp(i) = bmin + (i-1)*(bmax-bmin)/nb
  enddo
 
 elseif(btype=='exp' .or. btype=='Exp' .or. btype=='EXP') then
 
  spac=0d0
  do i = 1, nb/2
   spac=spac+1.1**i
  enddo
   spac=(bmax-bmin)/(2d0*spac)
 
  btmp(nb/2)=0d0
   j = 0
  do i = nb/2+1, nb-1
    j = j + 1
    btmp(i) = btmp(i-1)+1.1**j*spac
    btmp(nb-i) = -btmp(i)
  enddo
    btmp(nb) = bmax
 
 else
   write(*,*)"No other b grids implemented"
   stop
 endif
 
 call ezfio_set_cippres_tgrid(ttmp)
 call ezfio_set_cippres_zgrid(ztmp)
 call ezfio_set_cippres_bgrid(btmp)

 deallocate(ztmp,btmp,ttmp,ch_pcenter)

end subroutine setup_coll



