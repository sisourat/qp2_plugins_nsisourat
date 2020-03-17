! DO NOT MODIFY BY HAND
! Created by $QP_ROOT/scripts/ezfio_interface/ei_handler.py
! from file /home/nico/Workspace/qp2/src/cippres/EZFIO.cfg


BEGIN_PROVIDER [ character*(32), finput_coll  ]
  implicit none
  BEGIN_DOC
! xml file containing collision input info
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_cippres_finput_coll(has)
    if (has) then
      write(6,'(A)') '.. >>>>> [ IO READ: finput_coll ] <<<<< ..'
      call ezfio_get_cippres_finput_coll(finput_coll)
    else
      print *, 'cippres/finput_coll not found in EZFIO file'
      stop 1
    endif
  endif
  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( finput_coll, 1*32, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read finput_coll with MPI'
    endif
  IRP_ENDIF

  call write_time(6)

END_PROVIDER

BEGIN_PROVIDER [ integer, i_stj_job  ]
  implicit none
  BEGIN_DOC
! if i_stj_job is 1 then cippres_stieltjes performs Fano calculations (for ifanosta state), if it is 2 cippres_stieltjes performs dipole calculations (for idipsta state)
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_cippres_i_stj_job(has)
    if (has) then
      write(6,'(A)') '.. >>>>> [ IO READ: i_stj_job ] <<<<< ..'
      call ezfio_get_cippres_i_stj_job(i_stj_job)
    else
      print *, 'cippres/i_stj_job not found in EZFIO file'
      stop 1
    endif
  endif
  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( i_stj_job, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read i_stj_job with MPI'
    endif
  IRP_ENDIF

  call write_time(6)

END_PROVIDER

BEGIN_PROVIDER [ integer, ifcsf  ]
  implicit none
  BEGIN_DOC
! if ifcsf eq 1 then some CSFs were generated before, otherwise one should run cippres to generate the lists of CSFs first
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_cippres_ifcsf(has)
    if (has) then
      write(6,'(A)') '.. >>>>> [ IO READ: ifcsf ] <<<<< ..'
      call ezfio_get_cippres_ifcsf(ifcsf)
    else
      print *, 'cippres/ifcsf not found in EZFIO file'
      stop 1
    endif
  endif
  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( ifcsf, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read ifcsf with MPI'
    endif
  IRP_ENDIF

  call write_time(6)

END_PROVIDER

BEGIN_PROVIDER [ character*(32), finput_cippres  ]
  implicit none
  BEGIN_DOC
! xml file containing ormas info
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    
    call ezfio_has_cippres_finput_cippres(has)
    if (has) then
      write(6,'(A)') '.. >>>>> [ IO READ: finput_cippres ] <<<<< ..'
      call ezfio_get_cippres_finput_cippres(finput_cippres)
    else
      print *, 'cippres/finput_cippres not found in EZFIO file'
      stop 1
    endif
  endif
  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( finput_cippres, 1*32, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read finput_cippres with MPI'
    endif
  IRP_ENDIF

  call write_time(6)

END_PROVIDER
