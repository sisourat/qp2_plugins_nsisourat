module general
implicit none

integer, parameter :: lenmax = 200
double precision, parameter :: pi = dacos(-1d0)
double complex, parameter :: imag = dcmplx(0d0,1d0)

end module general

module SlaterDeterminant
implicit none

type Sdeterminant
 integer :: nunpairel, nalpha, nbeta
 integer, dimension(:), allocatable :: alpha, beta
end type Sdeterminant

end module SlaterDeterminant
