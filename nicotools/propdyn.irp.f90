module splineinterp

integer :: ksave

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      REAL*8 yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=10000)
      INTEGER i,k
      REAL*8 p,qn,sig,un,u(NMAX)
      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+   &
     1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig* &
     u(i-1))/p
11    continue
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      END SUBROUTINE spline 
!  (C) Copr. 1986-92 Numerical Recipes Software #>.)@1.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE splint(klo,khi,xa,rya,cya,ry2a,cy2a,n,x,ry,cy)
      INTEGER n
      REAL*8 x,ry,cy,xa(n),ry2a(n),rya(n),cy2a(n),cya(n)
      INTEGER k,khi,klo
      REAL*8 a,b,h

      h=xa(khi)-xa(klo)
      if (h.eq.0.) pause 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      ry=a*rya(klo)+b*rya(khi)+((a**3-a)*ry2a(klo)+(b**3-b)*ry2a(khi))*(h**2)/6.
      cy=a*cya(klo)+b*cya(khi)+((a**3-a)*cy2a(klo)+(b**3-b)*cy2a(khi))*(h**2)/6.
     return

      END SUBROUTINE splint
!  (C) Copr. 1986-92 Numerical Recipes Software #>.)@1.

end module splineinterp


module linearinterp
implicit none

   contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine interp(i,xa,ya,n,x,ry,cy)
implicit none

integer, intent(in) :: n
double precision, dimension(n), intent(in) :: xa
double complex, dimension(n), intent(in) :: ya
double precision :: x, ry, cy 

double precision :: a
integer :: i, j

!! Linear interpolation

! searching for two closer points of x in xa
! xa must be ascending sorted

 if(x .lt. xa(2) .or. x .gt. xa(n-1)) then
  write(*,*)'x out of range in interp, I stop', x
  stop
 endif

! call hunt(xa,n,x,i)
! call locate(xa,n,x,i)

 a = (real(ya(i+1)) - real(ya(i)))/(xa(i+1) - xa(i))
 ry = real(ya(i)) + a*(x-xa(i))

 a = (aimag(ya(i+1)) - aimag(ya(i)))/(xa(i+1) - xa(i))
 cy = aimag(ya(i)) + a*(x-xa(i))

end subroutine interp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE hunt(xx,n,x,jlo)
      INTEGER jlo,n
      REAL*8 x,xx(n)
      INTEGER inc,jhi,jm
      LOGICAL ascnd
      ascnd=xx(n).gt.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
        jlo=0
        jhi=n+1
        goto 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then
1       jhi=jlo+inc
        if(jhi.gt.n)then
          jhi=n+1
        else if(x.ge.xx(jhi).eqv.ascnd)then
          jlo=jhi
          inc=inc+inc
          goto 1
        endif
      else
        jhi=jlo
2       jlo=jhi-inc
        if(jlo.lt.1)then
          jlo=0
        else if(x.lt.xx(jlo).eqv.ascnd)then
          jhi=jlo
          inc=inc+inc
          goto 2
        endif
      endif
3     if(jhi-jlo.eq.1)return
      jm=(jhi+jlo)/2
      if(x.gt.xx(jm).eqv.ascnd)then
        jlo=jm
      else
        jhi=jm
      endif
      goto 3
      END SUBROUTINE hunt
!  (C) Copr. 1986-92 Numerical Recipes Software #>.)@1.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE locate(xx,n,x,j)
      INTEGER j,n
      REAL*8 x,xx(n)
      INTEGER jl,jm,ju
      jl=0
      ju=n+1
10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      goto 10
      endif
      j=jl
      return
      END SUBROUTINE locate
!  (C) Copr. 1986-92 Numerical Recipes Software #>.)@1.

end module linearinterp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module propdyn
use splineinterp
use linearinterp
implicit none

double complex, dimension(:,:), allocatable :: matintrp
double precision, dimension(:,:,:), allocatable :: rmat2intrp, cmat2intrp
double complex, dimension(:), allocatable :: psi
double precision, dimension(:,:,:), allocatable ::  mat
double complex, dimension(:,:,:), allocatable ::  mcoup

integer :: ntime, ntotsta, nsi, ndi, nsta
double precision, dimension(:), allocatable :: timegrid
double precision, dimension(:), allocatable :: esta

double precision, dimension(:,:,:), allocatable :: g_si_intrp, g_ddi_intrp
double precision, dimension(:,:,:), allocatable :: g_si, g_ddi
double precision, dimension(:,:), allocatable :: g_sdi_intrp
double precision, dimension(:,:), allocatable :: g_sdi
double precision :: p_si, p_ddi, p_sdi

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine hpsit(time,psiin,psiout,psidim,lhpsi,ihpsi,rhpsi,chpsi)
implicit none

logical    lhpsi(*)
integer    psidim, ihpsi(*)
real*8     rhpsi(*), time
complex*16 psiin(psidim), psiout(psidim), chpsi(*)

double precision :: rmat, cmat
double complex, parameter :: imag = dcmplx(0d0,1d0)

integer :: i, j, ista, jsta
integer :: klo, khi, k

! --- INTERPOLATION OF MCOUP MATRIX
      if(timegrid(ksave).lt. time .and. timegrid(ksave+1).gt. time) then
       klo=ksave
       khi=ksave+1
      else
       klo=1
       khi=ntime
1      if (khi-klo.gt.1) then
         k=(khi+klo)/2
         if(timegrid(k).gt.time)then
           khi=k
         else
           klo=k
         endif
       goto 1
       endif
       ksave = klo
      endif

matintrp(:,:)=0d0
!!!$OMP PARALLEL DO COLLAPSE(2)  PRIVATE(i,j,rmat,cmat)
!$OMP PARALLEL DO PRIVATE(i,j,rmat,cmat)
! SCHEDULE(DYNAMIC) 
do i = 1, nsta
 do j = 1, nsta
!   call interp(timegrid%a,mcoup(:,j,i),timegrid%na,time,rmat,cmat)
   call splint(klo,khi,timegrid,real(mcoup(:,j,i)),aimag(mcoup(:,j,i)),rmat2intrp(:,j,i),cmat2intrp(:,j,i),ntime,time,rmat,cmat)
   matintrp(j,i) = matintrp(j,i) +dcmplx(rmat,cmat)*exp(-imag*(esta(i)-esta(j))*time)

!   call splint(klo,khi,timegrid,g_si(:,j,i),g_si(:,j,i),g_si_intrp(:,j,i),g_si_intrp(:,j,i),ntime,time,rmat,cmat)
!   matintrp(j,i) = matintrp(j,i) - imag*rmat
!   call splint(klo,khi,timegrid,g_ddi(:,j,i),g_ddi(:,j,i),g_ddi_intrp(:,j,i),g_ddi_intrp(:,j,i),ntime,time,rmat,cmat)
!   matintrp(j,i) = matintrp(j,i) - imag*rmat
!   matintrp(j,i) = matintrp(j,i)*exp(-imag*(esta(i)-esta(j))*time)

 enddo
   call splint(klo,khi,timegrid,g_si(:,i,i),g_si(:,i,i),g_si_intrp(:,i,i),g_si_intrp(:,i,i),ntime,time,rmat,cmat)
   matintrp(i,i) = matintrp(i,i) - imag*rmat
enddo
!$OMP END PARALLEL DO

! --- COMPUTE ACTION OF HAMILTONIAN ON WAVEFUNCTION ---
 psiout = -dcmplx(0d0,1d0)*matmul(matintrp,psiin)

! write(60,'(60(f15.6,1X))')time,mat(1,1)*dcmplx(0d0,-1d0)/vproj
! write(*,'(60(f15.6,1X))')time,(cdabs(psiin(i))**2,i=1,nsta)


end subroutine hpsit


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine dyn

double complex, dimension(:), allocatable :: psit, dtpsit
double precision :: machprec
integer :: psidim

double precision :: time

double precision :: IntPeriod,AbsTime,LargeStep,TolError,ActualLargeStep,NextLargeStep, Step
integer :: IntOrder,SmallSteps,ErrorCode
external  AbsBSError,PolyExtrapol

double complex, dimension(:,:), allocatable :: AuxPsi !(psidim,intorder+2)
double precision, dimension(:), allocatable :: RData
double complex, dimension(:), allocatable :: CData
integer, dimension(:), allocatable :: IData
logical, dimension(:), allocatable :: LData

logical :: RestartABM
integer :: Steps, RepeatedSteps
double precision :: InitStep
External   AbsABMError

integer :: i, j

 psidim = nsta
 IntOrder = 6
 TolError = 1d-09

 allocate(Psit(nsta),dtPsit(nsta))
 allocate(RData(nsta),CData(nsta),IData(nsta),LData(nsta))
 allocate(AuxPsi(psidim,IntOrder+2))

! for spline interpolation
ksave = 1
do i = 1, nsta
 do j = 1, nsta
   call spline(timegrid,real(mcoup(:,j,i)),ntime,0d0,0d0,rmat2intrp(:,j,i))
   call spline(timegrid,aimag(mcoup(:,j,i)),ntime,0d0,0d0,cmat2intrp(:,j,i))
   call spline(timegrid,g_si(:,j,i),ntime,0d0,0d0,g_si_intrp(:,j,i))
   call spline(timegrid,g_ddi(:,j,i),ntime,0d0,0d0,g_ddi_intrp(:,j,i))
 enddo
   call spline(timegrid,g_sdi(:,i),ntime,0d0,0d0,g_sdi_intrp(:,i))
enddo

! psi is declared in module general and must be given in main with the proper initial conditions
 Psit(:) = psi(:)

 RestartABM = .true.
 InitStep = abs(timegrid(2)-timegrid(1))
 Abstime = timegrid(2)
 IntPeriod = (timegrid(ntime)-timegrid(5))

 call hpsit(Abstime,Psit,DtPsit,PsiDim,LData,IData,RData,CData)
 call  ABM (Psit,DtPsit,PsiDim,IntPeriod,AbsTime,IntOrder,              &
                           InitStep,TolError,RestartABM,Steps,          &
                           RepeatedSteps,ErrorCode,AuxPsi,              &
                           hpsit,AbsABMError,CData,RData,               &
                           IData,LData)

 !call rk4(hpsit, Abstime, Psit, 0.005d0, PsiDim, timegrid(ntime-3), lData,iData,rData,cData) 

 Abstime = Abstime + IntPeriod

 if(ErrorCode /= 0) then
   write(*,*)'Error in ABM, I stop', ErrorCode
   stop
 endif

 psi(:) = Psit(:)

! deallocate(mat)
 deallocate(RData,CData,IData,LData)
 deallocate(AuxPsi)
 deallocate(Psit,dtPsit)

end subroutine dyn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine rk4 (f,t,x,h,n,tmax,lhpsi,ihpsi,rhpsi,chpsi)   
        integer, intent (in) :: n
        double complex, dimension(n) :: f1, f2, f3, f4, ta
        double precision :: t, tmax, h
        double complex, dimension(n) , intent(inout) :: x

        logical    lhpsi(*)
        integer    psidim, ihpsi(*)
        real*8     rhpsi(*)
        complex*16 chpsi(*)
        
        external f

        do while(t<tmax)
           call f(t, x, f1, n, lhpsi,ihpsi,rhpsi,chpsi)       
           call f(t + h/2.0, x + 0.5*h*f1, f2, n, lhpsi,ihpsi,rhpsi,chpsi)   
           call f(t + h/2.0, x + 0.5*h*f2, f3, n, lhpsi,ihpsi,rhpsi,chpsi)   
           call f(t + h, x + h*f3, f4, n, lhpsi,ihpsi,rhpsi,chpsi)
           x = x + h*(f1 + 2*f2  +2*f3 + f4)/6.0
           t = t + h
        end do

end subroutine rk4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module propdyn

