! the subroutine nonlinearSolver() solves nonlinear system of functions F(u)=Lu+N(u)-f=0
! by two steps iteration:
! 1: solve v_n in F'(u_n)v_n=-F(u_n)+r_n given |r_n|/|F(u_n)|<=\eta_n
! 2: update u_{n+1}=u_n+v_n
! where u_n is the solution at each step, r_n is the residual, F'=L+N'(u) the jocobian matrix

subroutine nonlinearSolver
use bicg
use pbedata
use comdata
implicit none
real*8, allocatable, dimension(:) ::  resdu, du
real*8 errInner, tol, eta, alpha, residueNorm
integer iErr, nDim, iDisplay, nItr 

! allocate residual and charge of u
allocate(resdu(nx*ny*nz),du(nx*ny*nz),stat=ierr)
if (ierr .ne. 0) then
    write(*,*) 'Error allocating resdu, du'
    stop
endif

nDim=nx*ny*nz
eta=1.d-2     ! "forcing term" for solving F'(u_n)v_n=-F(u_n)+r_n, |r_n|/|F(u_n)|<=\eta_n
tol=1.d-6     ! tolerance for inexact newton
errInner=1.d0 ! initial error
iDisplay=1    ! output error at each Newton iteration (1: Yes; 0: No)    
nItr=0        ! number of iteration for Newtwon

do while(errInner.GT.tol)
    ! resdu --> du --> updated du thus u=u+du --> resdu until tol is met
    resdu=0.d0; du=0.d0
    nItr=nItr+1
    call seekResidue(resdu) ! compute residual of the Jacobian linear sys. r_n=F'(u_n)v_n+F(u_n)
    !implement alternative TST(r^n) in algorithm 7 [Holst94]
    !residueNorm=maxval(abs(resdu))
    !eta=0.01d0*residueNorm*residueNorm
    !print *,eta,residueNorm
    !if (eta>residueNorm) eta=residueNorm 
    call seekDirection(ndim,resdu,du,eta)
    alpha=1.d0 !xxxx
    !call lineSearch(alpha, resdu, du)
    errInner=MAXVAL(ABS(du))
    if(iDisplay.EQ.1) then
        print *,'---Inner Error: ', nItr, errInner
    endif
    u=u+alpha*du
enddo

deallocate(resdu, du, stat=ierr)
if (ierr .ne. 0) then
    write(*,*) 'Error deallocating resdu, du'
    stop
endif

end subroutine nonlinearSolver

!------------------------------------  
! compute residual: r=b-Lu-N(u), essentially r=F(u)
! here N(u)=-cappa^2*sinh(u)
subroutine seekResidue(residue)
use bicg
use pbedata
use comdata
implicit double precision(a-h, o-z)

real*8 :: residue(nx*ny*nz), Au(nx*ny*nz)
integer :: i,j,k,ijk,ndim

!print *,'computing residual...'
ndim=nx*ny*nz
Au=0.d0
call atimes(ndim, u, Au, 0)
!print *,'b-Au', maxval(abs(bftc-Au))

residue=bftc-Au
!print *,'finishing r=b-Au ...'

do k=2,nz-1     !handle the nonlinear term
    do j=2,ny-1
        do i=2,nx-1
            if(io(i,j,k).EQ.1) then
                ijk=id3d(i,j,k)
                residue(ijk)=residue(ijk)+cappa_nl*sinh(u(ijk))
                !print *,i,j,k,real(u(ijk)),real(sinh(u(ijk)))
            endif
        enddo
    enddo
enddo
!print *,'finished computing residual...'
end subroutine seekResidue

!--------------------------------------------
! solve v_n in F'(u_n)v_n=-F(u_n)+r_n
! here sol=v_n
! F'(u_n)=A+N'(u_n)
! bvector=-F(u_n)+r_n
subroutine seekDirection(ndim,bvector,sol,tol)
use bicg
use pbedata
use comdata
implicit none 
integer, intent(in) :: ndim
real*8, intent(in) :: tol
real*8, intent(in), dimension(ndim) :: bvector
real*8, intent(out),dimension(ndim) :: sol

integer iter, itol, itmax, ijk, ix, jy, kz
real*8 :: err

do ijk=1,ndim   
    call line2cube(ix,jy,kz,ijk)
    if(ix.NE.1.AND.ix.NE.nx.AND.jy.NE.1.AND.jy.NE.ny.AND.kz.NE.1.AND.kz.NE.nz) then
        if(io(ix,jy,kz).EQ.1) sa(ijk)=sa(ijk)-cappa_nl*cosh(u(ijk))
    endif
enddo

sol=0.D0
itol=2          ! converge criterion measured by 2-norm
itmax = 1000    ! max iteration number allowed

call linbcg(ndim,bvector,sol,itol,tol,itmax,iter,err)
print *,iter,err
sa=sb

end subroutine seekDirection

!--------------------------------------------
subroutine lineSearch(alpha,residue,deltau)
use bicg
use comdata
use pbedata
implicit none 
real*8 :: residue(nx*ny*nz), deltau(nx*ny*nz), u_tmp(nx*ny*nz)
integer NNN, i
real*8 dh, alpha, residue_norm, residue_norm_new, tol


tol=5.d-10
residue_norm=MAXVAL(ABS(residue))
residue=0.D0
NNN=100;dh=1.D0/NNN
do i=NNN,1,-1
    alpha=dh*i
    u_tmp=u
    u=u_tmp+alpha*deltau
    call seekResidue(residue)
    residue_norm_new=MAXVAL(ABS(residue))
    if(residue_norm_new.LT.residue_norm+tol) then
        print *,'alpha= ',alpha
        goto 1234
    endif
enddo
print *,'Fail to find the steplength',residue_norm_new,residue_norm
stop
!1234 u=u_tmp
1234 u=u_tmp 
end subroutine lineSearch
