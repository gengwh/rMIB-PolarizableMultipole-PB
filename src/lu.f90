!
!--------------------------------------------------------------------------------
!
SUBROUTINE ludcmp(a,n,np,indx,d)
IMPLICIT DOUBLE PRECISION(A-H,O-Z)
DIMENSION indx(n),a(np,np)
PARAMETER (TINY=1.0e-20)
!-------------------------------------------------------------------------
!Given a matrix a(1:n,1:n), with physical dimension np by np, this routine
!replaces it by the LU decomposition of a rowwise permutation of itself.
!a and n are input. a is output, arranged as in equation (2.3.14) above;
!indx(1:n) is an output vector that records the row permutation e^Kffected
!by the partial pivoting; d is output as ^F+-1 depending on whether the
!number of row interchanges was even or odd, respectively. This routine is
!used in combination with lubksb to solve linear equations or invert a matrix.
!----------------------------------------------------------------------------
DIMENSION vv(n)
d=1.D0
do i=1,n
    aamax=0.D0
    do j=1,n
        if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
    enddo
    if (aamax.eq.0.) then 
		print *,transpose(a)
		write(*,*) 'singular matrix in ludcmp'
	endif
    vv(i)=1.D0/aamax
enddo
do j=1,n
    do i=1,j-1
        sum=a(i,j)
        do k=1,i-1
            sum=sum-a(i,k)*a(k,j)
        enddo
        a(i,j)=sum
    enddo
    aamax=0.D0
    do i=j,n
        sum=a(i,j)
        do k=1,j-1
            sum=sum-a(i,k)*a(k,j)
        enddo
        a(i,j)=sum
        dum=vv(i)*abs(sum)
        if (dum.ge.aamax) then
             imax=i
             aamax=dum
        endif
    enddo
    if (j.ne.imax)then
        do k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
        enddo
        d=-d
        vv(imax)=vv(j)
    endif
    indx(j)=imax
    if(a(j,j).eq.0.)a(j,j)=TINY
    if(j.ne.n)then
        dum=1.D0/a(j,j)
        do i=j+1,n
            a(i,j)=a(i,j)*dum
        enddo
    endif
enddo
return
END
!
!----------------------------------------------------------------------
!
SUBROUTINE lubksb(a,n,np,indx,b)
IMPLICIT DOUBLE PRECISION(A-H,O-Z)
DIMENSION indx(n),a(np,np),b(n)
!------------------------------------------------------------------------
!Solves the set of n linear equations AX = B. Hereais input, not as the
!matrix A but rather as its LU decomposition, determined by the routine
!ludcmp. indx is input as the permutation vector returned by ludcmp.
!b(1:n) is input as the right-hand side vector B, and returns with the
!solution vector X. a, n, np, and indx are not modi^Led by this routine
!and can be left in place for successive calls with di^Kerent right-hand
!sides b. This routine takes into account the possibility that b will
!begin with many zero elements, so it is e^Ncient for use in matrix inversion.
!--------------------------------------------------------------------------
ii=0
do i=1,n
    ll=indx(i)
    sum=b(ll)
    b(ll)=b(i)
    if (ii.ne.0)then
        do j=ii,i-1
            sum=sum-a(i,j)*b(j)
        enddo
    else if (sum.ne.0.) then
        ii=i
    endif
    b(i)=sum
enddo
do i=n,1,-1
    sum=b(i)
    do j=i+1,n
        sum=sum-a(i,j)*b(j)
    enddo
    b(i)=sum/a(i,i)
enddo
return
END
