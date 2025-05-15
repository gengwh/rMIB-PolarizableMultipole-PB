 
Function fn(x,y,z)
Implicit Double Precision(A-H, O-Z)
r2=x**2+y**2+z**2
fn=0.d0
Return
End
!------------------------------------------------------------------------------
Function fp(xx,yy,zz)
use pbedata
Implicit Double Precision(A-H, O-Z)
r2=xx**2+yy**2+zz**2
r=sqrt(r2)
fp=0.d0

if (inl==1) then
    fp=-cappa_nl*sinh(1.d0/eps/r)
endif


Return
End
!------------------------------------------------------------------------------
Function fj(x,y,z)
Implicit Double Precision(A-H, O-Z)
fj=fp(x,y,z)-fn(x,y,z)
Return
End

!---------------------------------------------------------------------------
function Up(xx,yy,zz) !outside
use comdata
use molecule
use pbedata
Implicit Double Precision(A-H, O-Z)
r2=xx**2+yy**2+zz**2; r=sqrt(r2)
up = 0.d0
if (abs(cappa)>1.d-10 .and. inl==0) then !linear PB
    up=1.d0/eps1/(1.d0+cappa*rds)/r*exp(-cappa*(r-rds))
else !Poisson or NLPB
!     do i=1,nchr
!         dxx=(xx-atmpos(1,i))
!         dyy=(yy-atmpos(2,i))
!         dzz=(zz-atmpos(3,i))
!         r2=dxx**2+dyy**2+dzz**2
!         r=sqrt(r2)

!         up=1.d0/eps1/r
! !         up=up + 3.d0/(2.d0*eps1+1.d0)*(dxx*di_mom(1,1)+dyy*di_mom(2,1)+dzz*di_mom(3,1))/(r**3)   !yang dipole
! !         up=up + 5.d0/(3.d0*eps1+2.d0)/(r**5)*(quad_mom(1,1,1)*dxx*dxx &
! !             +2.d0*quad_mom(2,1,1)*dxx*dyy+2.d0*quad_mom(3,1,1)*dxx*dzz+quad_mom(2,2,1)*dyy*dyy &
! !             +2.d0*quad_mom(3,2,1)*dyy*dzz+quad_mom(3,3,1)*dzz*dzz) !yang quadrupole *3/3.
!     enddo

    ! TABLE III IN SCHINIEDER??? YANG
    if (ipm==0) then
        up=1.d0/eps1/r !yang monopole charge=1, one charge
    elseif (ipm==1) then
        up=1.d0/eps1/r
        up=up + 3.d0/(2.d0*eps1+1.d0)*(xx*di_mom(1,1)+yy*di_mom(2,1)+zz*di_mom(3,1))/(r**3)   !yang dipole
        up=up + 5.d0/(3.d0*eps1+2.d0)/(r**5)*(quad_mom(1,1,1)*xx*xx &
            +2.d0*quad_mom(2,1,1)*xx*yy+2.d0*quad_mom(3,1,1)*xx*zz+quad_mom(2,2,1)*yy*yy &
            +2.d0*quad_mom(3,2,1)*yy*zz+quad_mom(3,3,1)*zz*zz) !yang quadrupole *3/3.

    endif

endif

Return
End

!------------------------------------------------------------------
function Un(xx,yy,zz)
use comdata
use molecule
use pbedata
Implicit Double Precision(A-H, O-Z)
r2=xx**2+yy**2+zz**2
r=sqrt(r2)
un = 0.d0

if (abs(cappa)>1.d-10 .and. inl==0) then !linear PB
    un=(1.d0/eps1/(1.d0+cappa*rds)/rds-1.d0/rds/eps0+1.d0/r/eps0)
else !Poisson or NLPB
!     do i=1,nchr
!         dxx=(xx-atmpos(1,i))
!         dyy=(yy-atmpos(2,i))
!         dzz=(zz-atmpos(3,i))
!         r2=dxx**2+dyy**2+dzz**2
!         r=sqrt(r2)

!         un=1.d0/r/eps0 -(1.d0-1.d0/eps)/rds/eps0 !yang monopole
! !         dtmp=(dxx*di_mom(1,i)+dyy*di_mom(2,i)+dzz*di_mom(3,i)) !dipole 
! !         un=un+dtmp/(r**3)-2.d0*(eps1-1.d0)/(2.d0*eps1+1.d0)/(rds**3)*dtmp !yang
! !         qtmp=(quad_mom(1,1,i)*dxx*dxx+2.d0*quad_mom(2,1,i)*dxx*dyy+2.d0*quad_mom(3,1,i)*dxx*dzz &
! !             +quad_mom(2,2,i)*dyy*dyy+2.d0*quad_mom(3,2,i)*dyy*dzz+quad_mom(3,3,i)*dzz*dzz)  !yang 
! !         un=un+1.d0/(r**5)*qtmp-3.d0*(eps1-1.d0)/(eps1*3.d0+2.d0)/(rds**5)*qtmp !quadrupole *3./3
!     enddo

    if (ipm==0) then
        un=1.d0/r/eps0 -(1.d0-1.d0/eps)/rds/eps0 !yang monopole
    elseif (ipm==1) then 
        un=1.d0/r/eps0 -(1.d0-1.d0/eps)/rds/eps0 !yang monopole
        dtmp=(xx*di_mom(1,1)+yy*di_mom(2,1)+zz*di_mom(3,1)) !dipole 
        un=un+dtmp/(r**3)-2.d0*(eps1-1.d0)/(2.d0*eps1+1.d0)/(rds**3)*dtmp !yang
        qtmp=(quad_mom(1,1,1)*xx*xx+2.d0*quad_mom(2,1,1)*xx*yy+2.d0*quad_mom(3,1,1)*xx*zz &
            +quad_mom(2,2,1)*yy*yy+2.d0*quad_mom(3,2,1)*yy*zz+quad_mom(3,3,1)*zz*zz)  !yang 
        un=un+1.d0/(r**5)*qtmp-3.d0*(eps1-1.d0)/(eps1*3.d0+2.d0)/(rds**5)*qtmp !quadrupole *3/3.d0

    endif


endif

Return
End

!----------------------------------------------------------------------------
function Uj(xx,yy,zz)
use pbedata
Implicit Double Precision(A-H, O-Z)
if (ibd==2 .or. ibd==3) then
	Uj=0.d0
else
	Uj=Up(xx,yy,zz)-Un(xx,yy,zz)
!     print *, "un at uj is ", Un(xx,yy,zz)
endif
Return
End



!------------------------------------------------------------------------------
function Udxp(xx,yy,zz)
use comdata
use pbedata
Implicit Double Precision(A-H, O-Z)
r2=xx**2+yy**2+zz**2
r=sqrt(r2)


if (abs(cappa)>1.d-10 .and. inl==0) then !linear PB
    Udxp=-exp(-cappa*(r-rds))*xx*(1.d0+cappa*r)/eps1/(1+cappa*rds)/r**3
else !Poisson or NLPB
    Udxp=-1.d0/eps1*xx/r**3
endif

Return
End

function Udxn(xx,yy,zz)
use comdata
use pbedata
Implicit Double Precision(A-H, O-Z)
r2=xx**2+yy**2+zz**2
r=sqrt(r2)
Udxn=-1.d0*xx/r**3/eps0 ! The same for Poisson and PB
Return
End

function Udxj(x,y,z)
Implicit Double Precision(A-H, O-Z)
Udxj=Udxp(x,y,z)-Udxn(x,y,z)
Return
End


function Udyp(xx,yy,zz)
use comdata
use pbedata
Implicit Double Precision(A-H, O-Z)
r2=xx**2+yy**2+zz**2
r=sqrt(r2)


if (abs(cappa)>1.d-10 .and. inl==0) then !linear PB
    Udyp=-exp(-cappa*(r-rds))*yy*(1.d0+cappa*r)/eps1/(1+cappa*rds)/r**3
else !Poisson or NLPB
    Udyp=-1.d0/eps1*yy/r**3
endif

Return
End

function Udyn(xx,yy,zz)
use comdata
use pbedata
Implicit Double Precision(A-H, O-Z)
r2=xx**2+yy**2+zz**2
r=sqrt(r2)
Udyn=-1.d0*yy/r**3/eps0 ! The same for Poisson and PB
Return
End

function Udyj(x,y,z)
Implicit Double Precision(A-H, O-Z)
Udyj=Udyp(x,y,z)-Udyn(x,y,z)
Return
End


function Udzp(xx,yy,zz)
use comdata
use pbedata
Implicit Double Precision(A-H, O-Z)
r2=xx**2+yy**2+zz**2
r=sqrt(r2)


if (abs(cappa)>1.d-10 .and. inl==0) then !linear PB
    Udzp=-exp(-cappa*(r-rds))*zz*(1.d0+cappa*r)/eps1/(1+cappa*rds)/r**3
else !Poisson or NLPB
    Udzp=-1.d0/eps1*zz/r**3
endif


Return
End

function Udzn(xx,yy,zz)
use comdata
use pbedata
Implicit Double Precision(A-H, O-Z)
r2=xx**2+yy**2+zz**2
r=sqrt(r2)
Udzn=-1.d0*zz/r**3/eps0 ! The same for Poisson and PB
Return
End

function Udzj(x,y,z)
Implicit Double Precision(A-H, O-Z)
Udzj=Udzp(x,y,z)-Udzn(x,y,z)
Return
End
