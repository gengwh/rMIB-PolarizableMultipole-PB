
!------------------------------------------------------------------------------
Function cappap(x,y,z)
use pbedata
Implicit Double Precision(A-H, O-Z)
cappap=-cappa**2*eps1
Return
End

Function cappan(x,y,z)
Implicit Double Precision(A-H, O-Z)
cappan=0.d0
Return
End

Function cappaj(x,y,z)
Implicit Double Precision(A-H, O-Z)
cappaj=cappap(x,y,z)-cappan(x,y,z)
Return
End

!------------------------------------------------------------------------------
Function betap(xx,yy,zz)
use comdata
use pbedata
use molecule
Implicit Double Precision(A-H, O-Z)
Betap=eps1
Return
End

Function betan(xx,yy,zz)
use comdata
use pbedata
use molecule
Implicit Double Precision(A-H, O-Z)
Betan=eps0
Return
End

Function beta(i,j,k)
Use pbedata
Implicit Double Precision(A-H, O-Z)
if (io(i,j,k)<=0) then
	beta=betan(x(i),y(j),z(k))
else
	beta=betap(x(i),y(j),z(k))
endif
End

Function betaj(x,y,z)
Implicit Double Precision(A-H, O-Z)
Betaj=Betap(x,y,z)-Betan(x,y,z)
Return
End

!Exact Solution for 1 charge in a ball
!--------------------------------------------------------------------------------------------
function upi(i,j,k)
use comdata
use pbedata
use molecule
implicit double precision(a-h,o-z)
real*8 pxyz(3)
pxyz=(/x(i),y(j),z(k)/)
call kirk_potential(1,pxyz,ptl,pmt,psv)
upi=ptl
end
!--------------------------------------------------------------------------------------------
function uni(i,j,k)
use comdata
use pbedata
use molecule
implicit double precision(a-h,o-z)
real*8 pxyz(3)
pxyz=(/x(i),y(j),z(k)/)
call kirk_potential(-1,pxyz,ptl,pmt,psv)
uni=pmt+psv
end

!-----------------------------------------------------------------------------
function bdcond(xx,yy,zz)
use comdata
use pbedata
use molecule

implicit double precision(a-h, o-z)
real*8 cp, aa
real*8 oneoverr_grad(3)	! 1st gradient 
real*8 oneoverr_grad_x(3),oneoverr_grad_y(3),oneoverr_grad_z(3) ! 2nd gradient

bdcond=0.d0

if (inl==1) then
    cp=sqrt(cappa_nl/eps1)
else
    cp=cappa
endif

aa = 1.d0

! just use debye-huckel at charge position
do ichg=1,natm
	dxx=xx-atmpos(1,ichg)
	dyy=yy-atmpos(2,ichg)
	dzz=zz-atmpos(3,ichg)
!     radius=sqrt((xx-atmpos(1,ichg))**2+(yy-atmpos(2,ichg))**2+(zz-atmpos(3,ichg))**2)
	radius=sqrt(dxx**2+dyy**2+dzz**2) ! center
! 	fardist = sqrt(atmpos(1,ichg)**2+atmpos(2,ichg)**2+atmpos(3,ichg)**2)

    if (ibd==2) then
!     	print *, "bdcond is ibd==2"
        bdcond=bdcond+atmchr(ichg)/radius
        if (ipm==1) then
			! dipole
			tmp_d=radius**(-3.d0)
			oneoverr_x=-dxx*tmp_d
			oneoverr_y=-dyy*tmp_d
			oneoverr_z=-dzz*tmp_d
        	oneoverr_grad=(/oneoverr_x,oneoverr_y,oneoverr_z/)
        	bdcond = bdcond - 3.d0*eps1/(2.d0*eps1+1.d0) * dot_product(di_mom(:,ichg),oneoverr_grad)

			! quadrupole
			tmp_q=radius**(-5.d0)	
			tmp_q1=3.d0*tmp_q
			oneoverr_xx=(dxx**2)*tmp_q1-tmp_d
			oneoverr_yy=(dyy**2)*tmp_q1-tmp_d
			oneoverr_zz=(dzz**2)*tmp_q1-tmp_d
			oneoverr_xy=dxx*dyy*tmp_q1
			oneoverr_xz=dxx*dzz*tmp_q1
			oneoverr_yz=dyy*dzz*tmp_q1

			oneoverr_grad_x=(/oneoverr_xx,oneoverr_xy,oneoverr_xz/)
			oneoverr_grad_y=(/oneoverr_xy,oneoverr_yy,oneoverr_yz/)
			oneoverr_grad_z=(/oneoverr_xz,oneoverr_yz,oneoverr_zz/)

			bdcond = bdcond + 5.d0*eps1/(3.d0*eps1+2.d0) * (dot_product(quad_mom(:,1,ichg),oneoverr_grad_x) &
			+dot_product(quad_mom(:,2,ichg),oneoverr_grad_y)+dot_product(quad_mom(:,3,ichg),oneoverr_grad_z))/3.d0
        endif

    elseif (ibd==3) then
        bdcond=bdcond+atmchr(ichg)*exp(-radius*cp)/radius
		if (ipm==1) then
			bdcond=bdcond+atmchr(ichg)/radius*exp(cp*(aa-radius))/(1.d0+cp*aa)
! 			bdcond=bdcond+atmchr(ichg)/radius*exp(cp*(atmrad(ichg)-radius))/(1.d0+cp*atmrad(ichg))
! 			! dipole
! 			tmp_d=radius**(-3.d0)
! 			oneoverr_x=-dxx*tmp_d*exp(-radius*cp)*(cp*radius+1.d0)
! 			oneoverr_y=-dyy*tmp_d*exp(-radius*cp)*(cp*radius+1.d0)
! 			oneoverr_z=-dzz*tmp_d*exp(-radius*cp)*(cp*radius+1.d0)
!         	oneoverr_grad=(/oneoverr_x,oneoverr_y,oneoverr_z/)
!         	bdcond = bdcond - dot_product(di_mom(:,ichg),oneoverr_grad)

! 			! quadrupole
! 			tmp_q=radius**(-5.d0)
! 			tmp_q1=3.d0*tmp_q
! 			tmp_q2=radius**(-4.d0)
! 			oneoverr_xx=exp(-radius*cp) *((dxx**2)*(cp**2)*tmp_d + (2.d0*(dxx**2)-(dyy**2)-(dzz**2))*(cp*tmp_q2+tmp_q))
! 			oneoverr_yy=exp(-radius*cp) *((dyy**2)*(cp**2)*tmp_d + (-(dxx**2)+2.d0*(dyy**2)-(dzz**2))*(cp*tmp_q2+tmp_q))
! 			oneoverr_zz=exp(-radius*cp) *((dzz**2)*(cp**2)*tmp_d + (-(dxx**2)-(dyy**2)+2.d0*(dzz**2))*(cp*tmp_q2+tmp_q))

! 			oneoverr_xy=dxx*dyy*exp(-radius*cp) * (cp*cp*tmp_d + 3.d0*cp*tmp_q2 + 3.d0*tmp_q)
! 			oneoverr_xz=dxx*dzz*exp(-radius*cp) * (cp*cp*tmp_d + 3.d0*cp*tmp_q2 + 3.d0*tmp_q)
! 			oneoverr_yz=dyy*dzz*exp(-radius*cp) * (cp*cp*tmp_d + 3.d0*cp*tmp_q2 + 3.d0*tmp_q)

! 			oneoverr_grad_x=(/oneoverr_xx,oneoverr_xy,oneoverr_xz/)
! 			oneoverr_grad_y=(/oneoverr_xy,oneoverr_yy,oneoverr_yz/)
! 			oneoverr_grad_z=(/oneoverr_xz,oneoverr_yz,oneoverr_zz/)
! 			bdcond = bdcond + 1.d0/3.d0*(dot_product(quad_mom(:,1,ichg),oneoverr_grad_x) &
! 			+dot_product(quad_mom(:,2,ichg),oneoverr_grad_y)+dot_product(quad_mom(:,3,ichg),oneoverr_grad_z))	

			! dipole
			tmp_d=radius**(-3.d0)
			oneoverr_x=-dxx*tmp_d
			oneoverr_y=-dyy*tmp_d
			oneoverr_z=-dzz*tmp_d
        	oneoverr_grad=(/oneoverr_x,oneoverr_y,oneoverr_z/)
        	bdcond = bdcond -dot_product(di_mom(:,ichg),oneoverr_grad)*3.d0*eps1*exp(cp*(aa-radius))*(1.d0+cp*radius) &
        	/(1.d0+cp*aa+eps1*(2.d0+2.d0*cp*aa+(cp*aa)**2))
!         	bdcond = bdcond -dot_product(di_mom(:,ichg),oneoverr_grad)*3.d0*eps1*exp(cp*(atmrad(ichg)-radius)) &
!         	*(1.d0+cp*radius) &
!         	/(1.d0+cp*atmrad(ichg)+eps1*(2.d0+2.d0*cp*atmrad(ichg)+(cp*atmrad(ichg))**2))

			! quadrupole
			tmp_q=radius**(-5.d0)	
			tmp_q1=3.d0*tmp_q
			oneoverr_xx=(dxx**2)*tmp_q1-tmp_d
			oneoverr_yy=(dyy**2)*tmp_q1-tmp_d
			oneoverr_zz=(dzz**2)*tmp_q1-tmp_d
			oneoverr_xy=dxx*dyy*tmp_q1
			oneoverr_xz=dxx*dzz*tmp_q1
			oneoverr_yz=dyy*dzz*tmp_q1

			oneoverr_grad_x=(/oneoverr_xx,oneoverr_xy,oneoverr_xz/)
			oneoverr_grad_y=(/oneoverr_xy,oneoverr_yy,oneoverr_yz/)
			oneoverr_grad_z=(/oneoverr_xz,oneoverr_yz,oneoverr_zz/)

			bdcond = bdcond + 5.d0*eps1*exp(cp*(aa-radius))*(3.d0+3.d0*cp*radius+(cp*raidus)**2) &
			/(2.d0*(3.d0+3.d0*cp*aa+(cp*aa)**2+eps1*(9.d0+9.d0*cp*aa+4.d0*(cp*aa)**2+(cp*aa)**3))) &
			* (dot_product(quad_mom(:,1,ichg),oneoverr_grad_x) &
			+dot_product(quad_mom(:,2,ichg),oneoverr_grad_y)+dot_product(quad_mom(:,3,ichg),oneoverr_grad_z))/3.d0
!  			bdcond = bdcond + 5.d0*eps1*exp(cp*(atmrad(ichg)-radius))*(3.d0+3.d0*cp*radius+(cp*raidus)**2) &
! 			/(2.d0*(3.d0+3.d0*cp*atmrad(ichg)+(cp*atmrad(ichg))**2+eps1*(9.d0+9.d0*cp*atmrad(ichg)+&
! 				4.d0*(cp*atmrad(ichg))**2+(cp*atmrad(ichg))**3))) &
! 			* (dot_product(quad_mom(:,1,ichg),oneoverr_grad_x) &
! 			+dot_product(quad_mom(:,2,ichg),oneoverr_grad_y)+dot_product(quad_mom(:,3,ichg),oneoverr_grad_z))/3.d0
 

		endif

    endif
enddo

! ! just use debye-huckel at charge position
! do ichg=1,natm
!     radius=sqrt((xx-atmpos(1,ichg))**2+(yy-atmpos(2,ichg))**2+(zz-atmpos(3,ichg))**2)
!     if (ibd==2) then
!         bdcond=bdcond+atmchr(ichg)/radius
!     elseif (ibd==3) then
!         bdcond=bdcond+atmchr(ichg)*exp(-radius*cp)/radius
!     endif
! enddo


bdcond=bdcond/eps1
End


!--------------------------------------------------------------------------------
subroutine chargedist(iatm)
use pbedata
use comdata
use molecule
implicit double precision(a-h,o-z)
dimension iloc(3),charge(8),loc_q(8,3),corlocq(8,3)
deltax=dcel;	deltay=dcel;	deltaz=dcel

charge=0.d0
iloc=1

x_q=atmpos(1,iatm)
y_q=atmpos(2,iatm)
z_q=atmpos(3,iatm)
q_q=atmchr(iatm)

i_q=inverx(x_q) ;	j_q=invery(y_q) ;	k_q=inverz(z_q)
x0=x(i_q)   ;		y0=y(j_q)   ;		z0=z(k_q)
xd1=x_q-x0  ;		yd1=y_q-y0  ;		zd1=z_q-z0

if(xd1==0.d0) iloc(1)=0
if(yd1==0.d0) iloc(2)=0
if(zd1==0.d0) iloc(3)=0

do i=0,1
	do j=0,1
		do k=0,1
			loc_q(k*4+j*2+i+1,:)=(/i_q+i,j_q+j,k_q+k/)
			corlocq(k*4+j*2+i+1,:)=(/x(i_q+i),y(j_q+j),z(k_q+k)/)
		enddo
	enddo
enddo

if(sum(iloc(:))==0)then
    charge(1)=1.d0
elseif(sum(iloc(:))==3)then
	do i=0,1
		do j=0,1
			do k=0,1
				xd=i*deltax-xd1
				yd=j*deltay-yd1
				zd=k*deltaz-zd1
				charge(4*k+2*j+i+1)=1.d0/abs(xd*yd*zd)
			enddo
		enddo
	enddo
elseif(sum(iloc(:))==2)then
    if(iloc(1)==0)then
		do j=0,1
			do k=0,1
				yd=j*deltay-yd1
				zd=k*deltaz-zd1
				charge(j*2+4*k+1)=1.d0/abs(yd*zd)
			enddo
		enddo
	elseif(iloc(2)==0)then
	    do i=0,1
		    do k=0,1
			    xd=i*deltax-xd1
				zd=k*deltaz-zd1
				charge(4*k+i+1)=1.d0/abs(xd*zd)
			enddo
		enddo
	elseif(iloc(3)==0)then
	    do i=0,1
		    do j=0,1
			    xd=i*deltax-xd1
				yd=j*deltay-yd1
				charge(i+2*j+1)=1.d0/abs(xd*yd)
			enddo
		enddo
	endif
elseif(sum(iloc(:))==1)then
    if(iloc(1)==1)then
	    charge(1)=1.d0/xd1
		charge(2)=1.d0/(deltax-xd1)
	elseif(iloc(2)==1)then
	    charge(1)=1.d0/yd1
		charge(3)=1.d0/(deltay-yd1)
    elseif(iloc(3)==1)then
	    charge(1)=1.d0/zd1
		charge(5)=1.d0/(deltaz-zd1)
	endif
endif

asum=sum(charge(:))
charge=q_q*charge/asum

loc_qt(iatm,:,:)=loc_q(:,:)
corlocqt(iatm,:,:)=corlocq(:,:)
charget(iatm,:)=charge(:)

return
end



!-------------------------------------------
!
function inverx(x)
use comdata
implicit double precision(a-h,o-z)

inverx=int((x-xleft)/dcel)+1

return
end
!
!-------------------------------------------
!
function invery(y)
use comdata
implicit double precision(a-h,o-z)

invery=int((y-yleft)/dcel)+1

return
end
!
!-------------------------------------------
!
function inverz(z)
use comdata
implicit double precision(a-h,o-z)

inverz=int((z-zleft)/dcel)+1

return
end
