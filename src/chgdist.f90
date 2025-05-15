!subroutine chgdist
!use comdata
! collection of subroutines related to charge distribution including:
! 1: Green's function and its derivative
! 2: subroutines for computing solvation energy
! yang: using multipole 
function phi_star(xx,yy,zz)
use comdata
use molecule
use pbedata

implicit double precision(a-h,o-z)
real*8 oneoverr_grad(3)	! 1st gradient 
real*8 oneoverr_grad_x(3),oneoverr_grad_y(3),oneoverr_grad_z(3) ! 2nd gradient

C=1.d0/eps0
phi_star=0.d0
if (ipm==1) then 
	phi_star_m=0.d0
	phi_star_d=0.d0
	phi_star_q=0.d0

	do i=1,nchr
		dxx=xx-chrpos(1,i)
		dyy=yy-chrpos(2,i)
		dzz=zz-chrpos(3,i)
		dist_sqr=dxx**2+dyy**2+dzz**2
		if (abs(dxx)<1.d-10 .and. abs(dyy)<1.d-10 .and. abs(dzz)<1.d-10) then 
			phi_star_m=phi_star_m+0.d0
			phi_star_d=phi_star_d+0.d0
			phi_star_q=phi_star_q+0.d0
		else
			! monopole
			phi_star_m=phi_star_m+atmchr(i)/sqrt(dist_sqr)

			! dipole
			tmp_d=dist_sqr**(-1.5d0)
			oneoverr_x=-dxx*tmp_d
			oneoverr_y=-dyy*tmp_d
			oneoverr_z=-dzz*tmp_d
			oneoverr_grad=(/oneoverr_x,oneoverr_y,oneoverr_z/)
			phi_star_d=phi_star_d-dot_product(di_mom(:,i),oneoverr_grad)

			! quadrupole
			tmp_q=dist_sqr**(-2.5d0)	
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

			! yang March 2024, add 1/3 constant to quadrupole component
			! the 3/2 constant count in the quadrupole moment if traceless
			phi_star_q=phi_star_q + (dot_product(quad_mom(:,1,i),oneoverr_grad_x) &
			+dot_product(quad_mom(:,2,i),oneoverr_grad_y)+dot_product(quad_mom(:,3,i),oneoverr_grad_z))/3.d0
		endif
	enddo
	phi_star=(phi_star_m+phi_star_d+phi_star_q)*C
! 	phi_star=(phi_star_d)*C

elseif(ipm==0) then
	do i=1,nchr
		if (abs(xx-chrpos(1,i))<1.d-10 .and. abs(yy-chrpos(2,i))<1.d-10 .and. abs(zz-chrpos(3,i))<1.d-10) then 
			phi_star=phi_star+0.d0
		else
			phi_star=phi_star+atmchr(i)/sqrt((xx-chrpos(1,i))**2+(yy-chrpos(2,i))**2+(zz-chrpos(3,i))**2)
		endif
	enddo
	phi_star=phi_star*C
endif

End


subroutine phi_star_grad_mul(xx,yy,zz,phi_star_grad_d,phi_star_grad_q)
use comdata
use molecule
use pbedata
implicit double precision(a-h,o-z)
real*8 oneoverr_grad_x(3),oneoverr_grad_y(3),oneoverr_grad_z(3),phi_star_grad_d(3)
real*8 oneoverr_grad_x_x(3),oneoverr_grad_x_y(3),oneoverr_grad_x_z(3),oneoverr_grad_y_y(3)
real*8 oneoverr_grad_y_z(3),oneoverr_grad_z_z(3),phi_star_grad_q(3)
C=1.d0/eps0
phi_star_x_d=0.d0 ! dipole
phi_star_y_d=0.d0
phi_star_z_d=0.d0

phi_star_x_q=0.d0 ! quadrupole
phi_star_y_q=0.d0
phi_star_z_q=0.d0

do i=1,nchr
	dxx=(xx-chrpos(1,i))
	dyy=(yy-chrpos(2,i))
	dzz=(zz-chrpos(3,i))

	if (abs(dxx)<1.d-10 .and. abs(dyy)<1.d-10 .and. abs(dzz)<1.d-10) then 
		phi_star_x_d=phi_star_x_d+0.d0
		phi_star_y_d=phi_star_y_d+0.d0
		phi_star_z_d=phi_star_z_d+0.d0
		phi_star_x_q=phi_star_x_q+0.d0
		phi_star_y_q=phi_star_y_q+0.d0
		phi_star_z_q=phi_star_z_q+0.d0
	else
		!dipole
		dist_sqr=dxx**2+dyy**2+dzz**2
		tmp_d=dist_sqr**(-1.5d0)
		tmp_q=dist_sqr**(-2.5d0)
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

		phi_star_x_d=phi_star_x_d-dot_product(di_mom(:,i),oneoverr_grad_x)
		phi_star_y_d=phi_star_y_d-dot_product(di_mom(:,i),oneoverr_grad_y)
		phi_star_z_d=phi_star_z_d-dot_product(di_mom(:,i),oneoverr_grad_z)

		!quadrupole
		tmp_q2=9.d0*tmp_q
		tmp_q3=15.d0*(dist_sqr**(-3.5d0))
		oneoverr_xxx=dxx*tmp_q2-tmp_q3*(dxx**3)
		oneoverr_yyy=dyy*tmp_q2-tmp_q3*(dyy**3)
		oneoverr_zzz=dzz*tmp_q2-tmp_q3*(dzz**3)
	
		oneoverr_xxy=dyy*tmp_q1-tmp_q3*(dxx**2)*dyy
		oneoverr_xxz=dzz*tmp_q1-tmp_q3*(dxx**2)*dzz
	
		oneoverr_xyy=dxx*tmp_q1-tmp_q3*(dyy**2)*dxx
		oneoverr_xzz=dxx*tmp_q1-tmp_q3*(dzz**2)*dxx
		oneoverr_xyz=-tmp_q3*dxx*dyy*dzz
	
		oneoverr_yyz=dzz*tmp_q1-tmp_q3*(dyy**2)*dzz
		oneoverr_yzz=dyy*tmp_q1-tmp_q3*(dzz**2)*dyy


		oneoverr_grad_x_x=(/oneoverr_xxx,oneoverr_xxy,oneoverr_xxz/)
		oneoverr_grad_x_y=(/oneoverr_xxy,oneoverr_xyy,oneoverr_xyz/)
		oneoverr_grad_x_z=(/oneoverr_xxz,oneoverr_xyz,oneoverr_xzz/)
		oneoverr_grad_y_y=(/oneoverr_xyy,oneoverr_yyy,oneoverr_yyz/)
		oneoverr_grad_y_z=(/oneoverr_xyz,oneoverr_yyz,oneoverr_yzz/)
		oneoverr_grad_z_z=(/oneoverr_xzz,oneoverr_yzz,oneoverr_zzz/)

	
		phi_star_x_q=phi_star_x_q+dot_product(oneoverr_grad_x_x,quad_mom(:,1,i))+dot_product(oneoverr_grad_x_y,quad_mom(:,2,i)) &
		+dot_product(oneoverr_grad_x_z,quad_mom(:,3,i))
		phi_star_y_q=phi_star_y_q+dot_product(oneoverr_grad_x_y,quad_mom(:,1,i))+dot_product(oneoverr_grad_y_y,quad_mom(:,2,i)) &
		+dot_product(oneoverr_grad_y_z,quad_mom(:,3,i))
		phi_star_z_q=phi_star_z_q+dot_product(oneoverr_grad_x_z,quad_mom(:,1,i))+dot_product(oneoverr_grad_y_z,quad_mom(:,2,i)) &
		+dot_product(oneoverr_grad_z_z,quad_mom(:,3,i))

	endif
enddo

phi_star_grad_d=(/phi_star_x_d,phi_star_y_d,phi_star_z_d/)*C
phi_star_grad_q=(/phi_star_x_q,phi_star_y_q,phi_star_z_q/)*C/3.d0! yang march 2024

End


function phis_x(xx,yy,zz)
use comdata
use molecule
use pbedata
implicit double precision(a-h,o-z)
C=1.d0/eps0
phis_x=0.d0
do i=1,nchr
	phis_x=phis_x-atmchr(i)*((xx-chrpos(1,i))**2+(yy-chrpos(2,i))**2+(zz-chrpos(3,i))**2) & 
	**(-1.5d0)*((xx-chrpos(1,i)))
enddo
phis_x=phis_x*C
End

function phis_y(xx,yy,zz)
use comdata
use molecule
use pbedata
implicit double precision(a-h,o-z)
C=1.d0/eps0
phis_y=0.d0
do i=1,nchr
	phis_y=phis_y-atmchr(i)*((xx-chrpos(1,i))**2+(yy-chrpos(2,i))**2+(zz-chrpos(3,i))**2) & 
	**(-1.5d0)*((yy-chrpos(2,i)))
enddo
phis_y=phis_y*C
End

function phis_z(xx,yy,zz)
use comdata
use molecule
use pbedata
implicit double precision(a-h,o-z)
C=1.d0/eps0
phis_z=0.d0
do i=1,nchr
	phis_z=phis_z-atmchr(i)*((xx-chrpos(1,i))**2+(yy-chrpos(2,i))**2+(zz-chrpos(3,i))**2) & 
	**(-1.5d0)*((zz-chrpos(3,i)))
enddo
phis_z=phis_z*C
End

!---------------------------------------------------------------------
subroutine solenergy(engy)	! yang
use comdata
use molecule
use pbedata

implicit double precision(a-h,o-z)
real*8 se(2),diff(3),engy 		
integer icood(3)
para1=332.0716d0
soleng=0.d0; couleng=0.d0


if (ipm==1) then ! point multipole
	do i=1,nchr
! 		! numerical approaches:
	    soleng=soleng + atmchr(i)*phibar_intp(chrpos(1,i),chrpos(2,i),chrpos(3,i))	! yang:monopole
		soleng=soleng+dot_product(di_mom(:,i),ugrad(:,i))	! yang:dipole
		soleng=soleng + (dot_product(quad_mom(:,1,i),ugradgrad(:,1,i))+dot_product(quad_mom(:,2,i),ugradgrad(:,2,i)) &
		    +dot_product(quad_mom(:,3,i),ugradgrad(:,3,i)))/3.0	! yang 1/3 for quadrupole component for traceless

	    do j=i+1,nchr
	        diff=chrpos(:,i)-chrpos(:,j)
	        dist=sqrt(dot_product(diff,diff))
	        couleng=couleng+1/eps0/dist*atmchr(i)*atmchr(j)
	    enddo
	    
! 	    checking using analytical solution
! 	    soleng=soleng+atmchr(i)*atmchr(i)*(1/rds*(1/eps-1))	!atmchr(i)=atmchr(i)*atmchr(i)=1 yang
! 	    soleng=soleng+(-2.d0*(eps-1.d0)/(1.d0+2.d0*eps))/(rds**3)*dot_product(di_mom(:,i),di_mom(:,i)) !yang
! 		soleng=soleng+(-3.d0*(eps-1.d0)/(2.d0+3.d0*eps))/(rds**5)*sum(quad_mom(:,:,i)**2)*2.d0/3.d0 !schnieder

	enddo
elseif(ipm==0) then  ! point charge
	do i=1,nchr
	    soleng=soleng+atmchr(i)*phibar_intp(chrpos(1,i),chrpos(2,i),chrpos(3,i))
	    do j=i+1,nchr
	        diff=chrpos(:,i)-chrpos(:,j)
	        dist=sqrt(dot_product(diff,diff))
	        couleng=couleng+1/eps0/dist*atmchr(i)*atmchr(j)
	    enddo
	    !checking using analytical solution
	    !soleng=soleng+atmchr(i)*(1/rds*(1/eps-1))
	enddo
endif

! print*, "ugrad is ", ugrad(:,2)
soleng=0.5d0*soleng*para1; couleng=couleng*para1
print *,'solvation energy =: ',soleng, ' kcal/mol.'
print *,'-----free energy =: ',couleng+soleng,' kcal/mol.'
engy=soleng

end



!----------------------------------------------------------------------------------------
subroutine irr_interface
use comdata
use pbedata
use molecule
implicit double precision(a-h, o-z)
integer idx(3), ix1(3), ix2(3)

allocate(irrintf(6,nirr), rrintf(6,nirr), onsrf(6,nirr))
irrintf=0;	rrintf=0.d0;	onsrf=0.d0
xacc=1.d-10


do irr=1,nirr
	i=irrxyz(1,irr); j=irrxyz(2,irr); k=irrxyz(3,irr)
	if (io(i,j,k)==0) then		! For grids on interface 
		continue	
	elseif (irrpts(i,j,k)==1 .and. io(i,j,k)==-1) then	! Irregular points
		
		do ik=1,6
			call idx6(i,j,k,ik,idx)					! idx stores coodinates of the fictitious points outside
			if (io(idx(1),idx(2),idx(3))==1) then	! fictitious point is needed in this direction
				irrintf(ik,irr)=1
				ix1=(/min(i,idx(1)),min(j,idx(2)),min(k,idx(3))/)
				ix2=(/max(i,idx(1)),max(j,idx(2)),max(k,idx(3))/)

!###################################################################################################				
!This code gets the surface information from dabao's code				
				if (isf==1 .or. isf==2 .or. isf==4) then
					irr0=irr
					if (ix1(1)==ix2(1) .and. ix1(2)==ix2(2)) then				! change in z direction
						xs=x(i);	ys=y(j);	
						irr0=indirr(i,j,k-1)
						if (mcz(irr) == -1) then								! ++++(==*---->  ==:DZR
							zs=z(k-1)+dzl(irr0)
							onsrf(ik,irr)=zs-z(k)								! Interface position as center
						else
							!if (mcz(irr0) ==1 .and. ik == 1 ) ! bug: irr0 could be zero
							if (ik == 1) then
								if (irr0==0) then 
									print *,'type II interface error',i,j,k
									goto 10
								endif
								zs=z(k-1)+dzl(irr0)								! +++++(==*===)++++ Both sides			
							else
								zs=z(k)+dzl(irr)								! ----*==)++++>  ==:DZL
							endif
							10 continue
							onsrf(ik,irr)=zs-z(k)
						endif
						rrintf(ik,irr)=zs
					elseif (ix1(1)==ix2(1) .and. ix1(3)==ix2(3)) then			! change in y direction
						xs=x(i);	zs=z(k);
						irr0=indirr(i,j-1,k)
						if (mcy(irr) == -1) then								! ++++(==*---->  ==:DYR	
							ys=y(j-1)+dyl(irr0)
							onsrf(ik,irr)=ys-y(j)	
						else
							!if (mcy(irr0)==1 .and. ik==2) then
							if (ik==2) then
								if (irr0==0) then 
									print *,'type II interface error',i,j,k
									goto 20
								endif
								ys=y(j-1)+dyl(irr0)
							else
								ys=y(j)+dyl(irr)
							endif
							20 continue
							onsrf(ik,irr)=ys-y(j)
						endif
						rrintf(ik,irr)=ys
					else														! change in x direction
						ys=y(j);	zs=z(k);
						irr0=indirr(i-1,j,k)
						if (mcx(irr) == -1) then								! ++++(==*---->  ==:DXR
							xs=x(i-1)+dxl(irr0)
							onsrf(ik,irr)=xs-x(i)	
						else
							!if (mcx(irr0)==1 .and. ik==3) then
							if (ik==3) then
								if (irr0==0) then 
									print *,'type II interface error',i,j,k
									goto 30
								endif
								xs=x(i-1)+dxl(irr0)
							else
								xs=x(i)+dxl(irr)
							endif
							30 continue
							onsrf(ik,irr)=xs-x(i)
						endif
						rrintf(ik,irr)=xs
					endif
				endif
!###################################################################################################	

!***************************************************************************************************
				if (isf==0) then
					if (ix1(1)==ix2(1) .and. ix1(2)==ix2(2)) then				! change in z direction
						!^^^^^^^^^^^^^^^^^^^^^^^^^
						!if (z(ix1(3))>1.d-10) then 
						!	rrintf(ik,irr)=1.d0
						!else
						!	rrintf(ik,irr)=-1.d0
						!endif
						rrintf(ik,irr)=rtbis(z(ix1(3)),z(ix2(3)),xacc,x(i),y(j),0.d0,3)
					elseif (ix1(1)==ix2(1) .and. ix1(3)==ix2(3)) then			! change in y direction
						!^^^^^^^^^^^^^^^^^^^^^^^^^
						!if (y(ix1(2))>1.d-10) then 
						!	rrintf(ik,irr)=1.d0
						!else
						!	rrintf(ik,irr)=-1.d0
						!endif
						rrintf(ik,irr)=rtbis(y(ix1(2)),y(ix2(2)),xacc,x(i),0.d0,z(k),2)  
					else														! change in x direction
						!^^^^^^^^^^^^^^^^^^^^^^^^^
						!if (x(ix1(1))>1.d-10) then 
						!	rrintf(ik,irr)=1.d0
						!else
						!	rrintf(ik,irr)=-1.d0
						!endif
						rrintf(ik,irr)=rtbis(x(ix1(1)),x(ix2(1)),xacc,0.d0,y(j),z(k),1)	
					endif
				endif
!***************************************************************************************************
			endif
		enddo
	endif
enddo



end


!-----------------------------------------------------------------
subroutine coodinates(position,icood)
use comdata
use molecule
use pbedata
implicit double precision(a-h,o-z)
real*8 position(3)
integer icood(3)

do k=1,nz
	if (abs(position(3)-z(k))<1.d-10) then
		icood(3)=k
		goto 100
	endif
enddo
100 continue

do j=1,ny
	if (abs(position(2)-y(j))<1.d-10) then
		icood(2)=j
		goto 200
	endif
enddo
200 continue

do i=1,nx
	if (abs(position(1)-x(i))<1.d-10) then
		icood(1)=i
		goto 300
	endif
enddo
300 continue
End

!-----------------------------------------
function phibar_intp(xx,yy,zz)
use comdata
use molecule
use pbedata
implicit double precision(a-h,o-z)
real*8 xwt(-1:1),ywt(-1:1),zwt(-1:1),aa(3),wt(27),value(27),w,ax(nx),by(ny),cz(nz),bb(3),cc(3)
integer ixyz(1)

ixyz=minloc(abs(x-xx))
ix=ixyz(1)
ixyz=minloc(abs(y-yy))
iy=ixyz(1)
ixyz=minloc(abs(z-zz))
iz=ixyz(1)

aa=x(ix-1:ix+1)
call weights(xx,aa,2,2,0,xwt)
aa=y(iy-1:iy+1)
call weights(yy,aa,2,2,0,ywt)
aa=z(iz-1:iz+1)
call weights(zz,aa,2,2,0,zwt)

ijk=0
do k=-1,1
	do j=-1,1
		do i=-1,1
			ijk=ijk+1
			wt(ijk)=xwt(i)*ywt(j)*zwt(k)
			if (io(ix+i,iy+j,iz+k)>0) then
				if (ix+i>=4) then
					if (io(ix+i-1,iy+j,iz+k)<=0 .and. io(ix+i-2,iy+j,iz+k)<=0 .and. io(ix+i-3,iy+j,iz+k)<=0) then
						aa=x(ix+i-1:ix+i-3:-1)
						call weights(x(ix),aa,2,2,0,bb)
						Do ii=1,3
							cc(ii)=phi0(id3d(ix+i-ii,iy+j,iz+k))+u(id3d(ix+i-ii,iy+j,iz+k))
						enddo
						value(ijk)=dot_product(bb,cc)
						goto 100
					endif
				elseif(ix+i<=nx-3) then
					if(io(ix+i+1,iy+j,iz+k)<=0 .and. io(ix+i+2,iy+j,iz+k)<=0 .and. io(ix+i+3,iy+j,iz+k)<=0) then
						aa=x(ix+i+1:ix+i+3)
						call weights(x(ix),aa,2,2,0,bb)
						Do ii=1,3
							cc(ii)=phi0(id3d(ix+i+ii,iy+j,iz+k))+u(id3d(ix+i+ii,iy+j,iz+k))
						enddo
						value(ijk)=dot_product(bb,cc)
						goto 100
					endif
				elseif(iy+j>=4) then
					if(io(ix+i,iy+j-3,iz+k)<=0 .and. io(ix+i,iy+j-2,iz+k)<=0 .and. io(ix+i,iy+j-1,iz+k)<=0) then
						aa=y(iy+j-1:iy+j-3:-1)
						call weights(y(iy),aa,2,2,0,bb)
						Do ii=1,3
							cc(ii)=phi0(id3d(ix+i,iy+j-ii,iz+k))+u(id3d(ix+i,iy+j-ii,iz+k))
						enddo
						value(ijk)=dot_product(bb,cc)
						goto 100
					endif
				elseif(iy+j<=ny-3) then
					if(io(ix+i,iy+j+1,iz+k)<=0 .and. io(ix+i,iy+j+2,iz+k)<=0 .and. io(ix+i,iy+j+3,iz+k)<=0) then
						aa=y(iy+j+1:iy+j+3)
						call weights(y(iy),aa,2,2,0,bb)
						Do ii=1,3
							cc(ii)=phi0(id3d(ix+i,iy+j+ii,iz+k))+u(id3d(ix+i,iy+j+ii,iz+k))
						enddo
						value(ijk)=dot_product(bb,cc)
						goto 100
					endif
				elseif(iz+k>=4) then
					if(io(ix+i,iy+j,iz+k-3)<=0 .and. io(ix+i,iy+j,iz+k-2)<=0 .and. io(ix+i,iy+j,iz+k-1)<=0) then
						aa=z(iz+k-1:iz+k-3:-1)
						call weights(z(iz),aa,2,2,0,bb)
						Do ii=1,3
							cc(ii)=phi0(id3d(ix+i,iy+j,iz+k-ii))+u(id3d(ix+i,iy+j,iz+k-ii))
						enddo
						value(ijk)=dot_product(bb,cc)
						goto 100
					endif
				elseif(iz+k<=nz-3) then
					if(io(ix+i,iy+j,iz+k+1)<=0 .and. io(ix+i,iy+j,iz+k+2)<=0 .and. io(ix+i,iy+j,iz+k+3)<=0) then
						aa=z(iz+k+1:iz+k+3)
						call weights(z(iz),aa,2,2,0,bb)
						Do ii=1,3
							cc(ii)=phi0(id3d(ix+i,iy+j,iz+k+ii))+u(id3d(ix+i,iy+j,iz+k+ii))
						enddo
						value(ijk)=dot_product(bb,cc)
						goto 100
					endif
				else
					print *,i,j,k,'failed to interpolate'
				endif
			else
				value(ijk)=phi0(id3d(ix+i,iy+j,iz+k))+u(id3d(ix+i,iy+j,iz+k))
			endif
100 continue
		enddo
	enddo
enddo
phibar_intp=dot_product(wt,value)

End

! !--------------------------------------------------------------------------------------------

! ! Feb 2024 yang

! function test_function(xx,yy,zz)
! use comdata
! use molecule
! use pbedata
! implicit double precision(a-h,o-z)

! ! test_function=cos(xx)*sin(yy)*exp(zz)
! test_function = 1.d0/sqrt(xx**2+yy**2+zz**2)
! End
! !--------------------------------------------------------------------------------------------

! subroutine test_intp(xx,yy,zz,test_intp_grad0,test_intp_grad1,test_intp_grad2)
! use comdata
! use molecule
! use pbedata
! implicit double precision(a-h,o-z)
! real*8 xwt(-1:1,0:2),ywt(-1:1,0:2),zwt(-1:1,0:2),aa(3),wt(27),value(27,0:2),w,ax(nx),by(ny),cz(nz),bb(3,0:2),cc(3) !, xwtt(-1:1,0:1)
! integer ixyz(1)
! real*8 test_intp_grad0
! real*8 test_intp_grad1(3), test_intp_grad2(3,3),wtx(27),wty(27),wtz(27),wtxx(27),wtxy(27),wtxz(27)
! real*8 wtyy(27),wtyz(27),wtzz(27)
! real*8 mesh
! integer num
! num = 10000 ! 100, 1000, 10000
! mesh = 0.002d0 ! 0.2, 0.02,0.002
! print *, "mesh is ", mesh
! allocate(xt(num),yt(num),zt(num))
! ! xt=0.d0; yt=0.d0; zt=0.d0
! do i=1,num
! 	xt(i)=0.d0+(i-1)*mesh
! enddo
! do i=1,num
! 	yt(i)=0.d0+(i-1)*mesh
! enddo
! do i=1,num
! 	zt(i)=0.d0+(i-1)*mesh
! enddo
! print *, "xt is ", xt(1:3)
! gradx=0.d0
! grady=0.d0
! gradz=0.d0
! gradgradxx=0.d0
! gradgradxy=0.d0
! gradgradxz=0.d0
! gradgradyy=0.d0
! gradgradyz=0.d0
! gradgradzz=0.d0


! ixyz=minloc(abs(xt-xx))
! ix=ixyz(1)
! ixyz=minloc(abs(yt-yy))
! iy=ixyz(1)
! ixyz=minloc(abs(zt-zz))
! iz=ixyz(1)

! aa=xt(ix-1:ix+1)
! call weights(xx,aa,2,2,2,xwt)
! aa=yt(iy-1:iy+1)
! call weights(yy,aa,2,2,2,ywt)
! aa=zt(iz-1:iz+1)
! call weights(zz,aa,2,2,2,zwt)

! ijk=0
! do k=-1,1
! 	do j=-1,1
! 		do i=-1,1
! 			ijk=ijk+1
! 			wt(ijk)=xwt(i,0)*ywt(j,0)*zwt(k,0)

! 			wtx(ijk)=xwt(i,1)*ywt(j,0)*zwt(k,0)
! 			wty(ijk)=xwt(i,0)*ywt(j,1)*zwt(k,0)
! 			wtz(ijk)=xwt(i,0)*ywt(j,0)*zwt(k,1)

! 			wtxx(ijk)=xwt(i,2)*ywt(j,0)*zwt(k,0)
! 			wtxy(ijk)=xwt(i,1)*ywt(j,1)*zwt(k,0)
! 			wtxz(ijk)=xwt(i,1)*ywt(j,0)*zwt(k,1)

! 			wtyy(ijk)=xwt(i,0)*ywt(j,2)*zwt(k,0)
! 			wtyz(ijk)=xwt(i,0)*ywt(j,1)*zwt(k,1)
! 			wtzz(ijk)=xwt(i,0)*ywt(j,0)*zwt(k,2)

! 			value(ijk,0)=test_function(xt(ix+i),yt(iy+j),zt(iz+k))
! 		enddo
! 	enddo
! enddo
! test_intp_grad0=dot_product(wt,value(:,0))

! gradx=dot_product(wtx,value(:,0))
! grady=dot_product(wty,value(:,0))
! gradz=dot_product(wtz,value(:,0))

! test_intp_grad1 = (/gradx,grady,gradz/)

! gradgradxx=dot_product(wtxx,value(:,0))
! gradgradxy=dot_product(wtxy,value(:,0))
! gradgradxz=dot_product(wtxz,value(:,0))
! gradgradyy=dot_product(wtyy,value(:,0))
! gradgradyz=dot_product(wtyz,value(:,0))
! gradgradzz=dot_product(wtzz,value(:,0))
! ! phibar_intp_grad2 = (/gradgradxx,gradgradxy,gradgradxz,gradgradyy,gradgradyz,gradgradzz/)

! test_intp_grad2(:,1)=(/gradgradxx,gradgradxy,gradgradxz/)
! test_intp_grad2(:,2)=(/gradgradxy,gradgradyy,gradgradyz/)
! test_intp_grad2(:,3)=(/gradgradxz,gradgradyz,gradgradzz/)

! End subroutine test_intp


