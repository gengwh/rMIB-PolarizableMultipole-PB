!This program is updated based on the following contents
!1) MIBPB-I (Y. Zhou et al.) with interface treatment 
!2) MIBPB-II(S. Yu et al.) with the shape interface treatment
!3) Since rMIB has distinct advantage in singular charge treatment, we will remove charge treatment in MIBPB-I-II and the three components green function based charge treatment in MIBPB-III

program mibpb_newReg
use comdata
use molecule
use pbedata
use bicg

implicit double precision(a-h,o-z)
! real*8 test_intp_grad0, test_intp_grad1(3), test_intp_grad2(3,3)
real*8 engy
character(7) cdens
real*8 rslt, cappa2
common /pi/ pi
character(100) fhead

call cpu_time(cpu0)
cpu1=cpu0
pi=acos(-1.d0)

!??????????????????????????????????????????????????????????????
!PARAMETERS: now can be read from usrdata.in file

!PB equation 
!eps0=1.d0;            !the dielectric constant in molecule 
!eps1=80.d0;           !the dielectric constant in solvent
!bulk_strength=0.15d0  !ion_strength with units (M)$I=\sum\limits_{i=1}^nc_iz_i^2$

open(101,file="usrdata.in")
READ(101,*,IOSTAT = MEOF) fhead, fname
READ(101,*,IOSTAT = MEOF) fhead, dcel 
READ(101,*,IOSTAT = MEOF) fhead, den 
READ(101,*,IOSTAT = MEOF) fhead, eps0 
READ(101,*,IOSTAT = MEOF) fhead, eps1 
READ(101,*,IOSTAT = MEOF) fhead, bulk_strength 
READ(101,*,IOSTAT = MEOF) fhead, icg 
READ(101,*,IOSTAT = MEOF) fhead, isf 
READ(101,*,IOSTAT = MEOF) fhead, ibd 
READ(101,*,IOSTAT = MEOF) fhead, inl 
READ(101,*,IOSTAT = MEOF) fhead, ipm
close(101)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!icg,isf,ibd gives the different handing of charge, interface and boundary
!icg=2;	0: no charge;		2: green function 
!isf=0;	0: exact interface	1: msms interface;	4: eses interface		
!ibd=0;	0: exact solution	1: kirkwood solution;	2: 1/r decay	3:exp decay
!inl=1; 0: linear PB; 		1: nonlinear PB 	
!ipm=1; 0: point charge     1: point multipole (monopole+dipole+quadrupole)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!reading the path if batch calculation
pathname='test_proteins/'
lenpath = len(pathname)
do while (pathname(lenpath:lenpath) .eq. ' ')
	lenpath = lenpath - 1
enddo  

!set up parameters: radius, dielectric constants, screening strength .../
prds=1.4d0
!eps0=1.d0	
!eps1=80.d0
!cappa=1.d0
cappa2=8.430325455*bulk_strength/eps1     !cappa2 in 300K
cappa=dsqrt(cappa2)                       !here the cappa is Debye–Hückel parameter 

!note: in mibpb_para.f90 func capppa() sets cappap=-cappa**2*eps1 as the negative squred modified DH 
eps=eps1/eps0
para1=332.0716d0
nfc=17
if (inl==1) then ! for nonlinear, setup matrix for Poisson eqn, save cappa into cappa_nl 
    cappa_nl=cappa**2
    cappa_nl=cappa_nl*eps1    !here cappa_nl is the squre of "modified" Debye–Hückel parameter
    cappa=0 !zero cappa as it is treated with nonlinearity 
endif

!################################################################################
! analytical surface
if (isf==0) then
	dx=dcel;	dy=dcel;	dz=dcel
	sphbd=4.d0;	sphbdyz=4.d0
! 	sphbd=4.6d0;	sphbdyz=4.6d0
	xleft=-sphbd;	xright=sphbd
	yleft=-sphbdyz;	yright=sphbdyz
	zleft=-sphbdyz;	zright=sphbdyz
	NX = INT( (XRIGHT - XLEFT)/DCEL) + 1
	NY = INT( (YRIGHT - YLEFT)/DCEL) + 1
	NZ = INT( (ZRIGHT - ZLEFT)/DCEL) + 1
	
    natm=1
    if (natm==1) then !sphere
        ALLOCATE(ATMPOS(3,natm), ATMRAD(natm))  
        rds=2.d0
        ATMRAD=rds
		atmpos(:,1)=(/0.d0,0.d0,0.d0/)

        write(*,'(a,3f10.4)') ' Atom positions:',atmpos
!         if (ipm==1) then
        	allocate(di_mom(3,natm),quad_mom(3,3,natm)) !yang natm==1
			di_mom(:,natm)=(/-0.32728d0, 0.00000d0, 0.10264d0/)	! real case (ln 2781)
			quad_mom(:,1,natm)=(/-0.33730d0, 0.00000d0, 0.41510d0/) !real case
			quad_mom(:,2,natm)=(/0.00000d0, 0.25185d0, 0.00000d0/)
			quad_mom(:,3,natm)=(/0.41510d0, 0.00000d0, 0.08545d0/)
	! 		di_mom(:,natm)=(/0.0d0,1.0d0,0.0d0/)
	! 		quad_mom(1,:,natm)=(/1.0d0, 0.0d0, 0.0d0/) !real case
	! 		quad_mom(2,:,natm)=(/0.0d0, -1.0d0, 0.0d0/)
	! 		quad_mom(3,:,natm)=(/0.0d0, 0.0d0, 0.0d0/)	
! 		endif	
    else 
        ALLOCATE(ATMPOS(3,natm), ATMRAD(natm))  
        rds=2.d0
        ATMRAD=rds
        atmpos(:,1)=(/xa,0.d0,0.d0/)
        atmpos(:,2)=(/xb,0.d0,0.d0/)
        write(*,'(a,6f10.4)') ' Atom positions:',atmpos
    endif 
endif


!read msms interface information by using dabao's code
if (isf==1)then
    call readin
!     if (natm==1) then 
! 		rds=atmrad(natm)
!     endif
    call readatom_gridgen
    call setirrpt2(imsmserr)
!     rslt=system('rm '//pathname(1:lenpath)//fname(1:lenfname)//'.xyzr')
!     rslt=system('rm '//pathname(1:lenpath)//fname(1:lenfname)//'.vert')
!     rslt=system('rm '//pathname(1:lenpath)//fname(1:lenfname)//'.face')
elseif (isf == 4) then
	call eses_readin
endif

! setup charge position (in topology.f90) 
call setChgPos
if (isf == 4) then
	allocate(irrpts(nx,ny,nz))
	irrpts=0;
else
	allocate(x(nx),y(ny),z(nz),io(nx,ny,nz),irrpts(nx,ny,nz))
	x=0.d0;		y=0.d0;		z=0.d0;		io=0;		irrpts=0;
endif
! allocation and initialization
allocate(bftc(nx*ny*nz),b(nx*ny*nz))
allocate(u(nx*ny*nz),uexa(nx*ny*nz),u0(nx*ny*nz),phi0(nx*ny*nz), phi0e_d(3,6,nx*ny*nz))
allocate(phi_s(nx*ny*nz),phi_v(nx*ny*nz),phi_rxn(nx*ny*nz))
bftc=0.d0;	phi0=0.d0;	b=0.d0;		phi0e_d=0.d0
u0=0.d0;	u=0.d0;		uexa=0.d0;	
!****************************************************
if (isf==0) then
	do i=1,nx
		x(i)=xleft+(i-1)*dx
	enddo
	do i=1,ny
		y(i)=yleft+(i-1)*dy
	enddo
	do i=1,nz
		z(i)=zleft+(i-1)*dz
	enddo
	xright=x(nx);	yright=y(ny);	zright=z(nz)
	call findio
	call findirrpts
	nirr=sum(irrpts)
	print *,'nirr= ',nirr
	call setindx	!set index for the array
endif

!****************************************************
allocate(ftpts(nirr,6,nfc),ftc(nirr,6))
allocate(iftpts(nirr,6,3*nfc),iftwrg(nirr,6))
allocate(uirr(nirr),ueirr(nirr),irrdrc(6,nirr),isrf(0:6,nirr))
ftc=0.d0;	ftpts=0.d0;		!cftc=0.d0;	
iftpts=0;	iftwrg=0;		isrf=0;			
uirr=0.d0;	ueirr=0.d0;		irrdrc=0;	

!####################################################
!Interface for getting information from daobao's code
if (isf==1 .or. isf==2) then
	io=-lptio;	
	irrpts=irrpt
	x=xc;		y=yc;			z=zc
	call readin_dist_angle
	!call output_coodinate
	deallocate(lptio,irrpt,xc,yc,zc)
elseif (isf == 4) then
	irrpts=irrpt
	call readin_dist_angle
	deallocate(irrpt)
endif 

!####################################################
write(*,'(a,2f10.4,a,2f10.4,a,2f10.4,a)') ' Domain: [',xleft,xright,'] [',yleft,yright,'] [',zleft,zright,']'
write(*,*) 'Dimentions:',nx,ny,nz
! allocate storage space for fictitious points used for surface point
call testsurf
allocate(srfpts(nsrfdrn,nfc),isrfpts(nsrfdrn,3*nfc),srfftc(nsrfdrn))
srfpts=0.d0; isrfpts=0; srfftc=0.d0

call irr_interface

!****************kirkwood
if (ibd==1) then
! Set Parameters for kirkwood case; increase nt if charges are closed to surface
nt=10
call lagendini
print *,'Finishing initialization of spherical harmonics ...'
endif
!*****************
if (ibd<2) then
    call tsoln
endif

print *,'set the RHS of the PBE ...'
call Setb;	bftc=b 

print *,'mibpb scheme for fictitious point ...'
call setirrpts

if (nirrall*2 .ne. sum(irrdrc)) then 
	call fillmisspts
endif
call testftpts

print *,'form matrix for mibpb ...'
nreg=nx*ny*nz-nirr
nmax=nirr*nfc*3+nreg*7+nfc*nsrfdrn*3
allocate(sa(nmax),ijka(nmax))
call formmx

do i=1,nmax
	if (ijka(i)==0) then 
		nmax=i
		goto 100
	endif
enddo
100 continue


allocate(sb(nmax),ijkb(nmax))
sb=sa(1:nmax); ijkb=ijka(1:nmax)
deallocate(sa,ijka)
allocate(sa(nmax),ijka(nmax))
sa=sb;	ijka=ijkb
if (inl==0) deallocate(sb,ijkb)

call cpu_time(cpu2)

! parameter for biconjugate method
thresh=1.d-20
itol=2
itmax=1000000
tol=1.d-10

! call test_intp(3.333d0,5.555d0,7.777d0,test_intp_grad0,test_intp_grad1,test_intp_grad2)
! print *, "test_intp_grad0 is ", test_intp_grad0
! print *, "test_intp_grad1 is ", test_intp_grad1(1), test_intp_grad1(2), test_intp_grad1(3)
! print *, "test_intp_grad2 is ", test_intp_grad2(1,1), test_intp_grad2(1,2), test_intp_grad2(1,3)
! print *, "test_intp_grad2 is ", test_intp_grad2(2,2), test_intp_grad2(2,3), test_intp_grad2(3,3)
! print *, "test_function value at 3.333d0,5.555d0,7.777d0 is ", test_function(3.333d0,5.555d0,7.777d0)

! print *, "test_function grad value at 3.333d0,5.555d0,7.777d0 is ", -sin(3.333d0)*sin(5.555d0)*exp(7.777d0), &
! cos(3.333d0)*cos(5.555d0)*exp(7.777d0), cos(3.333d0)*sin(5.555d0)*exp(7.777d0)
! print *, "test_function gradgrad value at 3.333d0,3.333d0,3.333d0 is ",-cos(3.333d0)*sin(5.555d0)*exp(7.777d0), &
! -sin(3.333d0)*cos(5.555d0)*exp(7.777d0), -sin(3.333d0)*sin(5.555d0)*exp(7.777d0)
! print *, "test_function gradgrad value at 3.333d0,3.333d0,3.333d0 is ",-cos(3.333d0)*sin(5.555d0)*exp(7.777d0), &
! cos(3.333d0)*cos(5.555d0)*exp(7.777d0),cos(3.333d0)*sin(5.555d0)*exp(7.777d0)

! print *, "test_function grad value at 3.333d0,5.555d0,7.777d0 is ", -3.333d0*(test_function(3.333d0,5.555d0,7.777d0)**3), &
! -5.555d0*(test_function(3.333d0,5.555d0,7.777d0)**3), -7.777d0*(test_function(3.333d0,5.555d0,7.777d0)**3)
! print *, "test_function gradgrad value at 3.333d0,3.333d0,3.333d0 is ",&
! ((3.333d0)**(2.d0))*3.d0*(test_function(3.333d0,5.555d0,7.777d0)**5)-(test_function(3.333d0,5.555d0,7.777d0)**3), &
! 3.333d0*5.555d0*3.d0*(test_function(3.333d0,5.555d0,7.777d0)**5), &
! 3.333d0*7.777d0*3.d0*(test_function(3.333d0,5.555d0,7.777d0)**5)
! print *, "test_function gradgrad value at 3.333d0,3.333d0,3.333d0 is ", & 
! ((5.555d0)**(2.d0))*3.d0*(test_function(3.333d0,5.555d0,7.777d0)**5)-(test_function(3.333d0,5.555d0,7.777d0)**3), &
! 5.555d0*7.777d0*3.d0*(test_function(3.333d0,5.555d0,7.777d0)**5), &
! ((7.777d0)**(2.d0))*3.d0*(test_function(3.333d0,5.555d0,7.777d0)**5)-(test_function(3.333d0,5.555d0,7.777d0)**3)


call linbcg(nx*ny*nz,bftc,U,itol,tol,itmax,iter,err)
! print *, "U is ", U

if (ipm==1) then
	print *, "solving point multipole PB equation with AMOEBA ..."
	call ugrad_ext	!yang
endif

if (inl==1) then
    print *,'solving nonlinear PB equation ...'
    u=0.d0 !use zero initial condition
    call nonlinearSolver
endif

if (icg==2) then
! 	print *, "ipm is ", ipm
	call solenergy(engy)
	do k=1,nz
		do j=1,ny
			do i=1,nx
				ijk=id3d(i,j,k)
				phi_v(ijk)=phi_star(x(i),y(j),z(k))
				if (io(i,j,k)<=0) then
	                phi_rxn(ijk)=u(ijk)   
				    u(ijk)=phi_rxn(ijk)+phi_v(ijk)
				endif
			enddo
		enddo
	enddo
! 	print *, "u "
	phi_s=u;	
endif
	
call cpu_time(cpu3)

if (ibd==0 .or. ibd==1) then
! The output of the error
	do i=1,nirr
		ii=irrxyz(1,i)
		jj=irrxyz(2,i)
		kk=irrxyz(3,i)
		ueirr(i)=Uexa((kk-1)*nx*ny+(jj-1)*nx+ii)
		uirr(i)=U((kk-1)*nx*ny+(jj-1)*nx+ii)
	enddo
	call error1(uirr,ueirr,e1,e2,e3,nirr,dcel)
	print *,'interface error= ', e1,e2,e3
	!print *,'max interface error position', irrxyz(:,maxloc(abs(ueirr-uirr)))
	call ERROR(U,Uexa,CL_inf,CL_2,nx*ny*nz,dcel)
	Print *,'Total Error = ', CL_inf,CL_2
	!call line2cube(ii,jj,kk,maxloc(abs(u-uexa)))
	print *,'max error position', ii,jj,kk
	print *,x(ii),y(jj),z(kk)
	print *,io(ii,jj,kk),irrpts(ii,jj,kk)
endif


! Output potential in pot.dx file in apbs format
ipt=1
if (ipt==1) then
    u=u*para1/0.5961634386 
    open(102,file="pot.dx")
    write(102,'(a)') "# Data from MIBPB"
    write(102,'(a)') '#'
    write(102,'(a)') '# POTENTIAL (kT/e)'
    write(102,'(a)') '#'
    write(102,'(a,3i4)') 'object 1 class gridpositions counts ', nx, ny, nz
    write(102,'(a,3ES14.6E2)') 'origin ', minval(x),minval(y),minval(z)
    write(102,'(a,3ES14.6E2)') 'delta ', dx, 0.0, 0.0
    write(102,'(a,3ES14.6E2)') 'delta ', 0.0, dy, 0.0
    write(102,'(a,3ES14.6E2)') 'delta ', 0.0, 0.0, dz
    write(102,'(a,3i4)') 'object 2 class gridconnections counts ', nx, ny, nz
    write(102,'(a,i6,a)') 'object 3 class array type double rank 0 items ',nx*ny*nz, ' data follows'
    ijk=0; u1=0.d0; u2=0.d0; u3=0.d0
    do i=1,nx
        do j=1,ny
            do k=1,nz
                ijk=ijk+1
                u1=u2; u2=u3; u3=u((k-1)*nx*ny+(j-1)*nx+i)

                if (mod(ijk,3)==0) then
                    write(102,'(3ES14.6E2)') real(u1),real(u2),real(u3)
                endif
            enddo
        enddo
    enddo
    if (mod(ijk,3) == 2) then
        write(102,*) real(u2), real(u3)
    elseif (mod(ijk,3) == 1) then
        write(102,*) real(u3)
    endif
    close(102)
endif



call cpu_time(cpu4)
Print *,'Time for generating fict. pts:=', real(cpu2-cpu1)
Print *,'Time for solving pbe via MIB :=', real(cpu3-cpu2)
Print *,'Total Cpu Time	              :=', real(cpu4-cpu0)

open(103,file="output.txt",status='unknown',action='write',form='formatted',position="append")
write(103,*) engy, cpu4-cpu0, iter
close(103)


end program mibpb_newReg 
!--------------------------------------------------------------------------------------

!!!! yang
! compute the gradient(grad) and gradgrad of exact solution u
subroutine ugrad_ext
use comdata
use molecule
use pbedata
implicit double precision(a-h,o-z)
real*8 phibar_intp_grad1(3),phibar_intp_grad2(3,3)
allocate(ugrad(3,natm),ugradgrad(3,3,natm)) !yang

! if (ibd==0) then
! find the index of the exact center 

! 	ii=int((nx+1)/2)
! 	ij=int((ny+1)/2)
! 	ik=int((nz+1)/2)
! 	ugradx=(u(id3d(ii+1,ij,ik))-u(id3d(ii-1,ij,ik)))/(2.d0*dx)
! 	ugrady=(u(id3d(ii,ij+1,ik))-u(id3d(ii,ij-1,ik)))/(2.d0*dy)
! 	ugradz=(u(id3d(ii,ij,ik+1))-u(id3d(ii,ij,ik-1)))/(2.d0*dz)
! 	ugrad(:,1)=(/ugradx,ugrady,ugradz/)

! 	ugradxx=(u(id3d(ii+1,ij,ik))-2.d0*u(id3d(ii,ij,ik))+u(id3d(ii-1,ij,ik)))/(dx**2)
! 	ugradxy=((u(id3d(ii+1,ij+1,ik))-u(id3d(ii+1,ij-1,ik)))-(u(id3d(ii-1,ij+1,ik))-u(id3d(ii-1,ij-1,ik))))/(4.d0*dx*dy)
! 	ugradxz=((u(id3d(ii+1,ij,ik+1))-u(id3d(ii+1,ij,ik-1)))-(u(id3d(ii-1,ij,ik+1))-u(id3d(ii-1,ij,ik-1))))/(4.d0*dx*dz)
! 	ugradyy=(u(id3d(ii,ij+1,ik))-2.d0*u(id3d(ii,ij,ik))+u(id3d(ii,ij-1,ik)))/(dy**2)
! 	ugradyz=((u(id3d(ii,ij+1,ik+1))-u(id3d(ii,ij+1,ik-1)))-(u(id3d(ii,ij-1,ik+1))-u(id3d(ii,ij-1,ik-1))))/(4.d0*dy*dz)
! 	ugradzz=(u(id3d(ii,ij,ik+1))-2.d0*u(id3d(ii,ij,ik))+u(id3d(ii,ij,ik-1)))/(dz**2)

! 	! see bottcher and schnieders
! 	ugradxx=ugradxx/3.d0
! 	ugradxy=ugradxy/3.d0
! 	ugradxz=ugradxz/3.d0
! 	ugradyy=ugradyy/3.d0
! 	ugradyz=ugradyz/3.d0
! 	ugradzz=ugradzz/3.d0

! 	ugradgrad(:,1,1)=(/ugradxx,ugradxy,ugradxz/)
! 	ugradgrad(:,2,1)=(/ugradxy,ugradyy,ugradyz/)
! 	ugradgrad(:,3,1)=(/ugradxz,ugradyz,ugradzz/)

! 	print *, "ipm at ugrad is ", ipm
	do i=1,natm
		call phibar_intp_gradgrad(chrpos(1,i),chrpos(2,i),chrpos(3,i),phibar_intp_grad1,phibar_intp_grad2)
		ugradx = phibar_intp_grad1(1)
		ugrady = phibar_intp_grad1(2)
		ugradz = phibar_intp_grad1(3)

		ugrad(:,natm)=(/ugradx,ugrady,ugradz/)

		ugradxx=phibar_intp_grad2(1,1)
		ugradxy=phibar_intp_grad2(1,2)
		ugradxz=phibar_intp_grad2(1,3)
		ugradyy=phibar_intp_grad2(2,2)
		ugradyz=phibar_intp_grad2(2,3)
		ugradzz=phibar_intp_grad2(3,3)

		ugradgrad(:,1,natm)=(/ugradxx,ugradxy,ugradxz/)
		ugradgrad(:,2,natm)=(/ugradxy,ugradyy,ugradyz/)
		ugradgrad(:,3,natm)=(/ugradxz,ugradyz,ugradzz/)

	enddo
	print *, "ugrad(:,1) is ", ugrad(:,1)

end subroutine ugrad_ext


!-----------------------------------------

subroutine phibar_intp_gradgrad(xx,yy,zz,phibar_intp_grad1,phibar_intp_grad2)
use comdata
use molecule
use pbedata
implicit double precision(a-h,o-z)
real*8 xwt(-1:1,0:2),ywt(-1:1,0:2),zwt(-1:1,0:2),aa(3),wt(27),value(27,0:2),w,ax(nx),by(ny),cz(nz),bb(3,0:2),cc(3)
integer ixyz(1)
real*8 phibar_intp_grad1(3), phibar_intp_grad2(3,3),wtx(27),wty(27),wtz(27),wtxx(27),wtxy(27),wtxz(27)
real*8 wtyy(27),wtyz(27),wtzz(27)
real*8 gradx, grady, gradz, gradgradxx, gradgradxy, gradgradxz, gradgradyy, gradgradyz, gradgradzz
gradx=0.d0
grady=0.d0
gradz=0.d0
gradgradxx=0.d0
gradgradxy=0.d0
gradgradxz=0.d0
gradgradyy=0.d0
gradgradyz=0.d0
gradgradzz=0.d0

ixyz=minloc(abs(x-xx))
ix=ixyz(1)
ixyz=minloc(abs(y-yy))
iy=ixyz(1)
ixyz=minloc(abs(z-zz))
iz=ixyz(1)

aa=x(ix-1:ix+1)
call weights(xx,aa,2,2,2,xwt)
aa=y(iy-1:iy+1)
call weights(yy,aa,2,2,2,ywt)
aa=z(iz-1:iz+1)
call weights(zz,aa,2,2,2,zwt)

ijk=0
do k=-1,1
	do j=-1,1
		do i=-1,1
			ijk=ijk+1
			wt(ijk)=xwt(i,0)*ywt(j,0)*zwt(k,0)

			wtx(ijk)=xwt(i,1)*ywt(j,0)*zwt(k,0)
			wty(ijk)=xwt(i,0)*ywt(j,1)*zwt(k,0)
			wtz(ijk)=xwt(i,0)*ywt(j,0)*zwt(k,1)

			wtxx(ijk)=xwt(i,2)*ywt(j,0)*zwt(k,0)
			wtxy(ijk)=xwt(i,1)*ywt(j,1)*zwt(k,0)
			wtxz(ijk)=xwt(i,1)*ywt(j,0)*zwt(k,1)

			wtyy(ijk)=xwt(i,0)*ywt(j,2)*zwt(k,0)
			wtyz(ijk)=xwt(i,0)*ywt(j,1)*zwt(k,1)
			wtzz(ijk)=xwt(i,0)*ywt(j,0)*zwt(k,2)

			if (io(ix+i,iy+j,iz+k)>0) then
! 				print *, "special cases in phibar_intp_gradgrad"
				if (ix+i>=4) then
					if (io(ix+i-1,iy+j,iz+k)<=0 .and. io(ix+i-2,iy+j,iz+k)<=0 .and. io(ix+i-3,iy+j,iz+k)<=0) then
						aa=x(ix+i-1:ix+i-3:-1)
						call weights(x(ix),aa,2,2,2,bb)
						Do ii=1,3
							cc(ii)=phi0(id3d(ix+i-ii,iy+j,iz+k))+u(id3d(ix+i-ii,iy+j,iz+k))
						enddo
						value(ijk,0)=dot_product(bb(:,0),cc)
						goto 100
					endif
				elseif(ix+i<=nx-3) then
					if(io(ix+i+1,iy+j,iz+k)<=0 .and. io(ix+i+2,iy+j,iz+k)<=0 .and. io(ix+i+3,iy+j,iz+k)<=0) then
						aa=x(ix+i+1:ix+i+3)
						call weights(x(ix),aa,2,2,2,bb)
						Do ii=1,3
							cc(ii)=phi0(id3d(ix+i+ii,iy+j,iz+k))+u(id3d(ix+i+ii,iy+j,iz+k))
						enddo
						value(ijk,0)=dot_product(bb(:,0),cc)
						goto 100
					endif
				elseif(iy+j>=4) then
					if(io(ix+i,iy+j-3,iz+k)<=0 .and. io(ix+i,iy+j-2,iz+k)<=0 .and. io(ix+i,iy+j-1,iz+k)<=0) then
						aa=y(iy+j-1:iy+j-3:-1)
						call weights(y(iy),aa,2,2,2,bb)
						Do ii=1,3
							cc(ii)=phi0(id3d(ix+i,iy+j-ii,iz+k))+u(id3d(ix+i,iy+j-ii,iz+k))
						enddo
						value(ijk,0)=dot_product(bb(:,0),cc)
						goto 100
					endif
				elseif(iy+j<=ny-3) then
					if(io(ix+i,iy+j+1,iz+k)<=0 .and. io(ix+i,iy+j+2,iz+k)<=0 .and. io(ix+i,iy+j+3,iz+k)<=0) then
						aa=y(iy+j+1:iy+j+3)
						call weights(y(iy),aa,2,2,2,bb)
						Do ii=1,3
							cc(ii)=phi0(id3d(ix+i,iy+j+ii,iz+k))+u(id3d(ix+i,iy+j+ii,iz+k))
						enddo
						value(ijk,0)=dot_product(bb(:,0),cc)
						goto 100
					endif
				elseif(iz+k>=4) then
					if(io(ix+i,iy+j,iz+k-3)<=0 .and. io(ix+i,iy+j,iz+k-2)<=0 .and. io(ix+i,iy+j,iz+k-1)<=0) then
						aa=z(iz+k-1:iz+k-3:-1)
						call weights(z(iz),aa,2,2,2,bb)
						Do ii=1,3
							cc(ii)=phi0(id3d(ix+i,iy+j,iz+k-ii))+u(id3d(ix+i,iy+j,iz+k-ii))
						enddo
						value(ijk,0)=dot_product(bb(:,0),cc)
						goto 100
					endif
				elseif(iz+k<=nz-3) then
					if(io(ix+i,iy+j,iz+k+1)<=0 .and. io(ix+i,iy+j,iz+k+2)<=0 .and. io(ix+i,iy+j,iz+k+3)<=0) then
						aa=z(iz+k+1:iz+k+3)
						call weights(z(iz),aa,2,2,2,bb)
						Do ii=1,3
							cc(ii)=phi0(id3d(ix+i,iy+j,iz+k+ii))+u(id3d(ix+i,iy+j,iz+k+ii))
						enddo
						value(ijk,0)=dot_product(bb(:,0),cc)
! 						value(ijk,1)=dot_product(bb(:,1),cc)
! 						value(ijk,2)=dot_product(bb(:,2),cc)
						goto 100
					endif
				else
					print *,i,j,k,'failed to interpolate'
				endif
			else
				value(ijk,0)=phi0(id3d(ix+i,iy+j,iz+k))+u(id3d(ix+i,iy+j,iz+k))
			endif
100 continue
		enddo
	enddo
enddo
! phibar_intp=dot_product(wt,value)

gradx=dot_product(wtx,value(:,0))
grady=dot_product(wty,value(:,0))
gradz=dot_product(wtz,value(:,0))

phibar_intp_grad1 = (/gradx,grady,gradz/)

gradgradxx=dot_product(wtxx,value(:,0))
gradgradxy=dot_product(wtxy,value(:,0))
gradgradxz=dot_product(wtxz,value(:,0))
gradgradyy=dot_product(wtyy,value(:,0))
gradgradyz=dot_product(wtyz,value(:,0))
gradgradzz=dot_product(wtzz,value(:,0))
! phibar_intp_grad2 = (/gradgradxx,gradgradxy,gradgradxz,gradgradyy,gradgradyz,gradgradzz/)

phibar_intp_grad2(:,1)=(/gradgradxx,gradgradxy,gradgradxz/)
phibar_intp_grad2(:,2)=(/gradgradxy,gradgradyy,gradgradyz/)
phibar_intp_grad2(:,3)=(/gradgradxz,gradgradyz,gradgradzz/)

End subroutine phibar_intp_gradgrad


!--------------------------------------------------------------------------------------
subroutine output_coodinate
use comdata
use molecule
use pbedata
implicit double precision(a-h,o-z)
OPEN(11,FILE="cood.dat")  
do i=1,nx
	do j=1,ny
		do k=1,nz
			write(11,'(3i4,2i2,i8,3f10.4)') i,j,k,io(i,j,k),irrpts(i,j,k),indirr(i,j,k),x(i),y(j),z(k)
		enddo
	enddo
enddo
Close(11)
OPEN(11,FILE="srfc.dat")  
do irr=1,nirr
	write(11,'(i8,3i4)') irr, irrxyz(:,irr)
enddo
Close(11)

OPEN(11,FILE="dist_angle.dat")  
do irr=1,nirr
	do ik=1,6
		write(11,'(i8,i2,5f12.8)') irr,ik,dist_angle(1:5,ik,irr)
	enddo
enddo
Close(11)


End
!---------------------------------------------------------------------------------------
subroutine setindx
use comdata
use pbedata
use molecule
implicit double precision(a-h,o-z)
allocate(irrxyz(3,nirr),indirr(nx,ny,nz))
irrxyz=0;	indirr=0;	icount=0
do i=1,nx
	do j=1,ny
		do k=1,nz
			if (irrpts(i,j,k)==1) then
				icount=icount+1
				indirr(i,j,k)=icount
				irrxyz(:,icount)=(/i,j,k/)				
			endif
		enddo
	enddo
enddo

end

!---------------------------------------------------------------------------------------
subroutine findio
use comdata
use pbedata
use molecule
implicit double precision(a-h, o-z)
do i=1,nx
	do j=1,ny
		do k=1,nz
			if (abs(varphi(x(i),y(j),z(k)))<=1.d-10) then
				io(i,j,k)=0
			elseif (varphi(x(i),y(j),z(k))+1.d-10<=0.d0) then
				io(i,j,k)=-1
			else
				io(i,j,k)=1
			endif
		enddo
	enddo
enddo

end subroutine findio
!-----------------------------------------------------------------------------------------
subroutine findirrpts
use comdata
use pbedata
use molecule
implicit double precision(a-h, o-z)
m=1
isum=0
irrpts=0
do i=1+m,nx-m
	do j=1+m,ny-m
		do k=1+m,nz-m
			if (io(i,j,k)==0 .or. (min(minval(io(i-m:i+m,j,k)),minval(io(i,j-m:j+m,k)), minval(io(i,j,k-m:k+m))) &
			*max(maxval(io(i-m:i+m,j,k)),maxval(io(i,j-m:j+m,k)),maxval(io(i,j,k-m:k+m)))<0)) then 
				irrpts(i,j,k)=1
				isum=isum+1
			else
				irrpts(i,j,k)=0
			endif
		enddo	
	enddo
enddo
end subroutine findirrpts

!-----------------------------------------------------------------
subroutine testsurf
use comdata
use pbedata
use molecule
implicit double precision(a-h,o-z)

do irr=1,nirr
	i=irrxyz(1,irr);	j=irrxyz(2,irr);	k=irrxyz(3,irr)
	if (io(i,j,k)==0 .or. irrpts(i,j,k)==-1) then
		isum=io(i-1,j,k)+io(i+1,j,k)+io(i,j-1,k)+io(i,j+1,k)+io(i,j,k-1)+io(i,j,k+1)
		if (isum==2) then 
			isrf(0,irr)=1;		nsrfdrn=nsrfdrn+2
		elseif (isum==4) then
			isrf(0,irr)=1;		nsrfdrn=nsrfdrn+1
		elseif (isum==-2) then
			isrf(0,irr)=-1;		nsrfdrn=nsrfdrn+2
		elseif (isum==-4) then
			isrf(0,irr)=-1;		nsrfdrn=nsrfdrn+1
		else
			isrf(0,irr)=-1;		nsrfdrn=nsrfdrn+3
		endif
	endif
enddo

End

!---------------------------------------------------------------------------------------
function xv(ix)
use comdata
implicit double precision(a-h,o-z)
xv=xleft+(ix-1)*(xright-xleft)/dx
end

function yv(iy)
use comdata
implicit double precision(a-h,o-z)
yv=yleft+(iy-1)*(yright-yleft)/dy
end

function zv(iz)
use comdata
implicit double precision(a-h,o-z)
zv=zleft+(iz-1)*(zright-zleft)/dz
end

!-----------------------------------------------------------------------------------------
SUBROUTINE ERROR(F,F_EXA, CL1,CL3,n,DELTA)
IMPLICIT DOUBLE PRECISION(A-H,O-Z)
DIMENSION F(1:n),F_EXA(1:n)
CL1=0.0d0
CL2=0.0d0
DO I=1,n
    CL=dABS(F(I)-F_EXA(I))
 	CL2=CL2+(CL)**2
    IF(CL1<CL)CL1=CL
END DO
CL3=(CL2/n)**0.5
RETURN
END
!-----------------------------------------------------------------------------
SUBROUTINE ERROR1(F,F_EXA, e1,e2,e3,n,DELTA)
IMPLICIT DOUBLE PRECISION(A-H,O-Z)
DIMENSION F(1:n),F_EXA(1:n)
e1=0.0d0
e2=0.0d0
e3=0.0d0
DO I=1,n
!	if (abs(f_exa(i)) < 1.d-15) print *,'relative error failed',f_exa(i)
	err0=ABS(F(I)-F_EXA(I))
	err1=ABS((F(I)-F_EXA(I))/f_exa(i))
	e3=e3+err1
    IF (e1<err0) e1=err0
	IF (e2<err1) e2=err1
END DO
e2=e2*100.d0
e3=e3/n*100.d0
RETURN

END

SUBROUTINE ERROR2(F,FEXA,e1,e2,N)
IMPLICIT DOUBLE PRECISION(A-H,O-Z)
DIMENSION F(N),FEXA(N)
e1=0.0
e2=0.0
DO I=1,N
    err1=ABS(F(I)-FEXA(I))
	err2=abs((f(i)-fexa(i))/fexa(i))
    if(e1<err1) e1=err1
    IF(e2<err2) e2=err2
END DO
RETURN
END


!------------------------------------------------------------------------------
Subroutine Tsoln
Use comdata
use pbedata
Implicit Double Precision(A-H, O-Z)
common /pi/ pi
real*8 pxyz(3)
Uexa=0.d0

do i=1,nx
    do j=1,ny
        do k=1,nz
            if (irrpts(i,j,k)==1) then !surface error only, comment if 3-d
                pxyz=(/x(i),y(j),z(k)/)
		        if (io(i,j,k) <= 0) then
			        if (ibd==1) then
				        call kirk_potential(io(i,j,k),pxyz,ptl,pmt,psv)
				        ! uexa(id3d(i,j,k))=pmt+psv !phi~
				        uexa(id3d(i,j,k))=ptl
			        else
				        uexa((k-1)*nx*ny+(j-1)*nx+i)=Un(x(i),y(j),z(k))
! 				        print *, "un at tsoln is "!, uexa((k-1)*nx*ny+(j-1)*nx+i)
! 				        print *,i,j,k,uexa((k-1)*nx*ny+(j-1)*nx+i)
			        endif
		        else
			        if (ibd==1) then
				        call kirk_potential(io(i,j,k),pxyz,ptl,pmt,psv)
				        uexa(id3d(i,j,k))=ptl
			        else
				        uexa((k-1)*nx*ny+(j-1)*nx+i)=Up(x(i),y(j),z(k))
! 				        print *, "up at tsoln is "!, uexa((k-1)*nx*ny+(j-1)*nx+i)
			        endif
		        endif
            endif
	    enddo
    enddo
enddo

End
!------------------------------------------------------------------------------
Subroutine Setb
Use comdata
Use pbedata
use molecule
Implicit Double Precision(A-H, O-Z)
common /pi/ pi
b=0.d0
do i=1,nx
	do j=1,ny
		do k=1,nz 
			if ((i==1) .or. (i==nx) .or. (j==1) .or. (j==ny) .or. (k==1) .or. (k==nz))then
				if (ibd==0) then
					b((k-1)*nx*ny+(j-1)*nx+i)=up(x(i),y(j),z(k))
				elseif(ibd==1) then
					b((k-1)*nx*ny+(j-1)*nx+i)=upi(i,j,k)
				else
					b((k-1)*nx*ny+(j-1)*nx+i)=bdcond(x(i),y(j),z(k))
				endif
			else
				if (ibd==0) then
					if (io(i,j,k)==-1 .or. io(i,j,k)==0) then
						b((k-1)*nx*ny+(j-1)*nx+i)=b((k-1)*nx*ny+(j-1)*nx+i)+fn(x(i),y(j),z(k))
					else
						b((k-1)*nx*ny+(j-1)*nx+i)=b((k-1)*nx*ny+(j-1)*nx+i)+fp(x(i),y(j),z(k))
					endif
				endif
			endif
		enddo
	enddo
enddo



End
!----------------------------------------------------------------------------------------
subroutine setirrpts
use comdata
use pbedata
use molecule
implicit double precision(a-h, o-z)
integer idx(3)
! trig(4) cos(theta), sin(theta), cos(phi), sin(phi)
xacc=1.d-10
intppts=0
sfpts=0.d0
nirrall=0
err_flux=0.d0


do irr=1,nirr
	i=irrxyz(1,irr); j=irrxyz(2,irr); k=irrxyz(3,irr)
	!if (io(i,j,k)==0 .or. irrpts(i,j,k)==-1) then		! For grids on interface 
	if (io(i,j,k)==0 ) then	
		call grdptwts(i,j,k)  
	elseif (irrpts(i,j,k)==1 .and. io(i,j,k)==-1) then	! Irregular points
		
		do ik=1,6
			call idx6(i,j,k,ik,idx)					! idx stores coodinates of the fictitious points outside
			if (io(idx(1),idx(2),idx(3))==1) then	! fictitious point is needed in this direction
			
				call findfctin(i,j,k,ik,0)
			endif
		enddo
	endif

enddo

!print *,'err_flux= ', err_flux
end

!------------------------------------------------------------------------------
subroutine findfctin(i,j,k,ik,iforce)
use comdata
use pbedata
use molecule
implicit double precision(a-h,o-z)
real*8	P(3,3), wts(2,16), phi, rjp(4), coiuas(2,2,2,3), flux(1:3)
real*8	pcts(4),h1(3),h2(3),p0p(3),p0n(3),p1p(3),p1n(3),s1(3),t1(3),xpts(2,2,3), trig(4)
integer iwts(16,3), ix1(3), ix2(3), intppts(2,6,3), intppts_temp(6,3), idx(3),idy(3),idz(3)
integer ixtp(2),ick(6),i0(3),iuas(2,2,3,3),idsn(3),idnsn(3,2)

	irr=indirr(i,j,k)
	irr0=irr
	xacc=1.d-10
	iwts=0			
	call idx6(i,j,k,ik,idx)				! idx stores coodinates of the fictitious points outside
	call idx6(i,j,k,7-ik,idy)			! idy stores the coodinates of the opposite of idx 

	select case (ik)
	case (1)
		pcts=z(k+1:k-2:-1)
		iwts(1:4,3)=(/k+1,k,k-1,k-2/)
	case (2)	
		pcts=y(j+1:j-2:-1)
		iwts(1:4,2)=(/j+1,j,j-1,j-2/)
	case (3)	
		pcts=x(i+1:i-2:-1)
		iwts(1:4,1)=(/i+1,i,i-1,i-2/)
	case (4)
		pcts=x(i-1:i+2:1)
		iwts(1:4,1)=(/i-1,i,i+1,i+2/)
	case (5)
		pcts=y(j-1:j+2:1)
		iwts(1:4,2)=(/j-1,j,j+1,j+2/)
	case (6)
		pcts=z(k-1:k+2:1)
		iwts(1:4,3)=(/k-1,k,k+1,k+2/)
	end select

		irrdrc(ik,indirr(i,j,k))=1
		irrdrc(7-ik,indirr(idx(1),idx(2),idx(3)))=1
		
		! Can't get fictitious pts in this direction		
		if (io(idy(1),idy(2),idy(3))==1 .and. iforce==0) then	
			iftwrg(indirr(i,j,k),7-ik)=1
			iftwrg(indirr(idx(1),idx(2),idx(3)),ik)=1
			iftwrg(indirr(i,j,k),ik)=1
			iftwrg(indirr(idy(1),idy(2),idy(3)),7-ik)=1
			return
		endif

		call idx6(idx(1),idx(2),idx(3),ik,idz)
		
		if (io(idz(1),idz(2),idz(3))==-1 .and. iforce==0) then	
			iftwrg(indirr(i,j,k),7-ik)=1
			iftwrg(indirr(idx(1),idx(2),idx(3)),ik)=1
			return
		endif
		ix1=(/min(i,idx(1)),min(j,idx(2)),min(k,idx(3))/)
		ix2=(/max(i,idx(1)),max(j,idx(2)),max(k,idx(3))/)

		if (ix1(1)==ix2(1) .and. ix1(2)==ix2(2)) then				! change in z direction
			xs=x(i);	ys=y(j);	zs=rrintf(ik,irr);	
			idsn=(/0,0,1/);	idnsn(:,1)=(/1,0,0/);	idnsn(:,2)=(/0,1,0/)
		elseif (ix1(1)==ix2(1) .and. ix1(3)==ix2(3)) then			! change in y direction
			xs=x(i);	zs=z(k);	ys=rrintf(ik,irr)
			idsn=(/0,1,0/);	idnsn(:,1)=(/0,0,1/);	idnsn(:,2)=(/1,0,0/)
		else														! change in x direction
			ys=y(j);	zs=z(k);	xs=rrintf(ik,irr)
			idsn=(/1,0,0/);	idnsn(:,1)=(/0,1,0/);	idnsn(:,2)=(/0,0,1/)	
		endif
				
		iside=0
		ijust=0
		call intpraxpts(i,j,k,xs,ys,zs,ix1,ix2,intppts,iside,ijust)				!???????????

		!if (ijust .ne. 1) print *,i,j,k,ijust,' 6 aux points wrong'
		call test_intppts(i,j,k,intppts,iwrong)

		!***************************************************************************
		if (iforce==4) then			! use another way to search auxilary points
			if (iwrong==1) then
				ia=1;ib=1
			elseif (iwrong==2) then
				ia=2;ib=2
			elseif(iwrong==3) then
				ia=1;ib=2
			endif

			do i4=ia,ib
				isn=idx(1)-i;	jsn=idx(2)-j;	ksn=idx(3)-k
				call searchauxi4(i4,i,j,k,isn,jsn,ksn,idsn,idnsn,iuas,coiuas,-1,ijust)
				
				if (ijust .ne. 1) call searchauxi4(i4,i,j,k,isn,jsn,ksn,idsn,idnsn,iuas,coiuas,1,ijust)
				do jj=1,2
					do kk=1,3
						if (idsn(2)==1) then
							intppts(3-i4,(jj-1)*3+kk,:)=iuas(i4,jj,kk,:)
						else
							intppts(i4,(jj-1)*3+kk,:)=iuas(i4,jj,kk,:)	
						endif
					enddo
				enddo
			enddo
			call test_intppts(i,j,k,intppts,iwrong)
			if (iwrong==1) print *,i,j,k,'reach the boundary'
		endif
		!**********************************************************************
		if (ijust==0 .and. iforce==1) then
			print *,i,j,k,ik, ' double troubles, need patch 3'
			iforce=2
			return
		endif

		if (ijust==0 .and. iforce==0) then
			iftwrg(indirr(i,j,k),7-ik)=1
			iftwrg(indirr(idx(1),idx(2),idx(3)),ik)=1
			print *,i,j,k,ik,'can not find all auxilary points'
			return
		endif
		
		if (iforce==0) nirrall=nirrall+1							! nirrall stores number of smooth irr pts

		do ii=1,2 !Inside or outside of Points extroplatable, ii is the two planes
			ixtp(ii)=io(intppts(ii,1,1),intppts(ii,1,2),intppts(ii,1,3))
		enddo

		if (ix1(2)==ix2(2) .and. ix1(3)==ix2(3)) then				! change in x direction
			if (isf==1 .or. isf==2 .or. isf==4) then 
				if (mcx(irr)==-1) irr0=indirr(i-1,j,k)
				trig=clocal(:,1,irr0)
			endif
			do ii=1,2
				xpts(ii,1,1:3)=x(intppts(ii,1:3,1))					! xpts?????????????
				xpts(ii,2,1:3)=x(intppts(ii,4:6,1))
			enddo
			
			s1=(/ys,y(intppts(1,1,2)),y(intppts(1,4,2))/)
			t1=(/zs,z(intppts(2,1,3)),z(intppts(2,4,3))/)
					
			call findfdwts(xs,ys,zs,pcts,xpts,p0p,p0n,p1p,p1n,s1,t1)	
			iwts(1:4,2)=j;	iwts(1:4,3)=k
			idn=1				
		elseif (ix1(1)==ix2(1) .and. ix1(3)==ix2(3)) then			! change in y direction
			if (isf==1 .or. isf==2 .or. isf==4) then
				if (mcy(irr)==-1) irr0=indirr(i,j-1,k) 
				trig=clocal(:,2,irr0)
			endif
			do ii=1,2	
				xpts(ii,1,1:3)=y(intppts(ii,1:3,2))
				xpts(ii,2,1:3)=y(intppts(ii,4:6,2))
			enddo
			s1=(/xs,x(intppts(1,1,1)),x(intppts(1,4,1))/)
			t1=(/zs,z(intppts(2,1,3)),z(intppts(2,4,3))/)
			
			call findfdwts(ys,xs,zs,pcts,xpts,p0p,p0n,p1p,p1n,s1,t1)	
			iwts(1:4,1)=i;	iwts(1:4,3)=k
			idn=3				
		else														! change in z direction
			if (isf==1 .or. isf==2 .or. isf==4) then
				if (mcz(irr)==-1) irr0=indirr(i,j,k-1)
				trig=clocal(:,3,irr0)
			endif
			do ii=1,2
				xpts(ii,1,1:3)=z(intppts(ii,1:3,3))
				xpts(ii,2,1:3)=z(intppts(ii,4:6,3))
			enddo
			s1=(/xs,x(intppts(1,1,1)),x(intppts(1,4,1))/)
			t1=(/ys,y(intppts(2,1,2)),y(intppts(2,4,2))/)
			call findfdwts(zs,xs,ys,pcts,xpts,p0p,p0n,p1p,p1n,s1,t1)
			iwts(1:4,1)=i;	iwts(1:4,2)=j
			idn=5				
		endif

		if (isf==0) then
			call Transfmx(xs,ys,zs, P)
			!call Transfmx_Sining(xs,ys,zs,P)
			!ctht=costheta(xs,ys,zs) 
			!stht=sintheta(xs,ys,zs)
			!cphi=cosphi(xs,ys,zs)
			!sphi=sinphi(xs,ys,zs)
		elseif (isf==1 .or. isf==2 .or. isf==4) then
			call Transfmx_Dabao(trig,P)
		elseif(isf==3) then
			sphi=dist_angle(1,ik,indirr(i,j,k))
			cphi=dist_angle(2,ik,indirr(i,j,k))
			ctht=dist_angle(3,ik,indirr(i,j,k))
			stht=dist_angle(4,ik,indirr(i,j,k))
			call coortrans(p,sphi,cphi,stht,ctht)
		endif
		
		if (icg==0) then				! Directly set the jump condition
		    call findjumps(xs,ys,zs,P,rjp)		! Use this if the true solution known
		elseif (icg==2) then
                    call findRegJumps(xs,ys,zs,P,rjp,eps0,ipm)
		endif
		iwts(5:10,1:3)=intppts(1,1:6,1:3);	iwts(11:16,1:3)=intppts(2,1:6,1:3)			
					
		do l=1,16
			iftpts(indirr(i,j,k),7-ik,3*l-2:3*l)=iwts(l,1:3)
			iftpts(indirr(idx(1),idx(2),idx(3)),ik,3*l-2:3*l)=iwts(l,1:3)
		enddo
		

		itrigtest=0
		
		if	((abs(p(1,1))<1.d-10 .and. abs(p(1,2))<1.d-10) .or.	&	!phi=0, z-axis 
			 (abs(p(1,1))<1.d-10 .and. abs(p(1,3))<1.d-10) .or.	&	!cphi=0, cth=0, y-axis 
			 (abs(p(1,2))<1.d-10 .and. abs(p(1,3))<1.d-10)) then	!cphi=0, sth=0, x-axis
			itrigtest=1
		endif
		
		if (itrigtest == 1) then	
			if (ik==1 .or. ik==2 .or. ik==3) then					!Ux <> Un when 1,2,3
				p1p=-p1p;	p1n=-p1n
			endif
			call findtgwts(idx,ik,i,j,k,xs,ys,zs,p0p,p0n,p1p,p1n,rjp(1),rjp(2))
		else
			call findprtcts(i,j,k,p,betan(xs,ys,zs),betap(xs,ys,zs),cpp,cpn,cs,ct,const,ixtp,idn,rjp) !????
			call findwts(idx,ik,i,j,k,wts,cpp,cpn,cs,ct,const,xpts,p0p,p0n,p1p,p1n,s1,t1,rjp(1),ixtp) !???
		endif

End
!------------------------------------------------------------------------------
subroutine test_intppts(i,j,k,intppts,iwrong)
use pbedata
implicit double precision(a-h,o-z)
integer intppts(2,6,3),isum(2)

iwrong=0
do ii=1,2
	itemp=io(intppts(ii,1,1),intppts(ii,1,2),intppts(ii,1,3))
	do jj=1,6
		if (io(intppts(ii,jj,1),intppts(ii,jj,2),intppts(ii,jj,3)) .ne. itemp) then 
			iwrong=iwrong+ii
			!print *,i,j,k,ii,jj, 'is wrong'
			goto 103
		endif
	enddo
103 continue
enddo


End

!------------------------------------------------------------------------------
subroutine findfdwts(ps,ss,ts,pcts,xpts,p0p,p0n,p1p,p1n,s1,t1)				
implicit double precision(a-h,o-z)				
real*8 ps,ss,pcts(4),xpts(2,2,3),p0p(3),p0n(3),p1p(3),p1n(3),s1(3), t1(3),tp(0:2),xx(0:2),c(0:2,0:1)
!label  FORMAT('First Descriptor')
										! ps, ss, ts are the coordinates of the interface point

z=ps;	xx(0:2)=pcts(1:3);	m=1			! Negative function and derivative in primary direction
call weights (z,xx,2,2,m,c)
p0n=c(:,0);	p1n=c(:,1)

z=ps;	xx(0:2)=pcts(4:2:-1);	m=1		! Positive function and derivative in Primary direction
call weights (z,xx,2,2,m,c)
p0p=c(:,0);	p1p=c(:,1)

z=ss;	xx(0:2)=s1(1:3);	m=1			! Derivative in the first secondary direction
call weights (z,xx,2,2,m,c)
s1=c(:,1)

z=ts;	xx(0:2)=t1(1:3);	m=1			! Derivative in the second secondary direction
call weights (z,xx,2,2,m,c)
t1=c(:,1)

z=ps;	xx(0:2)=xpts(1,1,:);	m=0		! Interpolation of the 1st auxilary point in the 1st secondary direction. 
call weights (z,xx,2,2,m,tp)
xpts(1,1,:)=tp

z=ps;	xx(0:2)=xpts(1,2,:);	m=0		! Interpolation of the 2nd auxilary point in the 1st secondary direction. 
call weights (z,xx,2,2,m,tp)
xpts(1,2,:)=tp

z=ps;	xx(0:2)=xpts(2,1,:);	m=0		! Interpolation of the 1st auxilary point in the 2nd secondary direction. 
call weights (z,xx,2,2,m,tp)
xpts(2,1,:)=tp

z=ps;	xx(0:2)=xpts(2,2,:);	m=0		! Interpolation of the 2nd auxilary point in the 2nd secondary direction. 
call weights (z,xx,2,2,m,tp)
xpts(2,2,:)=tp
end
!-----------------------------------------------------------------------------
Subroutine Transfmx(x,y,z,P) 
Implicit Double Precision(A-H, O-Z)
real*8 x, y, z, P(3,3), xi(3), eta(3), tao(3), cx(3), cy(3), cz(3)

cx=(/1.d0,0.d0,0.d0/)
cy=(/0.d0,1.d0,0.d0/)
cz=(/0.d0,0.d0,1.d0/)

px=varphix(x,y,z)
py=varphiy(x,y,z)
pz=varphiz(x,y,z)

xi=(/px, py, pz/)/sqrt(px**2+py**2+pz**2)

if (px**2+py**2 <= 1.d-10) then 
	eta=(/pz,0.d0,-px/)/sqrt(px**2+pz**2)
	tao=(/-px*py, px**2+pz**2, -py*pz/)/sqrt((-px*py)**2 + (px**2+pz**2)**2 + (-py*pz)**2)
else
	eta=(/py,-px,0.d0/)/sqrt(px**2+py**2)
	tao=(/px*pz, py*pz, -px**2-py**2/)/sqrt((px*pz)**2 + (py*pz)**2 + (-px**2-py**2)**2)
endif

!Construct the transform matrix such that (xi, eta, tao)'=P(x,y,z)'
p(1,1)=DrcCos(cx,xi,3)
p(1,2)=DrcCos(cy,xi,3)
p(1,3)=DrcCos(cz,xi,3)
p(2,1)=DrcCos(cx,eta,3)
p(2,2)=DrcCos(cy,eta,3)
p(2,3)=DrcCos(cz,eta,3)
p(3,1)=DrcCos(cx,tao,3)
p(3,2)=DrcCos(cy,tao,3)
p(3,3)=DrcCos(cz,tao,3)

Return
End

!-----------------------------------------------------------------------------
Subroutine Transfmx_Dabao(angles,P) 
Implicit Double Precision(A-H, O-Z)
real*8 angles(4), P(3,3)
!Angle=/sin(phi),cos(phi),cos(theta), sin(theta)/
!theta: azimuth, the angle from x-axis to the projection of the normal vector to x-y plane
!phi:	zenith, the angle from z-axis to the normal vector	
ctht=angles(3); stht=angles(4); cphi=angles(2); sphi=angles(1);
 
p(1,1)=sphi*ctht ; p(1,2)=sphi*stht ; p(1,3)=cphi
p(2,1)=-stht      ; p(2,2)=ctht     ; p(2,3)=0.d0
p(3,1)=-cphi*ctht ; p(3,2)=-cphi*stht; p(3,3)=sphi

Return
End


!------------------------------------------------------------------
Subroutine Transfmx_Sining(xs,ys,zs,P)
implicit double precision(a-h,o-z)
real*8 p(3,3)
ctht=costheta(xs,ys,zs) 
stht=sintheta(xs,ys,zs)
cphi=cosphi(xs,ys,zs)
sphi=sinphi(xs,ys,zs)
print *,real(ctht),real(stht),real(cphi),real(sphi)
call coortrans(p,sphi,cphi,stht,ctht)
return
end

!-----------------------------------------------------------------------------
subroutine findjumps(xs,ys,zs,P,rjp)
implicit double precision(a-h,o-z)
real*8 rjp(1:4),p(3,3),du(0:1,1:3), dup(3),dun(3)

du(0,1)=udxn(xs,ys,zs);	du(0,2)=udyn(xs,ys,zs); du(0,3)=udzn(xs,ys,zs)
du(1,1)=udxp(xs,ys,zs); du(1,2)=udyp(xs,ys,zs);	du(1,3)=udzp(xs,ys,zs)

dun(1:3)=matmul(P,du(0,:));	dup(1:3)=matmul(P,du(1,:))

rjp(1)=uj(xs,ys,zs)
rjp(2)=betap(xs,ys,zs)*dup(1)-betan(xs,ys,zs)*dun(1)
rjp(3)=dup(2)-dun(2)
rjp(4)=dup(3)-dun(3)

end

!-----------------------------------------------------------------------------
subroutine findRegJumps(xs,ys,zs,P,rjp,eps0,ipm)
implicit double precision(a-h,o-z)
real*8 rjp(1:4),p(3,3)
real*8 phi_star_grad(3), phi_star_grad_m(3), phi_star_grad_d(3), phi_star_grad_q(3)

rjp(1)=phi_star(xs,ys,zs)
! print *, "ipm is ", ipm
if (ipm==1) then
	call phi_star_grad_mul(xs,ys,zs,phi_star_grad_d,phi_star_grad_q) !yang
	phi_star_grad_m=(/phis_x(xs,ys,zs),phis_y(xs,ys,zs),phis_z(xs,ys,zs)/)
	phi_star_grad=phi_star_grad_m + phi_star_grad_d + phi_star_grad_q
! 	phi_star_grad=phi_star_grad_d
elseif (ipm==0) then
	phi_star_grad=(/phis_x(xs,ys,zs),phis_y(xs,ys,zs),phis_z(xs,ys,zs)/)
endif

rjp(2)=eps0*dot_product(phi_star_grad,P(1,1:3))
rjp(3)=dot_product(phi_star_grad,P(2,1:3))
rjp(4)=dot_product(phi_star_grad,P(3,1:3))

end




!------------------------------------------------------------------------------
Function Det3(a)
implicit double precision(a-h, o-z)
Real*8 a(3,3), det
det3=a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))-a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1))+a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
End
!-----------------------------------------------------------------------------
Function DrcCos(A,B,n)
Implicit Double precision(A-H, O-Z)
Real*8 A(n), B(n), DrcCos
DrcCos=0.d0
do i=1,n
	DrcCos=DrcCos+A(i)*B(i)
enddo

Return
End
!--------------------------------------------------------------------------------------
subroutine findprtcts(i,j,k,p,betan,betap,cpp,cpn,cs,ct,const,ixtp,idn,rjp)
implicit double precision(a-h,o-z)
real*8 rjp(4),p(3,3),b(3,2),c(0:5), d(2:4,5), beta(2),cc(5),dd(4),t(2,4) !d(1-5): p+, p-, s, t, contant term
integer ixtp(2)
cpp=0.d0;	cpn=0.d0;	cs=0.d0;	ct=0.d0;	const=0.d0;		
d(2:4,5)=rjp(2:4); 


if ((abs(p(1,1))<1.d-10) .or. (abs(p(1,2))<1.d-10) .or. (abs(p(1,3))<1.d-10)) then
	its=0
	if (abs(p(1,1))<1.d-10) then			!cth=0, y-z plane     
		t(1,4)=rjp(2);	t(2,4)=rjp(4);
		if (idn==3) then
			its=2
			if (ixtp(2)==1) then			!Cancel px, pz-
				beta1=betan;	beta2=betap 
			else
				beta1=betap;	beta2=-betan	!Cancel px, pz+
			endif
			pp=p(1,3);		qq=p(3,3); 
			t(1,1:3)=(/betap*p(1,2),-betan*p(1,2), beta2*p(1,3)/)
			t(2,1:3)=(/p(3,2),-p(3,2), ixtp(2)*p(3,3)/)
		else	!idn=5
			its=2
			if (ixtp(2)==1) then			!Cancel px, py-
				beta1=betan;	beta2=betap
			else
				beta1=betap;	beta2=-betan	!Cancel px, py+
			endif
			pp=p(1,2);		qq=p(3,2);
			t(1,1:3)=(/betap*p(1,3),-betan*p(1,3), beta2*p(1,2)/)
			t(2,1:3)=(/p(3,3),-p(3,3), ixtp(2)*p(3,2)/)
		endif

	elseif (abs(p(1,2))<1.d-10) then		!sth=0, x-z plane
		t(1,4)=rjp(2);	t(2,4)=rjp(4);
		if (idn==1) then
			its=2
			if (ixtp(2)==1) then			!Cancel py, pz-
				beta1=betan;	beta2=betap 
			else
				beta1=betap;	beta2=-betan	!Cancel px, pz+
			endif
			pp=p(1,3);		qq=p(3,3);
			t(1,1:3)=(/betap*p(1,1),-betan*p(1,1), beta2*p(1,3)/)
			t(2,1:3)=(/p(3,1),-p(3,1), ixtp(2)*p(3,3)/)

		else	!idn=5
			
			its=1
			if (ixtp(1)==1) then			!Cancel py, px-
				beta1=betan;	beta2=betap 
			else
				beta1=betap;	beta2=-betan	!Cancel px, px+
			endif
			pp=p(1,1);		qq=p(3,1);
			t(1,1:3)=(/betap*p(1,3),-betan*p(1,3), beta2*p(1,1)/)
			t(2,1:3)=(/p(3,3),-p(3,3), ixtp(1)*p(3,1)/)
		endif

	else									!cphi=0, x-y plane
		t(1,4)=rjp(2);	t(2,4)=rjp(3);
		if (idn==1) then
			its=1
			if (ixtp(1)==1) then			!Cancel pz, py-
				beta1=betan;	beta2=betap 
			else
				beta1=betap;	beta2=-betan	!Cancel pz, py+
			endif
			pp=p(1,2);		qq=p(2,2);
			t(1,1:3)=(/betap*p(1,1),-betan*p(1,1), beta2*p(1,2)/)
			t(2,1:3)=(/p(2,1),-p(2,1), ixtp(1)*p(2,2)/)
		else	!idn=3
			its=1
			if (ixtp(1)==1) then			!Cancel pz, px-
				beta1=betan;	beta2=betap 
			else
				beta1=betap;	beta2=-betan	!Cancel pz, px+
			endif
			pp=p(1,1);		qq=p(2,1);
			t(1,1:3)=(/betap*p(1,2),-betan*p(1,2), beta2*p(1,1)/)
			t(2,1:3)=(/p(2,2),-p(2,2), ixtp(1)*p(2,1)/)
		endif
	endif

dd=t(2,:)*pp*beta1-t(1,:)*qq
cpp=dd(1);	cpn=dd(2);	ct=dd(3);	const=dd(4)

if (its==1) then
	cs=dd(3);	ct=0.d0
elseif (its==2) then
	cs=0.d0;	ct=dd(3)
else
	print *,'sth. is wrong!!'
endif

return
endif

select case (idn)
case(1)
	b(:,1)=p(:,2);	b(:,2)=p(:,3);
	c=(/p(1,2), p(2,2), p(3,2), p(1,3), p(2,3), p(3,3)/)
	if (ixtp(1)==1) then
		if (ixtp(2)==1) then 		!Cancel py-,pz-
			beta=(/betan,betan/)
			d(2,1:4)=(/betap*p(1,1),-betan*p(1,1), betap*p(1,2), betap*p(1,3)/)
			d(3,1:4)=(/p(2,1),-p(2,1), p(2,2), p(2,3)/)
			d(4,1:4)=(/p(3,1),-p(3,1), p(3,2), p(3,3)/)
		else						!Cancel py-,pz+
			beta=(/betan,betap/)
			d(2,1:4)=(/betap*p(1,1),-betan*p(1,1), betap*p(1,2), -betan*p(1,3)/)
			d(3,1:4)=(/p(2,1),-p(2,1), p(2,2), -p(2,3)/)
			d(4,1:4)=(/p(3,1),-p(3,1), p(3,2), -p(3,3)/)
		endif
	else 
		if (ixtp(2)==1) then		!Cancel py+,pz-
			beta=(/betap,betan/)
			d(2,1:4)=(/betap*p(1,1),-betan*p(1,1), -betan*p(1,2), betap*p(1,3)/)
			d(3,1:4)=(/p(2,1),-p(2,1), -p(2,2), p(2,3)/)
			d(4,1:4)=(/p(3,1),-p(3,1), -p(3,2), p(3,3)/)
		else						!Cancel py+,pz+
			beta=(/betap,betap/)
			d(2,1:4)=(/betap*p(1,1),-betan*p(1,1), -betan*p(1,2), -betan*p(1,3)/)
			d(3,1:4)=(/p(2,1),-p(2,1), -p(2,2), -p(2,3)/)
			d(4,1:4)=(/p(3,1),-p(3,1), -p(3,2), -p(3,3)/)
		endif
	endif
case(3)
	b(:,1)=p(:,1);	b(:,2)=p(:,3);
	c=(/p(1,1), p(2,1), p(3,1), p(1,3), p(2,3), p(3,3)/)
	if (ixtp(1)==1) then
		if (ixtp(2)==1) then 		!Cancel px-,pz-
			beta=(/betan,betan/)
			d(2,1:4)=(/betap*p(1,2),-betan*p(1,2), betap*p(1,1), betap*p(1,3)/)
			d(3,1:4)=(/p(2,2),-p(2,2), p(2,1), p(2,3)/)
			d(4,1:4)=(/p(3,2),-p(3,2), p(3,1), p(3,3)/)
		else						!Cancel px-,pz+
			beta=(/betan,betap/)
			d(2,1:4)=(/betap*p(1,2),-betan*p(1,2), betap*p(1,1), -betan*p(1,3)/)
			d(3,1:4)=(/p(2,2),-p(2,2), p(2,1), -p(2,3)/)
			d(4,1:4)=(/p(3,2),-p(3,2), p(3,1), -p(3,3)/)
		endif
	else 
		if (ixtp(2)==1) then		!Cancel px+,pz-
			beta=(/betap,betan/)
			d(2,1:4)=(/betap*p(1,2),-betan*p(1,2), -betan*p(1,1), betap*p(1,3)/)
			d(3,1:4)=(/p(2,2),-p(2,2), -p(2,1), p(2,3)/)
			d(4,1:4)=(/p(3,2),-p(3,2), -p(3,1), p(3,3)/)
		else						!Cancel px+,pz+
			beta=(/betap,betap/)
			d(2,1:4)=(/betap*p(1,2),-betan*p(1,2), -betan*p(1,1), -betan*p(1,3)/)
			d(3,1:4)=(/p(2,2),-p(2,2), -p(2,1), -p(2,3)/)
			d(4,1:4)=(/p(3,2),-p(3,2), -p(3,1), -p(3,3)/)
		endif
	endif
case(5)
	b(:,1)=p(:,1);	b(:,2)=p(:,2);
	c=(/p(1,1), p(2,1), p(3,1), p(1,2), p(2,2), p(3,2)/)
	if (ixtp(1)==1) then
		if (ixtp(2)==1) then 		!Cancel px-,py-
			beta=(/betan,betan/)
			d(2,1:4)=(/betap*p(1,3),-betan*p(1,3), betap*p(1,1), betap*p(1,2)/)
			d(3,1:4)=(/p(2,3),-p(2,3), p(2,1), p(2,2)/)
			d(4,1:4)=(/p(3,3),-p(3,3), p(3,1), p(3,2)/)
		else						!Cancel px-,py+
			beta=(/betan,betap/)
			d(2,1:4)=(/betap*p(1,3),-betan*p(1,3), betap*p(1,1), -betan*p(1,2)/)
			d(3,1:4)=(/p(2,3),-p(2,3), p(2,1), -p(2,2)/)
			d(4,1:4)=(/p(3,3),-p(3,3), p(3,1), -p(3,2)/)
		endif
	else 
		if (ixtp(2)==1) then		!Cancel px+,py-
			beta=(/betap,betan/)
			d(2,1:4)=(/betap*p(1,3),-betan*p(1,3), -betan*p(1,1), betap*p(1,2)/)
			d(3,1:4)=(/p(2,3),-p(2,3), -p(2,1), p(2,2)/)
			d(4,1:4)=(/p(3,3),-p(3,3), -p(3,1), p(3,2)/)
		else						!Cancel px+,py+
			beta=(/betap,betap/)
			d(2,1:4)=(/betap*p(1,3),-betan*p(1,3), -betan*p(1,1), -betan*p(1,2)/)
			d(3,1:4)=(/p(2,3),-p(2,3), -p(2,1), -p(2,2)/)
			d(4,1:4)=(/p(3,3),-p(3,3), -p(3,1), -p(3,2)/)
		endif
	endif
end select
if (abs(b(2,1))<1.d-10) then
	cc	=(d(2,:)-beta(1)*b(1,1)/b(3,1)*d(4,:)) &
		-(d(3,:)-b(2,1)/b(3,1)*d(4,:)) * (beta(2)*b(1,2)-beta(1)*b(3,2)*b(1,1)/b(3,1)) / (b(2,2)-b(3,2)*b(2,1)/b(3,1))
else
	cc	=(d(2,:)-beta(1)*b(1,1)/b(2,1)*d(3,:)) &
		-(d(4,:)-b(3,1)/b(2,1)*d(3,:)) * (beta(2)*b(1,2)-beta(1)*b(2,2)*b(1,1)/b(2,1)) / (b(3,2)-b(2,2)*b(3,1)/b(2,1))
endif
cpp=cc(1);	cpn=cc(2);	cs=cc(3);	ct=cc(4);	const=cc(5)
end

!!------------------------------------------------------------------------------
subroutine intpraxpts(it,jt,kt, xs,ys, zs, ix1, ix2, intppts,iside,ijust)
use comdata
use pbedata
use molecule
implicit double precision(a-h, o-z)
integer ix1(3), ix2(3), intppts(2,6,3), ia(4,4,3), ns(4), np(4), ms(4), mp(4), ib(2,2), idn, ipn
real*8 xx(4),yy(4),zz(4)

if (ix1(1) .ne. ix2(1)) then
	idn=1	!x direction
elseif (ix1(2) .ne. ix2(2)) then
	idn=3	!y direction
else
	idn=5	!z direction
endif 

!mp=(/-1,0,-2,1/);	
mp=(/0,-1,1,-2/);

ms=(/-1,-2,1,2/)		

i=ix1(1);	j=ix1(2);	k=ix1(3)

ijust1=0;	ijust2=0
do ipn=idn,idn+1
	select case (ipn)
	case(1)		!X direction X-Y plain
		xx=xs;	yy=y(j+ms);	zz=zs		!(xx(1:4),yy(1:4),zz(1:4)) are coodinates of points to be intrapolated 
		call find6pts(it,jt,kt,i,j,k,xx,yy,zz,ipn,ib,mp,ms,iside,ijust1)
		do ijk=0,2
			intppts(1,1+ijk,:)=(/i+ib(1,1)+ijk,j+ib(1,2),k/)
			intppts(1,4+ijk,:)=(/i+ib(2,1)+ijk,j+ib(2,2),k/)
		enddo
	case(2)		!X direction X-Z plain
		xx=xs;	yy=ys;	zz=z(k+ms)
		call find6pts(it,jt,kt,i,j,k,xx,yy,zz,ipn,ib,mp,ms,iside,ijust2)
		do ijk=0,2
			intppts(2,1+ijk,:)=(/i+ib(1,1)+ijk,j,k+ib(1,2)/)
			intppts(2,4+ijk,:)=(/i+ib(2,1)+ijk,j,k+ib(2,2)/)
		enddo
	case(3)		!Y direction X-Y plain
		xx=x(i+ms);	yy=ys;	zz=zs
		call find6pts(it,jt,kt,i,j,k,xx,yy,zz,ipn,ib,mp,ms,iside,ijust1)
		do ijk=0,2
			intppts(1,1+ijk,:)=(/i+ib(1,2),j+ib(1,1)+ijk,k/)
			intppts(1,4+ijk,:)=(/i+ib(2,2),j+ib(2,1)+ijk,k/)
		enddo
	case(4)		!Y direction Y-Z plain
		xx=xs;	yy=ys;	zz=z(k+ms)
		call find6pts(it,jt,kt,i,j,k,xx,yy,zz,ipn,ib,mp,ms,iside,ijust2)
		do ijk=0,2
			intppts(2,1+ijk,:)=(/i,j+ib(1,1)+ijk,k+ib(1,2)/)
			intppts(2,4+ijk,:)=(/i,j+ib(2,1)+ijk,k+ib(2,2)/)
		enddo
	case(5)		!Z direction X-Z plain
		xx=x(i+ms);	yy=ys;	zz=zs
		call find6pts(it,jt,kt,i,j,k,xx,yy,zz,ipn,ib,mp,ms,iside,ijust1)
		do ijk=0,2
			intppts(1,1+ijk,:)=(/i+ib(1,2),j,k+ib(1,1)+ijk/)
			intppts(1,4+ijk,:)=(/i+ib(2,2),j,k+ib(2,1)+ijk/)
		enddo
	case(6)		!Z direction Y-Z plain
		xx=xs;	yy=y(j+ms);	zz=zs
		call find6pts(it,jt,kt,i,j,k,xx,yy,zz,ipn,ib,mp,ms,iside,ijust2)
		do ijk=0,2
			intppts(2,1+ijk,:)=(/i,j+ib(1,2),k+ib(1,1)+ijk/)
			intppts(2,4+ijk,:)=(/i,j+ib(2,2),k+ib(2,1)+ijk/)
		enddo
	end select
	

enddo

!0: both failed; 1: both ok; 2: 1st ok,2nd failed; 3: 2nd ok,1st failed; 
if (ijust1==1 .and. ijust2==1) then
	ijust = 1
elseif (ijust1==0 .and. ijust2 ==0) then
	ijust = 0
elseif (ijust==1 .and. ijust2 ==0) then
	ijust = 2 
else
	ijust = 3
endif

End

!------------------------------------------------------------------------------
subroutine ifind3pts(iaa,iok)
implicit double precision(a-h, o-z)
integer iaa(3), iaxpt

if ((iaa(1) .ne. iaa(2)) .or. (iaa(1) .ne. iaa(3)) .or. (iaa(2) .ne. iaa(3))) then
	iok=0
else
	iok=iaa(1)
endif
End

!-------------------------------------------------------------------------------
subroutine find6pts(it,jt,kt,i,j,k,xx,yy,zz,ipn,ib,mp,ms,iside,ijust)
use comdata
use pbedata
use molecule
implicit double precision(a-h,o-z)		
integer ms(4), mp(4), ib(2,2), iok(4,4), iaxpt(4),ia(4,4,3), iflag(2), iaa(3)
real*8 xx(4),yy(4),zz(4)			!xx,yy,zz are coordinates of the points to be extraplated
									
!mp=(/-1,0,-2,1/);	ms=(/-1,-2,1,2/)		
!call ioaxpt(xx,yy,zz,iaxpt)		!Test the inside or outside of points to be extrapolated
do ns=1,4
	do np=1,4
		select case (ipn)
		case(1) !x primary, y secondary
			ia(ns,np,:)=io(i+mp(np):i+mp(np)+2,j+ms(ns),k)
		case(2)	!x primary, z secondary
			ia(ns,np,:)=io(i+mp(np):i+mp(np)+2,j,k+ms(ns))
		case(3)	!y primary, x secondary
			ia(ns,np,:)=io(i+ms(ns),j+mp(np):j+mp(np)+2,k)
		case(4)	!y primary, z secondary
			ia(ns,np,:)=io(i,j+mp(np):j+mp(np)+2,k+ms(ns))
		case(5)	!z primary, x secondary
			ia(ns,np,:)=io(i+ms(ns),j,k+mp(np):k+mp(np)+2)
		case(6)	!z primary, y secondary
			ia(ns,np,:)=io(i,j+ms(ns),k+mp(np):k+mp(np)+2)
		end select
		iaa=ia(ns,np,1:3)
		call ifind3pts(iaa,iok(ns,np))
	enddo
enddo

	do isgn=1,-1,-2
		iflag=0
		do ns=1,3,2
			do np=1,4 
				if (iok(ns,np)==isgn) then
					iflag((ns+1)/2)=isgn
					ib((ns+1)/2,:)=(/mp(np),ms(ns)/)
					goto 100
				endif
			enddo
	100		continue
		enddo
		if (iflag(1)==isgn .and. iflag(2)==isgn) then
			goto 400
		endif		

	!enddo
	!do isgn=1,-1,-2
		iflag=0
		do ns=1,2
			do np=1,4
				if (iok(ns,np)==isgn) then
					iflag(ns)=isgn
					ib(ns,:)=(/mp(np),ms(ns)/)				!ib(2.2)	|p1 s1| 
															!			|p2 s2|											
					goto 200
				endif
			enddo
	200		continue
		enddo
		if (iflag(1)==isgn .and. iflag(2)==isgn) goto 400
		
	!enddo
	!do isgn=1,-1,-2
		iflag=0	
		do ns=3,4
			do np=1,4
				if (iok(ns,np)==isgn) then
					iflag(ns-2)=isgn
					ib(ns-2,:)=(/mp(np),ms(ns)/)
					goto 300
				endif
			enddo
	300		continue
		enddo		
		if (iflag(1)==isgn .and. iflag(2)==isgn) goto 400

	enddo

400 continue	

if ((iflag(1) == 1) .and. (iflag(2) == 1)) then 
	ijust=1
	continue
elseif ((iflag(1) == -1) .and. (iflag(2) == -1)) then 
	ijust=1
	!print *, i,j,k,' fail to get all interpolation from outside'
else
	!print *, it,jt,kt,'Could not find all points to interpolate'
	ib(1,:)=(/0,1/); ib(2,:)=(/0,2/)
endif

End
!-------------------------------------------------------------------------------
subroutine findtgwts(idx,ik,i,j,k,xs,ys,zs,w0p,w0n,w1p,w1n, phi,psi)
use pbedata
use comdata
implicit double precision(a-h, o-z)
real*8 w0p(3), w0n(3), w1p(3), w1n(3)
real*8 ck1(6), ck2(6), c1(6), c2(6)
integer idx(3)

a11=w0n(3)
a12=-w0p(3)
a21=betan(xs,ys,zs)*w1n(3)
a22=-betap(xs,ys,zs)*w1p(3)

ck1=(/	-w0n(1),	-w0n(2), w0p(2),	w0p(1),	-1.d0,	0.d0/)
ck2=(/	-betan(xs,ys,zs)*w1n(1),	-betan(xs,ys,zs)*w1n(2),	betap(xs,ys,zs)*w1p(2),	&
		betap(xs,ys,zs)*w1p(1),	0.d0,					-1.d0/)

C1=	(a22*ck1-a12*ck2)	/	(a11*a22-a12*a21)			!f^-
C2=	(a11*ck2-a21*ck1)	/	(a11*a22-a12*a21)			!f^+

ftpts(indirr(i,j,k),7-ik,1:4)=c2(1:4)					!f^+
ftc(indirr(i,j,k),7-ik)=c2(5)*phi+c2(6)*psi
ftpts(indirr(idx(1),idx(2),idx(3)),ik,1:4)=c1(1:4)		!f^-
ftc(indirr(idx(1),idx(2),idx(3)),ik)=(c1(5)*phi+c1(6)*psi)

end
!-------------------------------------------------------------------------------
subroutine idx6(i,j,k,kk,idx)
integer idx(3)
select case (kk)
case (1)	
	idx=(/i,j,k-1/)
case (2)
	idx=(/i,j-1,k/)
case (3)
	idx=(/i-1,j,k/)
case (4)
	idx=(/i+1,j,k/)
case (5)
	idx=(/i,j+1,k/)
case (6)
	idx=(/i,j,k+1/)
end select

End
!------------------------------------------------------------------------------
subroutine find1wts(xx,yy,ss,wn,wp)
implicit double precision(a-h, o-z)
real*8 ss,c(0:2,0:1), xx(0:2), yy(0:2), wn(3),wp(3), x1(0:2),y1(0:2)
x1=xx
y1=yy
call weights(ss,x1,2,2,1,c);
wn=c(:,1)
call weights(ss,y1,2,2,1,c)
wp=c(:,1)
End
!-------------------------------------------------------------------------------
subroutine findwts(idx,ik,i,j,k,wts,cpp,cpn,cs,ct,const,h,p0p,p0n,p1p,p1n,s1,t1,phi,ixtp)
use pbedata
use comdata
implicit double precision(a-h, o-z)
real*8 wts(2,16), h(2,2,3), p0p(3), p0n(3), p1p(3), p1n(3), s1(3), t1(3)
real*8 ck1(9), ck2(9), cc(2,9)
real*8 cpn, cpp, cs, ct, const
integer idx(3),ik,i,j,k,ixtp(2)

wts=0.d0

ck1=(/	p0n(1),	p0n(2), -p0p(2),	-p0p(1),	0.d0,	0.d0,	0.d0,	0.d0,	phi/)

ck2=(/	-cpn*p1n(1)-cs*s1(1)*p0n(1)-ct*t1(1)*p0n(1),	-cpn*p1n(2)-cs*s1(1)*p0n(2)-ct*t1(1)*p0n(2),&
		-cpp*p1p(2),	-cpp*p1p(1),	-cs*s1(2),	-cs*s1(3),	-ct*t1(2),	-ct*t1(3),&
		const-phi*(cs*s1(1)*(2.d0*ixtp(1)+2.d0)/4.d0+ct*t1(1)*(2.d0*ixtp(2)+2.d0)/4.d0) /)                      

a11=p0p(3);			a12=-p0n(3)
a21=cpp*p1p(3);		a22=cpn*p1n(3)+(cs*s1(1)+ct*t1(1))*p0n(3)

cc(1,:)=	(a22*ck1-a12*ck2)	/	(a11*a22-a12*a21)			!f^+
cc(2,:)=	(a11*ck2-a21*ck1)	/	(a11*a22-a12*a21)			!f^-

wts(:,1)	=cc(:,1)
wts(:,2)	=cc(:,2)
wts(:,3)	=cc(:,3)	
wts(:,4)	=cc(:,4)	
wts(:,5)	=cc(:,5) * h(1,1,1)
wts(:,6)	=cc(:,5) * h(1,1,2)
wts(:,7)	=cc(:,5) * h(1,1,3)
wts(:,8)	=cc(:,6) * h(1,2,1)
wts(:,9)	=cc(:,6) * h(1,2,2)
wts(:,10)	=cc(:,6) * h(1,2,3)
wts(:,11)	=cc(:,7) * h(2,1,1)
wts(:,12)	=cc(:,7) * h(2,1,2)
wts(:,13)	=cc(:,7) * h(2,1,3)
wts(:,14)	=cc(:,8) * h(2,2,1)
wts(:,15)	=cc(:,8) * h(2,2,2)
wts(:,16)	=cc(:,8) * h(2,2,3)


ftpts(indirr(i,j,k),7-ik,1:16)=wts(1,:)									!f^+
ftc(indirr(i,j,k),7-ik)=cc(1,9)

ftpts(indirr(idx(1),idx(2),idx(3)),ik,1:16)=wts(2,:)						!f^-
ftc(indirr(idx(1),idx(2),idx(3)),ik)=cc(2,9)

end


!-------------------------------------------------------------------------------
subroutine grdptwts(i,j,k) 
use pbedata
use comdata
Implicit Double precision(A-H, O-Z)
real*8 wts(16), p(3,3), c(0:2,0:1), xx(0:2), pcts(4), ck1(14), ck2(14), c1(14), c2(14),flux(3)
real*8 wxp(3), wxn(3), wyp(3), wyn(3), wzp(3), wzn(3), rjp(4), aa(2,2), trig(4)
integer iwts(16,3), idx(3), idy(3), it(6), is(6), id(6), temp(3)
iwts=0; wts=0.d0
ijk=(k-1)*nx*ny+(j-1)*nx+i
xs=x(i); ys=y(j); zs=z(k)
btp=betap(xs,ys,zs);	btn=betan(xs,ys,zs)
irr=indirr(i,j,k)

if (isf==0) then
	call Transfmx(xs,ys,zs,p)
	!call Transfmx_Sining(xs,ys,zs,P)
	!ctht=costheta(xs,ys,zs) 
	!stht=sintheta(xs,ys,zs)
	!cphi=cosphi(xs,ys,zs)
	!sphi=sinphi(xs,ys,zs)

	!^^^^^^^^^^^^^^^^^^^^^
	!call transfmx_cube(xs,ys,zs,p)
else
	!trig=clocal(:,1,irr)
	trig=dist_angle(1:4,1,indirr(i,j,k))
	call Transfmx_dabao(trig,p)			! Use dabao's interface
	if (sum(abs(trig))<1.d-8) then
		print *,i,j,k, 'all zero normal condition'
		return
	endif
endif

if (icg==0) then
	call findjumps(xs,ys,zs,P,rjp)		! Use this if the true solution known
elseif (icg==2) then
        call findRegJumps(xs,ys,zs,P,rjp,eps0,ipm)
endif

! For outside point whose 7 points schemes touches the interface
do ii=1,6
	call idx6(i,j,k,ii,idx)
	if (io(idx(1),idx(2),idx(3))==1) then
		ix3=(idx(3)-1)*(nx*ny)+(idx(2)-1)*nx+idx(1)
		cof=betap((x(idx(1))+x(i))/2,(y(idx(2))+y(j))/2,(z(idx(3))+z(k))/2)
		cof=cof/(abs(x(idx(1))-x(i))+abs(y(idx(2))-y(j))+abs(z(idx(3))-z(k)))**2
		bftc(ix3)=bftc(ix3)-cof*rjp(1)
	endif
enddo

ismgrds=io(i+1,j,k)+io(i-1,j,k)+io(i,j+1,k)+io(i,j-1,k)+io(i,j,k+1)+io(i,j,k-1)	


if (io(i+1,j,k)==0) ismgrds=ismgrds+1
if (io(i-1,j,k)==0) ismgrds=ismgrds+1
if (io(i,j+1,k)==0) ismgrds=ismgrds+1
if (io(i,j-1,k)==0) ismgrds=ismgrds+1
if (io(i,j,k+1)==0) ismgrds=ismgrds+1
if (io(i,j,k-1)==0)	ismgrds=ismgrds+1
!print *,i*10000+j*100+k
!print *,io(i-1:i+1,j,k)
!print *,io(i,j-1:j+1,k)
!print *,io(i,j,k-1:k+1)
if (ismgrds==4 .or. ismgrds==-4) then   !might have problems here
	isn=ismgrds/4						!isn postive if 4+ 2- and negative if 4- 2+
	call grd4p4n(isn,i,j,k,xs,ys,zs,rjp)
elseif (ismgrds==2 .or. ismgrds==-2) then
	isn=ismgrds/2
	!print *,p
	!print *,xs,ys,zs
	!write(*,*)
	call grd2p2n(isn,i,j,k,xs,ys,zs,rjp,btn,btp,p)
else			!three positve and three negative
	call grd3p3n(i,j,k,xs,ys,zs,iwts,rjp,btn,btp,p)


endif
End
!-------------------------------------------------------------------------------
subroutine grd2p2n(isn,i,j,k,xs,ys,zs,rjp,btn,btp,p)
use comdata
use pbedata
implicit double precision(a-h,o-z)	
real*8 wxp(3), wxn(3), wyp(3), wyn(3), wzp(3), wzn(3), rjp(4), cx(2,2), cy(2,2), cz(2,2)
real*8 ck1(14), ck2(14), c1(14), c2(14), aa(2,2), pp(3), p(3,3)
integer it(6), id(6), iwts(16,3), idx(3), idy(3)
	ijk=(k-1)*nx*ny+(j-1)*nx+i
	isrf(0,indirr(i,j,k))=isn
	it=0
	iwts(1:4,2)=j;		iwts(1:4,3)=k
	iwts(5:7,1)=i;		iwts(5:7,3)=k
	iwts(8:10,1)=i;		iwts(8:10,2)=j

	do ii=1,6
		call idx6(i,j,k,ii,idx)
		call idx6(i,j,k,7-ii,idy)
	enddo
	a1=p(2,1)*p(2,2)+p(2,2)*p(2,3)+p(2,1)*p(2,3)
	a2=p(3,1)*p(3,2)+p(3,2)*p(3,3)+p(3,1)*p(3,3)
	flag=a1*a2
	call findgdwts(-isn,i,j,k,xs,ys,zs,wxp,wxn,wyp,wyn,wzp,wzn,iwts,it,id)
	!print *,i*10000+j*100+k,flag
	if (abs(flag)<1.d-10) then							!Tangent
		if (isn==1) then							!4 positive 2 negative
			
			if (abs(p(1,1))<1.d-10) then
				aa(1,1)=-btp*p(1,2)*wyp(3);	aa(1,2)=-btp*p(1,3)*wzp(3)
				aa(2,1)=-p(3,2)*wyp(3);		aa(2,2)=-p(3,3)*wzp(3)
				
				ck1=(/	0.d0, btp*(p(1,2)*wyp(2)+p(1,3)*wzp(2))-btn*(p(1,2)*wyn(1)+p(1,3)*wzn(1)), &
						0.d0, 0.d0, btp*p(1,2)*wyp(1), -btn*p(1,2)*wyn(2), -btn*p(1,2)*wyn(3), &
						btp*p(1,3)*wzp(1), -btn*p(1,3)*wzn(2), -btn*p(1,3)*wzn(3), &
						btp*(p(1,2)*wyp(2)+p(1,3)*wzp(2)), -1.d0, 0.d0, 0.d0/)
				
				ck2=(/	0.d0, (p(3,2)*wyp(2)+p(3,3)*wzp(2))-(p(3,2)*wyn(1)+p(3,3)*wzn(1)), &
						0.d0, 0.d0, p(3,2)*wyp(1), -p(3,2)*wyn(2), -p(3,2)*wyn(3), &
						p(3,3)*wzp(1), -p(3,3)*wzn(2), -p(3,3)*wzn(3), &
						(p(3,2)*wyp(2)+p(3,3)*wzp(2)), 0.d0, 0.d0, -1.d0/)				
				
			elseif(abs(p(1,2))<1.d-10) then
				aa(1,1)=-btp*p(1,1)*wxp(3);	aa(1,2)=-btp*p(1,3)*wzp(3)
				aa(2,1)=-p(3,1)*wxp(3);		aa(2,2)=-p(3,3)*wzp(3)
				
				ck1=(/	btp*p(1,1)*wxp(1), btp*(p(1,1)*wxp(2)+p(1,3)*wzp(2))-btn*(p(1,1)*wxn(1)+p(1,3)*wzn(1)), &
						-btn*p(1,1)*wxn(2), -btn*p(1,1)*wxn(3), 0.d0, 0.d0, 0.d0, &
						btp*p(1,3)*wzp(1), -btn*p(1,3)*wzn(2), -btn*p(1,3)*wzn(3), &
						btp*(p(1,1)*wxp(2)+p(1,3)*wzp(2)), -1.d0, 0.d0, 0.d0/)
				
				ck2=(/	p(3,1)*wxp(1), (p(3,1)*wxp(2)+p(3,3)*wzp(2))-(p(3,1)*wxn(1)+p(3,3)*wzn(1)), &
						-p(3,1)*wxn(2), -p(3,1)*wxn(3), 0.d0, 0.d0, 0.d0, &
						p(3,3)*wzp(1), -p(3,3)*wzn(2), -p(3,3)*wzn(3), &
						(p(3,1)*wxp(2)+p(3,3)*wzp(2)), 0.d0, 0.d0, -1.d0/)		
			else
				aa(1,1)=-btp*p(1,1)*wxp(3);	aa(1,2)=-btp*p(1,2)*wyp(3)
				aa(2,1)=-p(2,1)*wxp(3);		aa(2,2)=-p(2,2)*wyp(3)
				ck1=(/	btp*p(1,1)*wxp(1), btp*(p(1,1)*wxp(2)+p(1,2)*wyp(2))-btn*(p(1,1)*wxn(1)+p(1,2)*wyn(1)), &
						-btn*p(1,1)*wxn(2), -btn*p(1,1)*wxn(3), btp*p(1,2)*wyp(1), -btn*p(1,2)*wyn(2), -btn*p(1,2)*wyn(3), &
						0.d0, 0.d0, 0.d0, &
						btp*(p(1,1)*wxp(2)+p(1,2)*wyp(2)), -1.d0, 0.d0, 0.d0/)
				
				ck2=(/	p(2,1)*wxp(1), (p(2,1)*wxp(2)+p(2,2)*wyp(2))-(p(2,1)*wxn(1)+p(2,2)*wyn(1)), &
						-p(2,1)*wxn(2), -p(2,1)*wxn(3), p(2,2)*wyp(1), -p(2,2)*wyn(2), -p(2,2)*wyn(3), &
						0.d0, 0.d0, 0.d0, &
						(p(2,1)*wxp(2)+p(2,2)*wyp(2)), 0.d0, -1.d0, 0.d0/)			

			endif
			bftc(ijk)=	bftc(ijk) + fj(xs,ys,zs) &
					+(betap(xs,ys+dy/2.d0,zs)+betap(xs,ys-dy/2.d0,zs)+betap(xs+dx/2.d0,ys,zs) &
					+betap(xs-dx/2.d0,ys,zs)+betap(xs,ys,zs-dz/2.d0)+betap(xs,ys,zs+dz/2.d0))/dcel**2*rjp(1) &
					-cappap(xs,ys,zs)*rjp(1)

		
		else										!2 positive 4 negative
			if (abs(p(1,1))<1.d-10) then
				aa(1,1)=btn*p(1,2)*wyn(3);	aa(1,2)=btn*p(1,3)*wzn(3)
				aa(2,1)=p(3,2)*wyn(3);		aa(2,2)=p(3,3)*wzn(3)
					
				ck1=(/	0.d0, btp*(p(1,2)*wyp(1)+p(1,3)*wzp(1))-btn*(p(1,2)*wyn(2)+p(1,3)*wzn(2)), &
						0.d0, 0.d0, -btn*p(1,2)*wyn(1), btp*p(1,2)*wyp(2), btp*p(1,2)*wyp(3), &
						-btn*p(1,3)*wzn(1), btp*p(1,3)*wzp(2), btp*p(1,3)*wzp(3), &
						btp*(p(1,2)*wyp(1)+p(1,3)*wzp(1)), -1.d0, 0.d0, 0.d0/)
				
				ck2=(/	0.d0, (p(3,2)*wyp(1)+p(3,3)*wzp(1))-(p(3,2)*wyn(2)+p(3,3)*wzn(2)), &
						0.d0, 0.d0, -p(3,2)*wyn(1), p(3,2)*wyp(2), p(3,2)*wyp(3), &
						-p(3,3)*wzn(1), p(3,3)*wzp(2), p(3,3)*wzp(3), &
						(p(3,2)*wyp(1)+p(3,3)*wzp(1)), 0.d0, 0.d0, -1.d0/)				
					
			elseif(abs(p(1,2))<1.d-10) then
				aa(1,1)=btn*p(1,1)*wxn(3);	aa(1,2)=btn*p(1,3)*wzn(3)
				aa(2,1)=p(3,1)*wxn(3);		aa(2,2)=p(3,3)*wzn(3)
					
				ck1=(/	-btn*p(1,1)*wxn(1), btp*(p(1,1)*wxp(1)+p(1,3)*wzp(1))-btn*(p(1,1)*wxn(2)+p(1,3)*wzn(2)), &
						btp*p(1,1)*wxp(2), btp*p(1,1)*wxp(3), 0.d0, 0.d0, 0.d0, &
						-btn*p(1,3)*wzn(1), btp*p(1,3)*wzp(2), btp*p(1,3)*wzp(3), &
						btp*(p(1,1)*wxp(1)+p(1,3)*wzp(1)), -1.d0, 0.d0, 0.d0/)
				
				ck2=(/	-p(3,1)*wxn(1), (p(3,1)*wxp(1)+p(3,3)*wzp(1))-(p(3,1)*wxn(2)+p(3,3)*wzn(2)), &
						p(3,1)*wxp(2), p(3,1)*wxp(3), 0.d0, 0.d0, 0.d0, &
						-p(3,3)*wzn(1), p(3,3)*wzp(2), p(3,3)*wzp(3), &
						(p(3,1)*wxp(1)+p(3,3)*wzp(1)), 0.d0, 0.d0, -1.d0/)		
			else
				aa(1,1)=btn*p(1,1)*wxn(3);	aa(1,2)=btn*p(1,2)*wyn(3)
				aa(2,1)=p(2,1)*wxn(3);		aa(2,2)=p(2,2)*wyn(3)
				ck1=(/	-btn*p(1,1)*wxn(1), btp*(p(1,1)*wxp(1)+p(1,2)*wyp(1))-btn*(p(1,1)*wxn(2)+p(1,2)*wyn(2)), &
						btp*p(1,1)*wxp(2), btp*p(1,1)*wxp(3), -btn*p(1,2)*wyn(1), btp*p(1,2)*wyp(2), btp*p(1,2)*wyp(3), &
						0.d0, 0.d0, 0.d0, &
						btp*(p(1,1)*wxp(1)+p(1,2)*wyp(1)), -1.d0, 0.d0, 0.d0/)
					
				ck2=(/	-p(2,1)*wxn(1), (p(2,1)*wxp(1)+p(2,2)*wyp(1))-(p(2,1)*wxn(2)+p(2,2)*wyn(2)), &
						p(2,1)*wxp(2), p(2,1)*wxp(3), -p(2,2)*wyn(1), p(2,2)*wyp(2), p(2,2)*wyp(3), &
						0.d0, 0.d0, 0.d0, &
						(p(2,1)*wxp(1)+p(2,2)*wyp(1)), 0.d0, -1.d0, 0.d0/)			
			endif
		endif

		dta=(aa(1,1)*aa(2,2)-aa(1,2)*aa(2,1))
		!print *,'llf',dta
		c1=(ck1*aa(2,2)-ck2*aa(1,2))/dta
		c2=(ck2*aa(1,1)-ck1*aa(2,1))/dta
		do l=1,6
			call idx6(i,j,k,l,idx)
			if (it(l)==1 .and. id(l) .ne. 1) then
				
				indsrf=indsrf+1
				isrf(l,indirr(i,j,k))=indsrf
				isrfpts(indsrf,1:28:3)=iwts(1:10,1)
				isrfpts(indsrf,2:29:3)=iwts(1:10,2)
				isrfpts(indsrf,3:30:3)=iwts(1:10,3)
					
				if (abs(p(1,1))<1.d-10) then
					if (l==1 .or. l==6) then
						srfpts(indsrf,1:10)=C2(1:10)
						srfftc(indsrf)=sum(c2(11:14)*rjp(1:4))
					else
						srfpts(indsrf,1:10)=C1(1:10)
						srfftc(indsrf)=sum(c1(11:14)*rjp(1:4))
					endif
				elseif(abs(p(1,2))<1.d-10) then
					if (l==1 .or. l==6) then
						srfpts(indsrf,1:10)=C2(1:10)
						srfftc(indsrf)=sum(c2(11:14)*rjp(1:4))
					else
						srfpts(indsrf,1:10)=C1(1:10)
						srfftc(indsrf)=sum(c1(11:14)*rjp(1:4))
					endif
				else
					if (l==2 .or. l==5) then
						srfpts(indsrf,1:10)=C2(1:10)
						srfftc(indsrf)=sum(c2(11:14)*rjp(1:4))
					else
						srfpts(indsrf,1:10)=C1(1:10)
						srfftc(indsrf)=sum(c1(11:14)*rjp(1:4))
					endif
				endif

			endif
		enddo
	else												!Not Tangent
		if (isn==1) then								!4 positive 2 negative
			do ii=1,3
				if (id(ii)==1) then !pp(2)=p(2,ii)

					pp(:)=p(:,ii)
					cx(1,1)=btp*p(1,1)*pp(2)-btn*p(2,1)*pp(1);	cx(1,2)=-(btn*p(1,1)*pp(2)-btn*p(2,1)*pp(1))
					cy(1,1)=btp*p(1,2)*pp(2)-btn*p(2,2)*pp(1);	cy(1,2)=-(btn*p(1,2)*pp(2)-btn*p(2,2)*pp(1))
					cz(1,1)=btp*p(1,3)*pp(2)-btn*p(2,1)*pp(1);	cz(1,2)=-(btn*p(1,3)*pp(2)-btn*p(2,3)*pp(1))
					cx(2,1)=btp*p(1,1)*pp(3)-btn*p(3,1)*pp(1);	cx(2,2)=-(btn*p(1,1)*pp(3)-btn*p(3,1)*pp(1))
					cy(2,1)=btp*p(1,2)*pp(3)-btn*p(3,2)*pp(1);	cy(2,2)=-(btn*p(1,2)*pp(3)-btn*p(3,2)*pp(1))
					cz(2,1)=btp*p(1,3)*pp(3)-btn*p(3,1)*pp(1);	cz(2,2)=-(btn*p(1,3)*pp(3)-btn*p(3,3)*pp(1))

					ck1=(/	cx(1,2)*wxn(1),	&
							cx(1,1)*wxp(1)+cy(1,1)*wyp(1)+cz(1,1)*wzp(1)+cx(1,2)*wxn(2)+cy(1,2)*wyn(2)+cz(1,2)*wzn(2), &
							cx(1,2)*wxp(2),	cx(1,2)*wxp(3),	cy(1,2)*wyn(1), cy(1,1)*wyp(2), cy(1,1)*wyp(3), cz(1,2)*wzn(3), &
							cz(1,1)*wzp(2), cz(1,1)*wzp(3), -(cx(1,1)*wxp(1)+cy(1,1)*wyp(1)+cz(1,1)*wzp(1)), -pp(2), pp(1), 0.d0/)
					ck2=(/	cx(2,2)*wxn(1),	&
							cx(2,1)*wxp(1)+cy(2,1)*wyp(1)+cz(2,1)*wzp(1)+cx(2,2)*wxn(2)+cy(2,2)*wyn(2)+cz(2,2)*wzn(2), &
							cx(2,2)*wxp(2),	cx(2,2)*wxp(3),	cy(2,2)*wyn(1), cy(2,1)*wyp(2), cy(2,1)*wyp(3), cz(2,2)*wzn(3), &
							cz(2,1)*wzp(2), cz(2,1)*wzp(3), -(cx(2,1)*wxp(1)+cy(2,1)*wyp(1)+cz(2,1)*wzp(1)), -pp(3), 0.d0, pp(1)/)
				
					if (id(1)==1) then
						aa(1,1)=-cy(1,2)*wyn(3);		aa(1,2)=-cz(1,2)*wzn(3)
						aa(2,1)=-cy(2,2)*wyn(3);		aa(1,2)=-cz(2,2)*wzn(3)
					elseif (id(2)==1) then
						aa(1,1)=-cx(1,2)*wxn(3);		aa(1,2)=-cz(1,2)*wzn(3)
						aa(2,1)=-cx(2,2)*wxn(3);		aa(1,2)=-cz(2,2)*wzn(3)
					else
						aa(1,1)=-cx(1,2)*wxn(3);		aa(1,2)=-cy(1,2)*wyn(3)
						aa(2,1)=-cx(2,2)*wxn(3);		aa(1,2)=-cy(2,2)*wyn(3)
					endif
					
					dta=(aa(1,1)*aa(2,2)-aa(1,2)*aa(2,1))
					if (abs(dta)<1.d-8) return
					c1=(ck1*aa(2,2)-ck2*aa(1,2))/dta
					c2=(ck2*aa(1,1)-ck1*aa(2,1))/dta
					
					do l=1,6
						call idx6(i,j,k,l,idx)
						if (it(l)==1 .and. id(l) .ne. 1) then
							indsrf=indsrf+1
							isrf(l,indirr(i,j,k))=indsrf
							
							isrfpts(indsrf,1:28:3)=iwts(1:10,1)
							isrfpts(indsrf,2:29:3)=iwts(1:10,2)
							isrfpts(indsrf,3:30:3)=iwts(1:10,3)
							if (id(1)==1) then
								if (l==1 .or. l==6) then
									srfpts(indsrf,1:10)=C2(1:10)
									srfftc(indsrf)=sum(c2(11:14)*rjp(1:4))
								else
									srfpts(indsrf,1:10)=C1(1:10)
									srfftc(indsrf)=sum(c1(11:14)*rjp(1:4))
								endif
							elseif(id(2)==1) then
								if (l==1 .or. l==6) then
									srfpts(indsrf,1:10)=C2(1:10)
									srfftc(indsrf)=sum(c2(11:14)*rjp(1:4))
								else
									srfpts(indsrf,1:10)=C1(1:10)
									srfftc(indsrf)=sum(c1(11:14)*rjp(1:4))
								endif
							else
								if (l==2 .or. l==5) then
									srfpts(indsrf,1:10)=C2(1:10)
									srfftc(indsrf)=sum(c2(11:14)*rjp(1:4))
								else
									srfpts(indsrf,1:10)=C1(1:10)
									srfftc(indsrf)=sum(c1(11:14)*rjp(1:4))
								endif
							endif
						endif
					enddo
				endif
			enddo
		else
			print *,'temperarily not availle'
		endif
	endif
End
!-------------------------------------------------------------------------------
subroutine grd3p3n(i,j,k,xs,ys,zs,iwts,rjp,btn,btp,p)
use comdata
use pbedata
implicit double precision(a-h,o-z)
real*8 wxn(3),wxp(3),wyn(3),wyp(3),wzn(3),wzp(3),rjp(4),p(3,3),a(3,3),ck1(14),ck2(14),ck3(14),c1(14),c2(14),c3(14)
integer it(6),id(6),iwts(16,3),idx(3)

	c1=0.d0;	c2=0.d0;	c3=0.d0
	m=1;	it=0;
	iwts(1:4,2)=j;		iwts(1:4,3)=k
	iwts(5:7,1)=i;		iwts(5:7,3)=k
	iwts(8:10,1)=i;		iwts(8:10,2)=j

	call findgdwts(1,i,j,k,xs,ys,zs,wxn,wxp,wyn,wyp,wzn,wzp,iwts,it,id)	! Note: reverse wn and wp

	a(1,1)=btn*wxn(3)*p(1,1);	a(1,2)=btn*wyn(3)*p(1,2);	a(1,3)=btn*wzn(3)*p(1,3)
	a(2,1)=wxn(3)*p(2,1);		a(2,2)=wyn(3)*p(2,2);		a(2,3)=wzn(3)*p(2,3)
	a(3,1)=wxn(3)*p(3,1);		a(3,2)=wyn(3)*p(3,2);		a(3,3)=wzn(3)*p(3,3)
	dta=a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))-a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1))+a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
	print *,'lsf 3',dta
	if (dta<1.d-10) return
	ck1=(/	-btn*wxn(1)*p(1,1),	&
			 btp*(wxp(1)*p(1,1)+wyp(1)*p(1,2)+wzp(1)*p(1,3))-btn*(wxn(2)*p(1,1)+wyn(2)*p(1,2)+wzn(2)*p(1,3)), &
			 btp*wxp(2)*p(1,1),	 btp*wxp(3)*p(1,1),	&
			-btn*wyn(1)*p(1,2),	 btp*wyp(2)*p(1,2),	&
			 btp*wyp(3)*p(1,2),	-btn*wzn(1)*p(1,3),	&
			 btp*wzp(2)*p(1,3),	 btp*wzp(3)*p(1,3),	&
			 btp*(wxp(1)*p(1,1)+wyp(1)*p(1,2)+wzp(1)*p(1,3)), -1.d0, 0.d0, 0.d0 /)

	ck2=(/	-wxn(1)*p(2,1),	&
			 (wxp(1)*p(2,1)+wyp(1)*p(2,2)+wzp(1)*p(2,3))-(wxn(2)*p(2,1)+wyn(2)*p(2,2)+wzn(2)*p(2,3)), &
			 wxp(2)*p(2,1),	 wxp(3)*p(2,1),	&
			-wyn(1)*p(2,2),	 wyp(2)*p(2,2),	&
			 wyp(3)*p(2,2),	-wzn(1)*p(2,3),	&
			 wzp(2)*p(2,3),	 wzp(3)*p(2,3),	&
			 (wxp(1)*p(2,1)+wyp(1)*p(2,2)+wzp(1)*p(2,3)), 0.d0, -1.d0, 0.d0 /)

	ck3=(/	-wxn(1)*p(3,1),	&
			 (wxp(1)*p(3,1)+wyp(1)*p(3,2)+wzp(1)*p(3,3))-(wxn(2)*p(3,1)+wyn(2)*p(3,2)+wzn(2)*p(3,3)), &
			 wxp(2)*p(3,1),	 wxp(3)*p(3,1),	&
			-wyn(1)*p(3,2),	 wyp(2)*p(3,2),	&
			 wyp(3)*p(3,2),	-wzn(1)*p(3,3),	&
			 wzp(2)*p(3,3),	 wzp(3)*p(3,3),	&
			 (wxp(1)*p(3,1)+ wyp(1)*p(3,2)+wzp(1)*p(3,3)), 0.d0, 0.d0, -1.d0 /)

	C1=(  ck1*(a(2,2)*a(3,3)-a(2,3)*a(3,2))-ck2*(a(1,2)*a(3,3)-a(1,3)*a(3,2))+ck3*(a(1,2)*a(2,3)-a(2,2)*a(1,3)) )/dta
	C2=( -ck1*(a(2,1)*a(3,3)-a(2,3)*a(3,1))+ck2*(a(1,1)*a(3,3)-a(3,1)*a(1,3))-ck3*(a(1,1)*a(2,3)-a(2,1)*a(1,3)) )/dta
	C3=(  ck1*(a(2,1)*a(3,2)-a(2,2)*a(3,1))-ck2*(a(1,1)*a(3,2)-a(3,1)*a(1,2))+ck3*(a(1,1)*a(2,2)-a(1,2)*a(2,1)) )/dta
	do l=1,6
		call idx6(i,j,k,l,idx)
		if (it(l)==1) then
			indsrf=indsrf+1
			isrf(l,indirr(i,j,k))=indsrf
			isrfpts(indsrf,1:28:3)=iwts(1:10,1)
			isrfpts(indsrf,2:29:3)=iwts(1:10,2)
			isrfpts(indsrf,3:30:3)=iwts(1:10,3)
			
			if (l==1 .or. l==6) then
				srfpts(indsrf,1:10)=C3(1:10)
				srfftc(indsrf)=sum(c3(11:14)*rjp(1:4))
			elseif (l==2 .or. l==5) then
				srfpts(indsrf,1:10)=C2(1:10)
				srfftc(indsrf)=sum(c2(11:14)*rjp(1:4))
			else
				srfpts(indsrf,1:10)=C1(1:10)
				srfftc(indsrf)=sum(c1(11:14)*rjp(1:4))
			endif
			
		endif
	enddo
end
!-------------------------------------------------------------------------------
subroutine grd4p4n(isn,i,j,k,xs,ys,zs,rjp)
use comdata
use pbedata
implicit double precision(a-h,o-z) 
real*8 pcts(4), rjp(4), wts(10)
integer idx(3)
	isrf(0,indirr(i,j,k))=isn
	indsrf=indsrf+1; 
	if (io(i,j,k-1)==-isn)			then
		isrfpts(indsrf,1:10:3)=i;				isrfpts(indsrf,2:11:3)=j			
		isrfpts(indsrf,3:12:3)=(/k-2,k-1,k,k+1/);		pcts=z(k-2:k+1);				ik=1
	
	elseif (io(i,j-1,k)==-isn)		then
		isrfpts(indsrf,1:10:3)=i;				isrfpts(indsrf,2:11:3)=(/j-2,j-1,j,j+1/)	
		isrfpts(indsrf,3:12:3)=k;				pcts=y(j-2:j+1);						ik=2
	
	elseif (io(i-1,j,k)==-isn)		then
		isrfpts(indsrf,1:10:3)=(/i-2,i-1,i,i+1/);		isrfpts(indsrf,2:11:3)=j			
		isrfpts(indsrf,3:12:3)=k;				pcts=x(i-2:i+1);						ik=3
	
	elseif (io(i+1,j,k)==-isn)		then
		isrfpts(indsrf,1:10:3)=(/i+2,i+1,i,i-1/);	isrfpts(indsrf,2:11:3)=j
		isrfpts(indsrf,3:12:3)=k;				pcts=x(i+2:i-1:-1);						ik=4
	
	elseif (io(i,j+1,k)==-isn)		then
		isrfpts(indsrf,1:10:3)=i;				isrfpts(indsrf,2:11:3)=(/j+2,j+1,j,j-1/)
		isrfpts(indsrf,3:12:3)=k;				pcts=y(j+2:j-1:-1);						ik=5;
	else 
		isrfpts(indsrf,1:10:3)=i;				isrfpts(indsrf,2:11:3)=j
		isrfpts(indsrf,3:12:3)=(/k+2,k+1,k,k-1/);	pcts=z(k+2:k-1:-1);					ik=6;
	endif
	
	
	isrf(ik,indirr(i,j,k))=indsrf
	call idx6(i,j,k,ik,idx)
	
	if (io(i,j,k-1)==-isn .or. io(i,j,k+1)==-isn) then
		call findtggdwts(isn, idx, ik, i,j,k, zs, xs, ys, zs, pcts, wts, rjp(1), rjp(2))
	endif

	if (io(i,j-1,k)==-isn .or. io(i,j+1,k)==-isn) then
		call findtggdwts(isn, idx, ik, i,j,k, ys, xs, ys, zs, pcts, wts, rjp(1), rjp(2))
	endif	
	
	if (io(i-1,j,k)==-isn .or. io(i+1,j,k)==-isn) then
		call findtggdwts(isn, idx, ik, i,j,k, xs, xs, ys, zs, pcts, wts, rjp(1), rjp(2))
	endif
End
!-------------------------------------------------------------------------------
subroutine findgdwts(ipn,i,j,k,xs,ys,zs,wxp,wxn,wyp,wyn,wzp,wzn,iwts,it,id)
use pbedata
implicit double precision(a-h,o-z)

real*8 wxp(3),wxn(3),wyp(3),wyn(3),wzp(3),wzn(3),c(0:2,0:1), xx(0:2)
integer iwts(16,3), id(6), idx(3), idy(3), it(6)	

id=0; m=1; it=0

wxp=0.d0;	wxn=0.d0
wyp=0.d0;	wyn=0.d0
wzp=0.d0;	wzn=0.d0

do itest=1,6
	call idx6(i,j,k,itest,idx)
	call idx6(i,j,k,7-itest,idy)
	if (io(idx(1),idx(2),idx(3))==ipn .and. io(idy(1),idy(2),idy(3))==-ipn) it(itest)=1
	if (io(idx(1),idx(2),idx(3))==1 .and. io(idy(1),idy(2),idy(3))==1 .and. itest <= 3)	then
		it(itest)=1;	id(itest)=1
	endif
enddo	

	if (it(1)==1) then 
		call find1wts(z(k:k-2:-1),z(k+1:k-1:-1),zs,wzn,wzp)
		iwts(8:10,3)=(/k+1,k-1,k-2/)
	endif

	if (it(2)==1) then 
		call find1wts(y(j:j-2:-1),y(j+1:j-1:-1),ys,wyn,wyp)
		iwts(5:7,2)=(/j+1,j-1,j-2/)
	endif

	if (it(3)==1) then
		call find1wts(x(i:i-2:-1),x(i+1:i-1:-1),xs,wxn,wxp)
		iwts(1:4,1)=(/i+1,i,i-1,i-2/)
	endif

	if (it(4)==1) then
		call find1wts(x(i:i+2),x(i-1:i+1),xs,wxn,wxp)
		iwts(1:4,1)=(/i-1,i,i+1,i+2/)
	endif

	if (it(5)==1) then
		call find1wts(y(j:j+2),y(j-1:j+1),ys,wyn,wyp)
		iwts(5:7,2)=(/j-1,j+1,j+2/)
	endif

	if (it(6)==1) then
		call find1wts(z(k:k+2),z(k-1:k+1),zs,wzn,wzp)
		iwts(8:10,3)=(/k-1,k+1,k+2/)
	endif
end

!------------------------------------------------------------------------------------------
subroutine findtggdwts(isn, idx, ik, i, j, k, ps, xs, ys, zs, pcts, wts, phi, psi)
use pbedata
use comdata
implicit double precision(a-h,o-z)
real*8 pcts(1:4), btp(2), bts(2), wts(10), c(0:2,0:1), xx(0:2), p1n(3), p1p(3)
integer idx(3)
ijk=(k-1)*nx*ny+(j-1)*nx+i
m=1; c=0.d0
isfc=isrf(ik,indirr(i,j,k))

if (isn==1) then 
	xx(0:2)=pcts(1:3)
	call weights(ps,xx,2,2,m,c)
	p1n=c(:,1)
	if (pcts(1)>pcts(2)) p1n=-p1n

	xx(0:2)=pcts(4:2:-1)
	call weights(ps,xx,2,2,m,c)
	p1p=c(:,1)
	if (pcts(1)>pcts(2)) p1p=-p1p

	srfpts(isfc,1)=betan(xs,ys,zs)*p1n(1)
	srfpts(isfc,2)=betan(xs,ys,zs)*p1n(2)
	srfpts(isfc,3)=(betan(xs,ys,zs)*p1n(3) -betap(xs,ys,zs)*p1p(2))
	srfpts(isfc,4)=-betap(xs,ys,zs)*p1p(1)
	srfpts(isfc,1:4)=srfpts(isfc,1:4)/betap(xs,ys,zs)/p1p(3)
	srfftc(isfc)=(psi-p1p(2)*phi*betap(xs,ys,zs))/betap(xs,ys,zs)/p1p(3)

	bftc(ijk)	=bftc(ijk) + fj(xs,ys,zs) &
				+(betap(xs,ys+dy/2.d0,zs)+betap(xs,ys-dy/2.d0,zs)+betap(xs+dx/2.d0,ys,zs) &
				+betap(xs-dx/2.d0,ys,zs)+betap(xs,ys,zs-dz/2.d0)+betap(xs,ys,zs+dz/2.d0))/dcel**2*phi &
				-cappap(xs,ys,zs)*phi
else
	xx(0:2)=pcts(1:3)
	call weights(ps,xx,2,2,m,c)
	p1p=c(:,1)
	if (pcts(1)<pcts(2)) p1p=-p1p

	xx(0:2)=pcts(4:2:-1)
	call weights(ps,xx,2,2,m,c)
	p1n=c(:,1)
	if (pcts(1)<pcts(2)) p1n=-p1n

	srfpts(isfc,1)=betap(xs,ys,zs)*p1p(1)
	srfpts(isfc,2)=betap(xs,ys,zs)*p1n(2)
	srfpts(isfc,3)=(betap(xs,ys,zs)*p1p(3) -betan(xs,ys,zs)*p1n(2))
	srfpts(isfc,4)=-betan(xs,ys,zs)*p1n(1)
	srfpts(isfc,1:4)=srfpts(isfc,1:4)/betan(xs,ys,zs)/p1n(3)
	srfftc(isfc)=(-psi+p1p(3)*phi*betap(xs,ys,zs))/betan(xs,ys,zs)/p1n(3)

endif
end		

!-------------------------------------------------------------------------------------------
subroutine sasort(n,nra,ra,nrdc)
! Input:	n --- # of array items
!			nra --- postion of items
!           ra  --- weights of items
! Output:	nrb --- poition of items (ordered and reduced)
!			rb ---- weights of itmes (ordered and reduced)
implicit double precision(a-h,o-z)
integer nra(n),nrb(n)
real*8  ra(n),rb(n)

call sort3(n,nra,ra)
nrb=0; rb=0.d0
nrdc=0; itemp=-1           ! Only work for nra >0 None zero
do i=1,n
	if (nra(i) .ne. itemp) then
		nrdc=nrdc+1
		nrb(nrdc)=nra(i)
		rb(nrdc)=ra(i)	
		itemp=nrb(nrdc)
	else
		rb(nrdc)=rb(nrdc)+ra(i)
	endif
enddo
nra=nrb
ra=rb

End
!-------------------------------------------------------------------------------------------
subroutine fillmisspts
use comdata
use pbedata
use molecule
implicit double precision(a-h,o-z)	
integer idx1(3),idx2(3),itp(4)
real*8 cc(1:2),xx(3)

max_patch=int(nirr/2)
allocate(itps(max_patch,2),itpsnew(max_patch,8),idre(3,7))
allocate(icpatch(max_patch,7,50),akpatch(max_patch,7,50), pchftc(max_patch,7))
allocate(idpatch(max_patch),idpatchreal(max_patch,3))

itps=0;	lexch=0; idre=0
iph=0; idpatch=0; itnew=0;	
icpatch=0; akpatch=0.d0

!!-------------------------------------------------------------------
!!test
!i0=25;	j0=13;	k0=14
!print *,i0,j0,k0
!print *,io(i0-1:i0+1,j0,k0),irrpts(i0-1:i0+1,j0,k0)
!print *,io(i0,j0-1:j0+1,k0),irrpts(i0,j0-1:j0+1,k0)
!print *,io(i0,j0,k0-1:k0+1),irrpts(i0,j0,k0-1:k0+1)
!print *,io(i0-2:i0+2,j0,k0)
!print *,irrpts(i0-2:i0+2,j0,k0)
!!--------------------------------------------------------------------
!   goto 100
!	utest=0.d0
!	ixt=31 ; jyt=30 ; kzt=21 ; inet=2
!	print *,io(31:34,30,21)
!	print *,irrpts(31:34,30,21)
!	print *,iftpts(indirr(ixt,jyt,kzt),1:6,1)
!	print *,ixt,jyt,kzt,io(ixt,jyt,kzt)
!	do i=1,16
!	   i0=iftpts(indirr(ixt,jyt,kzt),inet,i*3-2)
!	   j0=iftpts(indirr(ixt,jyt,kzt),inet,i*3-1)
!	   k0=iftpts(indirr(ixt,jyt,kzt),inet,i*3)
!	   ijk=id3d(i0,j0,k0)
!	   utest=utest+uexa(ijk)*ftpts(indirr(ixt,jyt,kzt),inet,i)   
!	enddo
!	   utest=utest+ftc(indirr(ixt,jyt,kzt),inet)
!	   if (io(ixt,jyt,kzt)==0 .or. io(ixt,jyt,kzt)==-1) then
!		  uexact=up(x(ixt),y(jyt),z(kzt))
!	   else
!	   	  uexact=un(x(ixt),y(jyt),z(kzt))
!	   endif
!	   print *,'fic value of point',real(uexact),real(utest),abs(real(uexact-utest))
!	100 continue
!--------------------------------------------------------------------------------------

!itnew=0
itps=0;	itpsnew=0;	icpatch=0;	idpatch=0;	idpatchreal=0
akpatch=0.d0

nf=50
isum=0
ls=3
idre(1:3,1)=(/0,0,-1/)
idre(1:3,2)=(/0,-1,0/)
idre(1:3,3)=(/-1,0,0/)
idre(1:3,4)=(/1,0,0/)
idre(1:3,5)=(/0,1,0/)
idre(1:3,6)=(/0,0,1/)
idre(1:3,7)=(/0,0,0/)
!lexch=(/3,2,1,5,4,7,6/)

!print *,'Total Directional Replacement Applied In This Application:', sum(irrdrc)-nirrall*2
!print *,'total incomplete points', sum(iftwrg)
itrbl=0

do i=1,nx
	do j=1,ny
		do k=1,nz
			if (io(i,j,k)==-1 .and. irrpts(i,j,k)==1) then
				do ik=1,6
					if (iftwrg(indirr(i,j,k),ik)==1 .and. io(i,j,k)==-1)  then
						itrbl=itrbl+1
						inext=7-ik
						call settlecoefrp(i,j,k,inext)	
! 						print *, "test March 2024"
					endif
				enddo
			endif
		enddo
	enddo
enddo

call getfictitious

!Print *,'Total number of patches used by Sining', itrbl
!print *,'Total Directional Replacement=: ',isum
end

!---------------------------------------------------------------------------------------------
subroutine testftpts
use comdata
use pbedata
use bicg
implicit double precision(a-h,o-z)
integer ipcherr(3),ifcterr(3)
ipcherr = 0
!print *,'iph=', iph
fcterr=0.d0
do i=1,nx
	do j=1,ny
		do k=1,nz
			if (irrpts(i,j,k)==1) then
				do m=1,6
					if (io(i,j,k)==-1 .or. io(i,j,k)==0) then
						if (ibd==0) then
							ufc=up(x(i),y(j),z(k))
						elseif (ibd==1) then
							ufc=upi(i,j,k)	
						endif
					else
						if (ibd==0) then
							ufc=un(x(i),y(j),z(k))
! 							print *, "un at main ufc is "!, ufc
						elseif (ibd==1) then
							ufc=uni(i,j,k)
						endif
					endif
					uintp=0.d0
				
					if (abs(sum(ftpts(indirr(i,j,k),m,1:16)))>1.d-10 .and. abs(ftpts(indirr(i,j,k),m,nfc)) < 1.d-10) then
						do l=1,16
							
							if (abs(ftpts(indirr(i,j,k),m,l))>1.d-10) then
								ii=iftpts(indirr(i,j,k),m,3*l-2)
								jj=iftpts(indirr(i,j,k),m,3*l-1)
								kk=iftpts(indirr(i,j,k),m,3*l)
								uintp=uintp+ftpts(indirr(i,j,k),m,l)*uexa((kk-1)*nx*ny+(jj-1)*nx+ii)
							endif
							
						enddo
						uintp=uintp+ftc(indirr(i,j,k),m)
						!if (abs(uintp-ufc) > 7.d-2) then
						!	print *,i,j,k,m,real(uintp-ufc),real(varphi(x(i),y(j),z(k)))
						!endif
						if (fcterr < abs(uintp-ufc)) then 
							fcterr=abs(uintp-ufc)
							ifcterr=(/i,j,k/)
						endif
					endif


					if (iftpts(indirr(i,j,k),m,1)==-1) then	
						!print *,i,j,k,'patch'
						uintp=0.d0
						do iii=1,iph
							jjj=idpatch(iii)
							call line2cube(ipp,jpp,kpp,jjj)
							!print *,ipp,jpp,kpp,jjj
							if (i==ipp .and. j==jpp .and. k==kpp) then
								ipch=iii
								goto 20
							endif
						enddo
						print *,'cannot find pos of patched pt:', ipp,jpp,kpp
						20 continue
						!print *,'ipch= ',ipch
						do l=1,50
							if (abs(akpatch(ipch,m,l))>1.d-10 .and. (icpatch(ipch,m,l) .ne. 0)) then
								call line2cube(i0,j0,k0,icpatch(ipch,m,l))
								!if (i==21 .and. j==12 .and. k==11 .and. m==5) then
								!	print *,icpatch(ipch,m,l),i0,j0,k0,akpatch(ipch,m,l)
								!	print *,pchftc(ipch,m)
								!endif
								uintp=uintp+akpatch(ipch,m,l)*uexa((k0-1)*nx*ny+(j0-1)*nx+i0)
							endif
						enddo
						uintp=uintp+pchftc(ipch,m)
						!print *,'pch err',i,j,k,m,abs(uintp-ufc)	
						if (pcherr<abs(uintp-ufc)) then
							pcherr=abs(uintp-ufc)
							ipcherr=(/i,j,k/)
						endif
					endif
				enddo
			endif
		enddo
	enddo
enddo

write(*,*) 'The fct err is', real(fcterr), ' at ', ifcterr
write(*,*) 'The pch err is', real(pcherr), ' at ', ipcherr 
End

!---------------------------------------------------------------------------------------
function func(xyz,x,y,z,ixyz)
Implicit Double Precision(A-H, O-Z)
integer ixyz
Common	/pi/ pi
!ixyz=1,2,3 change in x,y,z direction respectively
if (ixyz==1) then
	x1=xyz;	y1=y;	z1=z
elseif (ixyz==2) then
	x1=x;	y1=xyz;	z1=z
else
	x1=x;	y1=y;	z1=xyz
endif
!func=x1**2+y1**2+z1**2-4.d0
func=varphi(x1,y1,z1)
end

FUNCTION rtbis(x1,x2,xacc,x,y,z,ixyz)
Implicit Double Precision(A-H, O-Z)
INTEGER JMAX,ixyz
REAL*8 rtbis,x1,x2,xacc,x,y,z
PARAMETER (JMAX=320) !Maximum allowed number of bisections.
					!Using bisection, find the root of a function func known to lie 
					!between x1 and x2. The root, returned as rtbis, will be refined 
					!until its accuracy is ¡Àxacc.
INTEGER j
REAL*8 dx,f,fmid,xmid
fmid=func(x2,x,y,z,ixyz)
f=func(x1,x,y,z,ixyz)
if(f*fmid.ge.0.) write(*,*) 'root must be bracketed in rtbisy'
if(f.lt.0.)then		!Orient the search so that f>0 lies at x+dx.
	rtbis=x1
	dx=x2-x1
else
	rtbis=x2
	dx=x1-x2
endif
do  j=1,JMAX		!Bisection loop.
	dx=dx*.5d0
	xmid=rtbis+dx
	fmid=func(xmid,x,y,z,ixyz)
	if	(fmid.le.0.) rtbis=xmid
	if(abs(dx).lt.xacc .or. fmid.eq.0.) return
enddo
print *,ixyz,x1,x2,rtbis
write(*,*) 'too many bisections in rtbis'

END

!----------------------------------------------------------------------------------------------------------
subroutine test
use pbedata
use molecule
use comdata
parameter(m=5)
implicit double precision(a-h,o-z)
real*8 xx(0:5),yy(0:5),c6(0:5,0:2),cc6(0:5),pxyz(3), p(3,3), dxyz(3), xd(-m:m), yd(-m:m), zd(-m:m), ww(-m:m,0:1)
real*8 xxx(0:m), yyy(0:m), zzz(0:m), rrr(0:m), cc(0:m,0:2), c(0:m), phi00(0:m), vvx(-m:m), vvy(-m:m), vvz(-m:m)
i0=13; j0=9; k0=9
xx(0)=x(i0)
dxx=dx/16.d0
do ii=1,5
	xx(ii)=xx(ii-1)-dxx
enddo

call weights(xx(0),xx,5,5,2,c6)
cc6=c6(0:5,2)

do ii=0,5
	pxyz=(/xx(ii),y(j0),z(k0)/)
	call kirk_potential(io(i0-ii,j0,k0),pxyz,ptl,pmt,psv)
	yy(ii)=pmt+psv
enddo

xs=sqrt(2.d0)/2.d0;	ys=sqrt(3.d0);	zs=sqrt(2.d0)/2.d0
xxx(0)=0.d0
yyy(0)=0.d0
zzz(0)=0.d0
rrr(0)=0.d0

mm=1
do i=1,m
	xxx(i)=xxx(0)+xs*dble(i)/dble(m*mm)
	yyy(i)=yyy(0)+ys*dble(i)/dble(m*mm)
	zzz(i)=zzz(0)+zs*dble(i)/dble(m*mm)
	rrr(i)=sqrt(xxx(i)**2+yyy(i)**2+zzz(i)**2)
enddo

call weights(rrr(m),rrr(0:m),m,m,1,cc)
rrr=cc(:,1)
do i=0,m
	phi00(i)=gx(xxx(i),yyy(i),zzz(i))	
enddo

print *,real(phi00)
flux=eps0*dot_product(rrr,phi00)
print *,'flux1= ',flux

call Transfmx(xs,ys,zs,P)
dxyz=(/4.d0*xs,12.d0*ys**3,4.d0/)
print *,dxyz
flux=dot_product(dxyz,P(1,1:3))
print *,'flux2=',flux	

do i=-m,m
	xd(i)=xs+dble(i)/dble(m)
	yd(i)=ys+dble(i)/dble(m)
	zd(i)=zs+dble(i)/dble(m)
	vvx(i)=gx(xd(i),ys,zs)
	vvy(i)=gx(xs,yd(i),zs)
	vvz(i)=gx(xs,ys,zd(i))
enddo

call weights(xd(0),xd,2*m,2*m,1,ww)
px=dot_product(ww(:,1),vvx)

call weights(yd(0),yd,2*m,2*m,1,ww)
py=dot_product(ww(:,1),vvy)

call weights(zd(0),zd,2*m,2*m,1,ww)
pz=dot_product(ww(:,1),vvz)

dxyz=(/px,py,pz/)
flux=dot_product(dxyz,P(1,1:3))
print *,dxyz
print *,'flux3=',flux

End 


function gx(xx,yy,zz)
implicit double precision(a-h,o-z)
gx=2.d0*xx**2+3.d0*yy**4+4.d0*zz
end
