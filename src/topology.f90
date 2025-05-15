SUBROUTINE READATOM_GRIDGEN
     USE MOLECULE
     USE COMDATA

! read in the MSMS files
      
     IMPLICIT REAL*8(A-H,O-Z) 
	 REAL*8  :: POS(3), VECTOR(3), RADIUS
     INTEGER :: NIND(5)
     CHARACTER(6)   :: CHAR1, CHAR2, CHAR3
     CHARACTER(100) :: FHEAD

! read the surface points

      OPEN(2,FILE=pathname(1:lenpath)//FNAME(1:lenfname)//".vert")
      READ(2,*) FHEAD
      READ(2,*) FHEAD
      READ(2,*) NSPT, ppp, qqq, rrr

      ALLOCATE(SPTPOS(3,NSPT), SPTNRM(3,NSPT), NATMAFF(NSPT), NSFTYPE(NSPT))
	  SPTPOS=0.D0; SPTNRM=0.D0; NATMAFF=0; NSFTYPE=0;	

      DO I=1,NSPT
         READ(2,*) POS(1:3), VECTOR(1:3), KK, NAFF, NAFFT 
         SPTPOS(:,I) = POS;   SPTNRM(:,I) = VECTOR
         NATMAFF(I)  = NAFF;  NSFTYPE(I)  = NAFFT
      END DO

      CLOSE(2)

! read the surface triangulization
      OPEN(3,FILE=pathname(1:lenpath)//FNAME(1:lenfname)//".face")
      READ(3,*) FHEAD
      READ(3,*) FHEAD
      READ(3,*) NFACE, PPP, QQQ, RRR
      ALLOCATE(NVERT(3,NFACE), MFACE(NFACE))
	  NVERT=0; MFACE=0
      DO I=1,NFACE 
         READ(3,*) NIND(1:5) 
         NVERT(1:3,I) = NIND(1:3);  MFACE(I) = NIND(4)
      END DO

      CLOSE(3)


! the mesh size DCEL should be specified elsewhere

      dcel1=0.5d0
	  nlay=max(7,int(1.d0/dcel1))

      xmin   = Minval(sptpos(1,1:nspt))
	  xmax   = Maxval(sptpos(1,1:nspt))
	  XLEFT  = int(xmin-nlay*dcel1)
      XRIGHT = int(xmax+nlay*dcel1)
      
      ymin   = Minval(sptpos(2,1:nspt))
	  ymax   = Maxval(sptpos(2,1:nspt))

	  YLEFT  = int(ymin-nlay*dcel1)
      YRIGHT = int(ymax+nlay*dcel1)

      zmin   = Minval(sptpos(3,1:nspt))
	  zmax   = Maxval(sptpos(3,1:nspt))
	  ZLEFT  = int(zmin-nlay*dcel1)
      ZRIGHT = int(zmax+nlay*dcel1)

      NX = INT( (XRIGHT - XLEFT)/DCEL) + 1
      NY = INT( (YRIGHT - YLEFT)/DCEL) + 1
      NZ = INT( (ZRIGHT - ZLEFT)/DCEL) + 1

	  nx=nx+1-mod(nx,2)
	  ny=ny+1-mod(ny,2)
	  nz=nz+1-mod(nz,2)

	  xright=xleft+(nx-1)*dcel
	  yright=yleft+(ny-1)*dcel
	  zright=zleft+(nz-1)*dcel

      MAX_NIRR = NX*NY*NZ/9
      NTOT = NX*NY*NZ
     
      ALLOCATE(XC(NX), YC(NY), ZC(NZ))
      
      DO I=1,NX
         XC(I) = (I-1)*DCEL + XLEFT         
      END DO
      
      DO J=1,NY
         YC(J) = (J-1)*DCEL + YLEFT
      END DO
      
      DO K=1,NZ
         ZC(K) = (K-1)*DCEL + ZLEFT
      END DO 
     
      DX = DCEL
      DY = DCEL
      DZ = DCEL
!  memory allocation
      print*,"allocating for irregular points"
      print*,"max_nirr=",max_nirr
      ALLOCATE(LPTIO(NX,NY,NZ), IRRPT(NX,NY,NZ), INDREG(NX,NY,NZ), INDIRR(NX,NY,NZ))
      ALLOCATE(IRRXYZ(3,MAX_NIRR))
	  ALLOCATE(MCX(MAX_NIRR),MCY(MAX_NIRR),MCZ(MAX_NIRR))    
    
	  mcx=0; mcy=0; mcz=0; lptio=0; irrpt=0; indreg=0; indirr=0; irrxyz=0
     
      ALLOCATE(DXL(MAX_NIRR), DXR(MAX_NIRR), DYL(MAX_NIRR), DYR(MAX_NIRR), DZL(MAX_NIRR), DZR(MAX_NIRR))

      DXL = 0.D0
      DXR = 0.D0
      DYL = 0.D0

      DYR = 0.D0
      DZL = 0.D0
      DZR = 0.D0

! the normal direction (cphi, sphi, ctheta, stheta) at all 3 possible interface points

! related with each irregular point

      ALLOCATE(CLOCAL(4,3,MAX_NIRR))     
      CLOCAL = 0.D0

END SUBROUTINE READATOM_GRIDGEN

!----------------------------------------------
! This subroutine set up charge location
! charge location = atom location if non-kirkwood 
subroutine setChgPos
use comdata
use molecule
use pbedata
implicit none
real*8 a2

if (ibd==1) then 
    !set number of charges for kirkwood case
    nchr=1
    ! Setup charge positions and amount
    allocate(chrpos(3,nchr),atmchr(nchr))
    atmchr=1.d0
    if (nchr==1) then 
        chrpos(:,1)=(/0.d0,0.d0,0.0d0/)
    elseif (nchr==2) then
        a2=1.5d0
        chrpos(:,1)=(/a2,0.d0,0.d0/)
        chrpos(:,2)=(/0.d0,0.d0,a2/)
    elseif (nchr==6) then
        chrpos(:,1)=(/0.2d0,0.2d0,0.2d0/)
        chrpos(:,2)=(/0.5d0,0.5d0,0.5d0/)
        chrpos(:,3)=(/0.8d0,0.8d0,0.8d0/)
        chrpos(:,4)=(/-0.2d0,0.2d0,-0.2d0/)
        chrpos(:,5)=(/0.5d0,-0.5d0,0.5d0/)
        chrpos(:,6)=(/-0.8d0,-0.8d0,-0.8d0/)
    endif
else
    nchr=natm    
    print *,nchr
    if (isf==0) then
        allocate(atmchr(nchr), chrpos(3,nchr))
        atmchr=1.d0
        if (nchr==2) atmchr(2)=0.d0 !test for lenhoff
        chrpos=atmpos
    elseif (isf==1) then
        !allocate(chrpos(3,nchr))
        !chrpos=atmpos
    endif
endif
end subroutine setChgPos
