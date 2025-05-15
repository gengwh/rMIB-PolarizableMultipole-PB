subroutine SETIRRPT2(imsmserr)
USE MOLECULE
USE COMDATA
IMPLICIT REAL*8(A-H,O-Z)


!********************************************************************************************************************

! Description of the program:

!

!    Implementation of MSMS surface in the MIB method. This program calculates

!    the accessibilty of all grid points(from which the dielectric function can

!    be calculated), marks all the irregular points, finds all the intersect points

!    of the mesh lines with the MSMS surface as well as the normal direction of the

!    MSMS surface at these intersect points.

!

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

!

!  Input data:

!    All the input data are passing through module MOLECULE

!

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

!

!  Output data:

!    All the output data are passing through module COMDATA

!                LPTIO:  the accessibility of grid points

!                        -1: accessible

!                         1: inaccessible

!                IRRPT:  marker of irregularity

!                         1: irregular point

!                         0: regular point

!                 NIRR:  number of irregular points

!                 NREG:  number of regular points

!               INDIRR:  the index of an irregular point in the array of all irregular points

!               INDREG:  the index of a regular point in the array of all regular points

!               IRRXYZ:  (IX,IY,IZ) of an irregular points

!          MCX,MCY,MCZ:  the state of intersection of surface with mesh line

!                        MCX = 1:  the surface acrosses the x-mesh line at x+

!                        MCX =-1:  the surface acrosses the x-mesh line at x-

!                        MCY = 1:  the surface acrosses the y-mesh line at y+

!                        MCY =-1:  the surface acrosses the y-mesh line at y-

!                        MCZ = 1:  the surface acrosses the z-mesh line at z+

!                        MCZ =-1:  the surface acrosses the z-mesh line at z-   

!                  DXL:  distance between irregular point (ix,iy,iz) and the associated x+ intersect point

!                  DXR:  DX - DXL

!                  DYL:  distance between irregular point (ix,iy,iz) and the associated y+ intersect point

!                  DYR:  DY - DYL

!                  DZL:  distance between irregular point (ix,iy,iz) and the associated z+ intersect point

!                  DZR:  DZ - DZL

!     CLOCAL(:,3,NIRR):  the normal directions at three mesh line/surface intersect points.

!                        

!                        *  an irregular points (x,y,z) may have at most 3 associated intersect points of mesh lines

!                        *  with surface; One at x+; the second at y+; the third at z+ 

!                        

!                        CLOCAL(:,1,:):    the normal direction at the intersect point of x-mesh line and surface

!                        CLOCAL(:,2,:):    the normal direction at the intersect point of y-mesh line and surface

!                        CLOCAL(:,3,:):    the normal direction at the intersect point of z-mesh line and surface

!                        CLOCAL(1:4,:,:):  (cos(theta),sin(theta),cos(phi),sin(phi))

!                                           theta: the angle between x-axis and the project of normal vector on x-y plane

!                                           phi:   the angle between z-axis and the normal vector

!

!          for exmaple:  check LPTIO(ix,iy,iz) to see it is outside the surface(-1) or inside the surface 

!                        check IRRPT(ix,iy,iz) to see whether (ix,iy,iz) is regular

!                        if it is irregular, the index of this irregular is INDIRR(ix,iy,iz), let it be INDEX

!                        check MCX(INDEX): if it is 1, then surface acrosses x-mesh line at x+

!                                          the distance between the x-intersect point and (ix,iy,iz) is DXL(INDEX)

!                                          the normal direction of surface at the intersect point is CLOCAL(:,1,INDEX)

!                        check MCY(INDEX): if it is 1, then surface acrosses y-mesh line at y+

!                                          the distance between the y-intersect point and (ix,iy,iz) is DYL(INDEX)

!                                          the normal direction of surface at the intersect point is CLOCAL(:,2,INDEX)

!                        check MCZ(INDEX): if it is 1, then surface acrosses z-mesh line at z+

!                                          the distance between the z-intersect point and (ix,iy,iz) is DZL(INDEX)

!                                          the normal direction of surface at the intersect point is CLOCAL(:,3,INDEX)

!

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

!

!  Local real variables:

!       AP, BP, CP, DP:  A, B, C, D in the equation Ax + By + Cx + D = 0 of the plane

!                        spanned by a surface triangle

!                 DPTB:  1.D0/SQRT(AP**2 + BP**2 + CP**2)

!            PLANENORM:  the unit normal vector of a surface triangle

!  VPOS1, VPOS2, VPOS3:  the coordinates of three vertices of a surface triangle

!                INDXF:  the index of a surface triangle in the .face file

!

!                XYZPT:  the coordinates of all intersect points of a surface triangle with all mesh lines

!                CIRCM:

!               POSMID:

!  SIDE1, SIDE2, SIDE3:  three sides of a surface triangle. Saved as vectors

!                 ATMC: 

!                 CXYZ:

!                 AXIS:  

!               SFNORM:  the normal direction of a surface triangle (cos(theta),sin(theta),cos(phi),sin(phi))

!              SFNORM2:  the normal direction of a surface point    (cos(theta),sin(theta),cos(phi),sin(phi))

!               XYSURF:  the coordiate of intersect points of each z-mesh line with MSMS surface

!               YZSURF:  the coordiate of intersect points of each x-mesh line with MSMS surface

!               XZSURF:  the coordiate of intersect points of each y-mesh line with MSMS surface

!              DXYSURF:  the normal direction of surface at each intersect point in XYSURF

!              DYZSURF:  the normal direction of surface at each intersect point in YZSURF

!              DXZSURF:  the normal direction of surface at each intersect point in XZSURF

!              PTSLINE:  temporary array for the coodinate of all intersect points on a mesh line

!             DPTSLINE:  temporary array for the normal direction at all intersect points on a mesh line

!               TVNORM:  temporary array for the normal vector of a surface triangle

!               AVNORM:  temporary vector. The average of three normal vectors at 

!                        three vertices of a surface triangle

!             SURFNORM:  the unit normal vector of a surface triangle 

!

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

!  Local integer variables:

!			   MXYSURF:  the state of intersection at each intersect points in XYSURF

!			   MYZSURF:  the state of intersection at each intersect points in YZSURF

!			   MXZSURF:  the state of intersection at each intersect points in XZSURF

!               MARKPT:  temporary array for the state of intersect points

!                   NV:  index of three vertices of a surface triangle in .vert file

!                MODPT:  the state of intersect points of a surface triangle with mesh lines

!               IXYZPT:  the index of mesh line intersecting with a surface triangle

!              MSPTSEC:  the state of intersection of a surface point

!            MSFACESPT:  the state of intersection of a surface triangle

!

!********************************************************************************************************************

              
     
     !COMMON/PLANE/AP,BP,CP,DP,DPTB,PLANENORM,VPOS1,VPOS2,VPOS3,INDXF
     REAL*8  :: DCEL9, PI
     REAL*8  :: ONES(3), VPOS1(3), VPOS2(3), VPOS3(3), PLANENORM(3), XYZPT(3,1000),& ! S.Y.
	            CIRCM(3,4), POSMID(3)
     REAL*8  :: SIDE1(3), SIDE2(3), SIDE3(3), ATMC(3), CXYZ(3), AXIS(3), SFNORM(4), SFNORM2(4)
     REAL*8  :: XYSURF(31,NX,NY), YZSURF(31,NY,NZ), XZSURF(31,NX,NZ), PTSLINE(31)
     REAL*8  :: DXYSURF(4,30,NX,NY), DYZSURF(4,30,NY,NZ), DXZSURF(4,30,NX,NZ)
     REAL*8  :: TVNORM(3), AVNORM(3), GPTV1(3), GPTV2(3)     
     REAL*8  :: SURFNORM(8,NFACE), DPTSLINE(4,30)
     
     INTEGER :: MXYSURF(30,NX,NY), MYZSURF(30,NY,NZ), MXZSURF(30,NX,NZ), MARKPT(30)
     INTEGER :: NV(3), MODPT(1000), IXYZPT(2,1000), MSPTSEC(3,NSPT), NFACESPT(31,NSPT) ! S.Y.
     
     imsmserr=0
     AP=0.d0;BP=0.d0;CP=0.d0;DP=0.d0;DPTB=0.d0;PLANENORM=0.d0;VPOS1=0.d0;VPOS2=0.d0;VPOS3=0.d0;INDXF=0.d0
     ONES = 1.D0
     PI   = ACOS(-1.D0)  
     DCEL9 = 9.D0*DCEL    
     
     SURFNORM = 0.D0
     MSPTSEC  = 0
     NFACESPT = 0
     
     XYSURF = 0.D0
     XZSURF = 0.D0
     YZSURF = 0.D0

     DXYSURF = 0.D0
     DXZSURF = 0.D0
     DYZSURF = 0.D0   

     MXYSURF = 0 
     MYZSURF = 0 
     MXZSURF = 0

	 !Weihua


    !print *,'begin to execute trig to cartisian...'
    !print *,natm,nface,dcel,xleft,yleft,zleft,nx,ny,nz,dx,dy,dz
            
!********************************************************************
! find the intersect points of every mesh line with the MSMS surface
!********************************************************************
       
     DO IFACE=1,NFACE
        MODPT  = 0
        XYZPT  = 0.D0
        IXYZPT = 0
!       three vertex of the triangle        
        NV(1:3) = NVERT(:,IFACE)
        VPOS1(1:3) = SPTPOS(1:3,NV(1))
        VPOS2(1:3) = SPTPOS(1:3,NV(2))
        VPOS3(1:3) = SPTPOS(1:3,NV(3))
        PDIFF12 = SUM(ABS(VPOS1 - VPOS2))
        PDIFF23 = SUM(ABS(VPOS1 - VPOS3))
        PDIFF13 = SUM(ABS(VPOS1 - VPOS3))
        IF(PDIFF12.LT.1.D-3.OR.PDIFF23.LT.1.D-3.OR.PDIFF13.LT.1.D-3) THEN
           GOTO 100
        END IF
!       find xmax, ymax, zmax, xmin, ymin, zmin       
        XMAX = MAX(VPOS1(1),VPOS2(1),VPOS3(1))
        XMIN = MIN(VPOS1(1),VPOS2(1),VPOS3(1))
                
        YMAX = MAX(VPOS1(2),VPOS2(2),VPOS3(2))
        YMIN = MIN(VPOS1(2),VPOS2(2),VPOS3(2))
        
        ZMAX = MAX(VPOS1(3),VPOS2(3),VPOS3(3))  
        ZMIN = MIN(VPOS1(3),VPOS2(3),VPOS3(3)) 
                 
!       find the possible intersect points of the plane with mesh lines

        IX1 = FLOOR((XMIN-XLEFT)/DX)+1; IX2 = CEILING((XMAX-XLEFT)/DX)+1
        IY1 = FLOOR((YMIN-YLEFT)/DY)+1; IY2 = CEILING((YMAX-YLEFT)/DY)+1
        IZ1 = FLOOR((ZMIN-ZLEFT)/DZ)+1; IZ2 = CEILING((ZMAX-ZLEFT)/DZ)+1

        IX1 = MAX(IX1,1); IX2 = MIN(NX,IX2)
        IY1 = MAX(IY1,1); IY2 = MIN(NY,IY2)
        IZ1 = MAX(IZ1,1); IZ2 = MIN(NZ,IZ2)
                
        IF(IX1.EQ.IX2.AND.IY1.EQ.IY2.AND.IZ1.EQ.IZ2) GOTO 100        
                      
!       set up the equation of the plane: ap*x + bp*y + cp*z +dp = 0 

        CALL EQPLANE(AP,BP,CP,DP,VPOS1,VPOS2,VPOS3)        
        DPTB = 1.D0/SQRT(AP**2 + BP**2 + CP**2)
        
		!if(idSurf==2)then
		!	planenorm(1:3)=plnnorm(iface,:)
		!else

			CALL NMPLANE(PLANENORM,VPOS1,VPOS2,VPOS3)
			AVNORM = (SPTNRM(:,NV(1)) + SPTNRM(:,NV(2)) + SPTNRM(:,NV(3)))/3.D0
			VNORM  = SQRT(SUM(AVNORM**2))
			AVNORM = AVNORM/VNORM
			PROJP  = SUM(PLANENORM*AVNORM)
			IF(PROJP.LT.0.D0) PLANENORM =-PLANENORM
		!endif
       
        SURFNORM(1,  IFACE) = AP
        SURFNORM(2,  IFACE) = BP
        SURFNORM(3,  IFACE) = CP
        SURFNORM(4,  IFACE) = DP
        SURFNORM(5,  IFACE) = DPTB
        SURFNORM(6:8,IFACE) = PLANENORM
                                     
		PROJXY = SQRT( PLANENORM(1)**2 + PLANENORM(2)**2 )
        INDXF  = IFACE
        
        SFNORM(1)   = PROJXY
        SFNORM(2)   = PLANENORM(3)
        
        IF(PROJXY.LE.1.d-6) THEN
           SFNORM(3) = 1.D0
           SFNORM(4) = 0.D0
        ELSE
           SFNORM(3) = PLANENORM(1)/PROJXY
           SFNORM(4) = PLANENORM(2)/PROJXY      
        END IF

        NPINT = 0
        XYZPT = 0.D0
        IF(ABS(CP).GT.1.D-14) THEN 
        DO IPX=IX1,IX2
           XX = XC(IPX)
           DO IPY=IY1,IY2
              YY = YC(IPY)
              ZZ =-(DP+AP*XX+BP*YY)/CP
              NPINT = NPINT + 1

              XYZPT(1,NPINT)  = XX
              XYZPT(2,NPINT)  = YY
              XYZPT(3,NPINT)  = ZZ
              IXYZPT(1,NPINT) = IPX
              IXYZPT(2,NPINT) = IPY
              
              MODPT(NPINT)    = 1
           END DO
        END DO
        END IF
        
        IF(ABS(BP).GT.1.D-14) THEN
        DO IPX=IX1,IX2
           XX = XC(IPX)
           DO IPZ=IZ1,IZ2
              ZZ = ZC(IPZ)
              YY =-(DP+AP*XX+CP*ZZ)/BP
              NPINT = NPINT + 1

              XYZPT(1,NPINT)  = XX
              XYZPT(2,NPINT)  = YY
              XYZPT(3,NPINT)  = ZZ
              IXYZPT(1,NPINT) = IPZ
              IXYZPT(2,NPINT) = IPX
              MODPT(NPINT)    = 2
           END DO
        END DO
        END IF
        IF(ABS(AP).GT.1.D-14) THEN
        DO IPY=IY1,IY2
           YY = YC(IPY)
           DO IPZ=IZ1,IZ2
              ZZ = ZC(IPZ)
              XX =-(DP+BP*YY+CP*ZZ)/AP
              NPINT = NPINT + 1

              XYZPT(1,NPINT)  = XX
              XYZPT(2,NPINT)  = YY
              XYZPT(3,NPINT)  = ZZ
              IXYZPT(1,NPINT) = IPY
              IXYZPT(2,NPINT) = IPZ
              MODPT(NPINT)    = 3
           END DO
        END DO
        END IF

! NPINT: number of intersect points with the plane spanned by current triangle
! XYZPT: coordinates of the intersect points with the the plane spanned by current triangle

! to rule out these intersect points not located in the triangle
        DO IPT=1,NPINT
           
           SIDE1 = VPOS1 - XYZPT(:,IPT)
           SIDE2 = VPOS2 - XYZPT(:,IPT)
           SIDE3 = VPOS3 - XYZPT(:,IPT)   

           SIDE1L = DSQRT(SUM(SIDE1**2))
           SIDE2L = DSQRT(SUM(SIDE2**2))
           SIDE3L = DSQRT(SUM(SIDE3**2))  
               
!           IF(SIDE1L+SIDE2L+SIDE3L.GT.DCEL9)  THEN   
!              MODPT(IPT) = 0
!              GOTO 101
!           END IF

! a vertex is a intersect point

           IF(SIDE1L.LE.1.d-6) THEN
              MODPT(IPT) =-NV(1) - MODPT(IPT)*1000000
              GOTO 101
           ELSEIF(SIDE2L.LE.1.d-6) THEN
              MODPT(IPT) =-NV(2) - MODPT(IPT)*1000000
              GOTO 101
           ELSEIF(SIDE3L.LE.1.d-6) THEN
              MODPT(IPT) =-NV(3) - MODPT(IPT)*1000000
              GOTO 101                      
           END IF

           THETA1 = SUM(SIDE1*SIDE2)/(SIDE1L*SIDE2L)
           THETA2 = SUM(SIDE2*SIDE3)/(SIDE2L*SIDE3L)
           THETA3 = SUM(SIDE3*SIDE1)/(SIDE3L*SIDE1L)

           IF(ABS(THETA1).GT.1.D0) THETA1 = THETA1/ABS(THETA1)
           IF(ABS(THETA2).GT.1.D0) THETA2 = THETA2/ABS(THETA2)
           IF(ABS(THETA3).GT.1.D0) THETA3 = THETA3/ABS(THETA3)
   
           THETA1 = DACOS(THETA1)
           THETA2 = DACOS(THETA2)
           THETA3 = DACOS(THETA3)
           
           IF(ABS(THETA1 + THETA2 + THETA3 - 2.D0*PI).GT.1.d-6) MODPT(IPT) = 0

           IF(ABS(THETA1-PI).LE.1.d-6) THEN
              MODPT(IPT) = MODPT(IPT) + 10
           ELSEIF(ABS(THETA2-PI).LE.1.d-6) THEN
              MODPT(IPT) = MODPT(IPT) + 20
           ELSEIF(ABS(THETA3-PI).LE.1.d-6) THEN
              MODPT(IPT) = MODPT(IPT) + 30
           END IF

101        CONTINUE

           MODXYZ = 0
           MODPTC = MODPT(IPT)          

           IF(MODPTC.GT.0) THEN
              IF(MODPTC.GT.10) THEN   ! mesh line intersects with side of the triangle
                 IF(MODPTC.GT.30) THEN
                    MODXYZ  = MODPTC - 30
                 ELSEIF(MODPTC.GT.20) THEN
                    MODXYZ  = MODPTC - 20
                 ELSE
                    MODXYZ  = MODPTC - 10
                 END IF
              ELSE
                 MODXYZ = MODPTC
              END IF
           ELSEIF(MODPTC.LT.0) THEN   ! mesh line intersects with vertex of the triangle
              MODXYZ = FLOOR(-1.D0*MODPTC/1000000)
              NPOSPT = MOD(-MODPTC,1000000)
              TVNORM = SPTNRM(:,NPOSPT)
              PROJXY = SQRT( TVNORM(1)**2 + TVNORM(2)**2 )
        
              SFNORM2(1)   = PROJXY
              SFNORM2(2)   = TVNORM(3)
        
              IF(PROJXY.LE.1.D-14) THEN
                 SFNORM2(3) = 1.D0
                 SFNORM2(4) = 0.D0
              ELSE
                 SFNORM2(3) = TVNORM(1)/PROJXY
                 SFNORM2(4) = TVNORM(2)/PROJXY
              END IF
              
              NN = NFACESPT(1,NPOSPT)
              MREAP = 0
              I     = 1
              DO WHILE(MREAP.EQ.0.AND.I.LE.NN)
                 IF(IFACE.EQ.NFACESPT(I+1,NPOSPT)) THEN
                    MREAP = 1                    
                 END IF
                 I = I + 1
              END DO
              
              IF(MREAP.EQ.0) THEN
                 NN = NN + 1
                 NFACESPT(1,NPOSPT)    = NN
                 NFACESPT(NN+1,NPOSPT) = IFACE
              END IF
  
              IF(MSPTSEC(MODXYZ,NPOSPT).EQ.0) THEN
                 MSPTSEC(MODXYZ,NPOSPT) = 1
              ELSE
                 MODXYZ = 0
              END IF
           END IF
          
           IF(MODXYZ.NE.0) NSECT = NSECT + 1
                      
           SELECT CASE(MODXYZ)             
           CASE(1)   ! z-mesh line
              IRRX = IXYZPT(1,IPT)
              IRRY = IXYZPT(2,IPT)
              NSET = INT(XYSURF(1,IRRX,IRRY)) + 1
              IF(NSET.GT.30) THEN
		 print *,' More than 30 mesh line/triangle intesect points'
		 imsmserr=2
		 return
	      ENDIF
              XYSURF(1,IRRX,IRRY) = NSET
              XYSURF(NSET+1,IRRX,IRRY) = XYZPT(3,IPT)
              IF(MODPTC.GT.0) THEN
                 DXYSURF(:,NSET,IRRX,IRRY) = SFNORM
              ELSEIF(MODPTC.LT.0) THEN
                 DXYSURF(:,NSET,IRRX,IRRY) = SFNORM2
              END IF
              
              IF(MODPTC.GT.10) THEN
                 MXYSURF(NSET,IRRX,IRRY) = IFACE
              ELSEIF(MODPTC.LT.0) THEN
                 MXYSURF(NSET,IRRX,IRRY) =-NPOSPT
              END IF              
           CASE(2)   ! y-mesh line
              IRRZ = IXYZPT(1,IPT)
              IRRX = IXYZPT(2,IPT) 
              NSET = INT(XZSURF(1,IRRX,IRRZ)) + 1
              IF(NSET.GT.30) THEN
		 print *,' More than 30 mesh line/triangle intesect points'
		 imsmserr=2
		 return
              ENDIF
	      XZSURF(1,IRRX,IRRZ) = NSET
              XZSURF(NSET+1,IRRX,IRRZ) = XYZPT(2,IPT)
              IF(MODPTC.GT.0) THEN
                 DXZSURF(:,NSET,IRRX,IRRZ) = SFNORM
              ELSEIF(MODPTC.LT.0) THEN
                 DXZSURF(:,NSET,IRRX,IRRZ) = SFNORM2
              END IF
              
              IF(MODPTC.GT.10) THEN
                 MXZSURF(NSET,IRRX,IRRZ) = IFACE
              ELSEIF(MODPTC.LT.0) THEN
                 MXZSURF(NSET,IRRX,IRRZ) =-NPOSPT
             END IF
           CASE(3)    ! x-mesh line
              IRRY = IXYZPT(1,IPT)
              IRRZ = IXYZPT(2,IPT)
              NSET = INT(YZSURF(1,IRRY,IRRZ)) + 1
              IF(NSET.GT.30) THEN
		 print *,' More than 30 mesh line/triangle intesect points'
		 imsmserr=2
		 return
              ENDIF
              YZSURF(1,IRRY,IRRZ) = NSET
              YZSURF(NSET+1,IRRY,IRRZ) = XYZPT(1,IPT)
              IF(MODPTC.GT.0) THEN
                 DYZSURF(:,NSET,IRRY,IRRZ) = SFNORM
              ELSEIF(MODPTC.LT.0) THEN
                 DYZSURF(:,NSET,IRRY,IRRZ) = SFNORM2
              END IF
              
              IF(MODPTC.GT.10) THEN
                 MYZSURF(NSET,IRRY,IRRZ) = IFACE
              ELSEIF(MODPTC.LT.0) THEN
                 MYZSURF(NSET,IRRY,IRRZ) =-NPOSPT
              END IF              
           END SELECT
        END DO        
    
100     CONTINUE

     END DO  

!print *,'determine the nature of intersection of mesh line with grid points or sides'
!****************************************************************************
! determine the nature of intersection of mesh line with grid points or sides
!****************************************************************************
     
! Z-mesh line
     DO I=1,NX
        DO J=1,NY
 		   NSET = INT(XYSURF(1,I,J))         
 		   MARKPT = 0
           DO K=1,NSET
              IF(MXYSURF(K,I,J).GT.0.AND.MARKPT(K).EQ.0) THEN
                 IFACE1 = MXYSURF(K,I,J)
                 DO KK=1,NSET
                    IF(KK.NE.K.AND.ABS(XYSURF(KK+1,I,J)-XYSURF(K+1,I,J)).LE.1.D-6) THEN
                       MARKPT(KK) = 1
                       IFACE2 = MXYSURF(KK,I,J)
                       CALL PASSORTANGSID(THETA1,THETA2,IFACE1,IFACE2,SURFNORM,1)
                       IF(THETA1*THETA2.GT.0.D0) THEN
                          XYSURF(KK+1,I,J) = 1.D10
                          MARKPT(KK)       = 1
                       ELSEIF(THETA1*THETA2.EQ.0.D0) THEN
                          WRITE(*,'(2f12.5,a)') theta1, theta2,'  Ohhhhhhh z'
                       END IF                       
                    END IF
                 END DO
              ELSEIF(MXYSURF(K,I,J).LT.0) THEN
                 NPOSPT = MXYSURF(K,I,J)
                 CALL PASSORTANGSPT(MULTP,NPOSPT,NFACESPT,SURFNORM,1, AP,BP,CP,DP,DPTB,PLANENORM,VPOS1,VPOS2,VPOS3,INDXF)
				 IF(MULTP.EQ.1) THEN
					NN = INT(XYSURF(1,I,J))
                    XYSURF(1,I,J)       = NN + 1
                    XYSURF(NN+2,I,J)    = XYSURF(K+1,I,J)
                    DXYSURF(:,NN+1,I,J) = DXYSURF(:,K,I,J)
                 END IF
              END IF
           END DO
        END DO
     END DO
! Y-mesh line
     DO I=1,NX
        DO K=1,NZ
           NSET = INT(XZSURF(1,I,K))           
           MARKPT = 0
           
           DO J=1,NSET              
              IF(MXZSURF(J,I,K).GT.0.AND.MARKPT(J).EQ.0) THEN    ! pass a side
                 IFACE1 = MXZSURF(J,I,K)
                 DO JJ=1,NSET
                    IF(JJ.NE.J.AND.ABS(XZSURF(JJ+1,I,K)-XZSURF(J+1,I,K)).LE.1.D-6) THEN
                       MARKPT(JJ) = 1
                       IFACE2 = MXZSURF(JJ,I,K)
                       CALL PASSORTANGSID(THETA1,THETA2,IFACE1,IFACE2,SURFNORM,2)
                       IF(THETA1*THETA2.GT.0.D0) THEN
                          XZSURF(JJ+1,I,K) = 1.D10
                          MARKPT(JJ)       = 1
                       ELSEIF(THETA1*THETA2.EQ.0.D0) THEN
                          WRITE(*,'(2f12.5,a)') theta1, theta2,'  Ohhhhhhh y'
                       END IF                       
                    END IF
                 END DO
              ELSEIF(MXZSURF(J,I,K).LT.0) THEN   ! pass a vertex
                 NPOSPT = MXZSURF(J,I,K)
                 CALL PASSORTANGSPT(MULTP,NPOSPT,NFACESPT,SURFNORM,2, AP,BP,CP,DP,DPTB,PLANENORM,VPOS1,VPOS2,VPOS3,INDXF)
                 
                 IF(MULTP.EQ.1) THEN
                    NN = INT(XZSURF(1,I,K))
                    XZSURF(1,I,K)       = NN + 1
                    XZSURF(NN+2,I,K)    = XZSURF(J+1,I,K)
                    DXZSURF(:,NN+1,I,K) = DXZSURF(:,J,I,K)
                 END IF
              END IF
           END DO
           
        END DO
     END DO   
! X-mesh line

     DO J=1,NY
        DO K=1,NZ
           NSET = INT(YZSURF(1,J,K))
           MARKPT = 0

           DO I=1,NSET
              IF(MYZSURF(I,J,K).GT.0.AND.MARKPT(I).EQ.0) THEN
                 IFACE1 = MYZSURF(I,J,K)
                 DO II=1,NSET
                    IF(II.NE.I.AND.ABS(YZSURF(II+1,J,K)-YZSURF(I+1,J,K)).LE.1.D-6) THEN
                       MARKPT(II) = 1
                       IFACE2 = MYZSURF(II,J,K)
                       CALL PASSORTANGSID(THETA1,THETA2,IFACE1,IFACE2,SURFNORM,3)
                       IF(THETA1*THETA2.GT.0.D0) THEN
                          YZSURF(II+1,J,K) = 1.D10
                          MARKPT(II)       = 1
                       ELSEIF(THETA1*THETA2.EQ.0.D0) THEN
                          WRITE(*,'(2f12.5,a)') theta1, theta2,'  Ohhhhhhh x'
                       END IF                       
                    END IF
                 END DO
              ELSEIF(MYZSURF(I,J,K).LT.0) THEN
                 NPOSPT = MYZSURF(I,J,K)
                 CALL PASSORTANGSPT(MULTP,NPOSPT,NFACESPT,SURFNORM,3, AP,BP,CP,DP,DPTB,PLANENORM,VPOS1,VPOS2,VPOS3,INDXF)
                 IF(MULTP.EQ.1) THEN
                    NN = INT(YZSURF(1,J,K))
                    YZSURF(1,J,K)       = NN + 1
                    YZSURF(NN+2,J,K)    = YZSURF(I+1,J,K)  
                    DYZSURF(:,NN+1,J,K) = DYZSURF(:,I,J,K)                
                 END IF
              END IF
           END DO 
        END DO
     END DO                
!print *, 'determine the accessibility of grid points and locate the irregular points'
!****************************************************************************
! determine the accessibility of grid points and locate the irregular points
!****************************************************************************

 ! Z-mesh line
     DO I=1,NX
        DO J=1,NY
           PTSLINE  = XYSURF(:,I,J)
           DPTSLINE = DXYSURF(:,:,I,J)
		   CALL SORTPTS(PTSLINE,DPTSLINE,ZRIGHT)
           NSET = INT(INT(PTSLINE(1)))
           IF(MOD(NSET,2).NE.0) THEN
              WRITE(*,'(3I6,A)') I,J,NSET,' Z WRONG'
			  print *,ptsline(1)
              imsmserr=1
	      return 
           END IF

           IZ1  = 1
           IADD = 1 
           IF(NSET.GT.0) THEN 
              DO IS=1,NSET
                 IADD =-1*IADD
                 IZ2  = FLOOR((PTSLINE(IS+1)-ZLEFT)/DZ) + 1
                 IF(ABS(PTSLINE(IS+1)-ZC(IZ2+1)).LT.1.D-6) IZ2 = IZ2 + 1
                 IF(IZ2.GE.IZ1) LPTIO(I,J,IZ1:IZ2) = IADD
                 IF(ABS(ZC(IZ2)-PTSLINE(IS+1)).LE.1.d-6) LPTIO(I,J,IZ2) = 0
                 IZ1 = IZ2 + 1
              END DO 
              LPTIO(I,J,IZ1:NZ) =-1
              IS1 = 1  

              DO K=1,NZ-1
		 IF(LPTIO(I,J,K)*LPTIO(I,J,K+1).EQ.-1) THEN
                    IS = IS1
                    INDXSET = 0	
                    DO WHILE(IS.LE.NSET.AND.INDXSET.EQ.0)
                       IF(PTSLINE(IS+1).GT.ZC(K).AND.PTSLINE(IS+1).LT.ZC(K+1)) THEN
                          INDXSET = IS
                       END IF
                       IS = IS + 1
                    END DO
                    IS1 = IS1 + 1
                    
                    IF(INDXSET.EQ.0) WRITE(*,*) '  Can not find the intersest point for z-mesh line'

                    IF(IRRPT(I,J,K).EQ.0) THEN
                       IRRPT(I,J,K)   = 1
                       NIRR = NIRR + 1
                       INDIRR(I,J,K)  = NIRR
                       IRRXYZ(1,NIRR) = I
                       IRRXYZ(2,NIRR) = J
                       IRRXYZ(3,NIRR) = K
                    END IF
                    IIRR = INDIRR(I,J,K)
                    MCZ(IIRR) = 1
                    DZL(IIRR) = PTSLINE(INDXSET+1) - ZC(K)
                    DZR(IIRR) = ZC(K+1) - PTSLINE(INDXSET+1)
                    CLOCAL(:,3,IIRR) = DPTSLINE(:,INDXSET)

                    IF(IRRPT(I,J,K+1).EQ.0) THEN
                       IRRPT(I,J,K+1) = 1
                       NIRR = NIRR + 1
                       INDIRR(I,J,K+1) = NIRR
                       IRRXYZ(1,NIRR)  = I
                       IRRXYZ(2,NIRR)  = J
                       IRRXYZ(3,NIRR)  = K+1
                    END IF
                    IIRR = INDIRR(I,J,K+1)
                    MCZ(IIRR) =-1
                 ELSEIF(LPTIO(I,J,K).EQ.0) THEN
                    IS = IS1
                    INDXSET = 0	
                    DO WHILE(IS.LE.NSET.AND.INDXSET.EQ.0)
                       IF(ABS(PTSLINE(IS+1)-ZC(K)).LE.1.d-6) THEN
                          INDXSET = IS
                       END IF
                       IS = IS + 1
                    END DO
                    IS1 = IS1 + 1
                    
                    IF(INDXSET.EQ.0) WRITE(*,*) '  Can not find the intersest point for z-mesh line'

                    IF(IRRPT(I,J,K).EQ.0) THEN
                       NIRR = NIRR + 1
                       INDIRR(I,J,K)  = NIRR
                       IRRXYZ(1,NIRR) = I
                       IRRXYZ(2,NIRR) = J
                       IRRXYZ(3,NIRR) = K
                    END IF
                    IIRR = INDIRR(I,J,K)
                    IRRPT(I,J,K) =-1
                    CLOCAL(:,3,IIRR) = DPTSLINE(:,INDXSET)
                 END IF                  
              END DO
           ELSE
               LPTIO(I,J,:) =-1
           END IF   
          
        END DO
     END DO     

 ! Y-mesh line

     DO I=1,NX
        DO K=1,NZ
           PTSLINE = XZSURF(:,I,K)
           DPTSLINE = DXZSURF(:,:,I,K)
           CALL SORTPTS(PTSLINE,DPTSLINE,YRIGHT)
           NSET = INT(PTSLINE(1))
           IF(MOD(NSET,2).NE.0) THEN
              WRITE(*,'(3I6,A)') I,K,NSET,' Y WRONG'
              imsmserr=1
              return
	   END IF

           IY1  = 1
           IADD = 1 
           IF(NSET.GT.0) THEN 
              DO IS=1,NSET
                 IADD =-1*IADD
                 IY2  = FLOOR((PTSLINE(IS+1)-YLEFT)/DY) + 1
                 IF(ABS(PTSLINE(IS+1)-YC(IY2+1)).LT.1.d-6) IY2 = IY2 + 1
                 IF(IY2.GE.IY1) LPTIO(I,IY1:IY2,K) = IADD
                 IF(ABS(YC(IY2)-PTSLINE(IS+1)).LE.1.d-6) LPTIO(I,IY2,K) = 0
                 IY1 = IY2 + 1
              END DO 
              LPTIO(I,IY1:NY,K) =-1
              IS1 = 1

              DO J=1,NY-1
                 IF(LPTIO(I,J,K)*LPTIO(I,J+1,K).EQ.-1) THEN
                    IS = IS1
                    INDXSET = 0	 
                    DO WHILE(IS.LE.NSET.AND.INDXSET.EQ.0)
                       IF(PTSLINE(IS+1).GT.YC(J).AND.PTSLINE(IS+1).LT.YC(J+1)) THEN
                          INDXSET = IS
                       END IF
                       IS = IS + 1
                    END DO
                    IS1 = IS1 + 1
                    
                    IF(INDXSET.EQ.0) WRITE(*,*) '  Can not find the intersest point for y-mesh line'

                    IF(IRRPT(I,J,K).EQ.0) THEN
                       IRRPT(I,J,K)   = 1
                       NIRR = NIRR + 1
                       INDIRR(I,J,K)  = NIRR
                       IRRXYZ(1,NIRR) = I
                       IRRXYZ(2,NIRR) = J
                       IRRXYZ(3,NIRR) = K
                    END IF
                    IIRR = INDIRR(I,J,K)
                    MCY(IIRR) = 1
                    DYL(IIRR) = PTSLINE(INDXSET+1) - YC(J)
                    DYR(IIRR) = YC(J+1) - PTSLINE(INDXSET+1)
                    CLOCAL(:,2,IIRR) = DPTSLINE(:,INDXSET)
                    IF(IRRPT(I,J+1,K).EQ.0) THEN
                       IRRPT(I,J+1,K) = 1
                       NIRR = NIRR + 1
                       INDIRR(I,J+1,K) = NIRR
                       IRRXYZ(1,NIRR)  = I
                       IRRXYZ(2,NIRR)  = J+1
                       IRRXYZ(3,NIRR)  = K
                    END IF
                    IIRR = INDIRR(I,J+1,K)
                    MCY(IIRR) =-1                     
                 ELSEIF(LPTIO(I,J,K).EQ.0) THEN
                    IS = IS1
                    INDXSET = 0	 
                    DO WHILE(IS.LE.NSET.AND.INDXSET.EQ.0)
                       IF(ABS(PTSLINE(IS+1)-YC(J)).LE.1.d-6) THEN
                          INDXSET = IS
                       END IF
                       IS = IS + 1
                    END DO
                    IS1 = IS1 + 1
                    
                    IF(INDXSET.EQ.0) WRITE(*,*) '  Can not find the intersest point for y-mesh line'

                    IF(IRRPT(I,J,K).EQ.0) THEN
                       NIRR = NIRR + 1
                       INDIRR(I,J,K)  = NIRR
                       IRRXYZ(1,NIRR) = I
                       IRRXYZ(2,NIRR) = J
                       IRRXYZ(3,NIRR) = K
                    END IF
                    IIRR = INDIRR(I,J,K)
                    IRRPT(I,J,K) =-1
                    CLOCAL(:,2,IIRR) = DPTSLINE(:,INDXSET)                    
                 END IF                  
              END DO
           ELSE
              LPTIO(I,:,K) =-1
           END IF   

        END DO
     END DO

 ! X-mesh line

     DO J=1,NY
        DO K=1,NZ
        
           PTSLINE  = YZSURF(:,J,K)
           DPTSLINE = DYZSURF(:,:,J,K)

           CALL SORTPTS(PTSLINE,DPTSLINE,XRIGHT)
           
           NSET = INT(PTSLINE(1))

           IF(MOD(NSET,2).NE.0) THEN
              WRITE(*,'(3I6,A)') J,K,NSET,' X WRONG'
              imsmserr=1
              return
           END IF

           IX1  = 1
           IADD = 1 
           IF(NSET.GT.0) THEN 
              DO IS=1,NSET
                 IADD =-1*IADD
                 IX2  = FLOOR((PTSLINE(IS+1)-XLEFT)/DX) + 1
                 IF(ABS(PTSLINE(IS+1)-XC(IX2+1)).LT.1.d-6) IX2 = IX2 + 1
                 IF(IX2.GE.IX1) LPTIO(IX1:IX2,J,K) = IADD
                 IF(ABS(XC(IX2)-PTSLINE(IS+1)).LE.1.d-6) LPTIO(IX2,J,K) = 0
                 IX1 = IX2 + 1
              END DO 
              LPTIO(IX1:NX,J,K) =-1
              IS1 = 1 

              DO I=1,NX-1
                 IF(LPTIO(I,J,K)*LPTIO(I+1,J,K).EQ.-1) THEN
                    IS = IS1
                    INDXSET = 0	
                    DO WHILE(IS.LE.NSET.AND.INDXSET.EQ.0)
                       IF(PTSLINE(IS+1).GT.XC(I).AND.PTSLINE(IS+1).LT.XC(I+1)) THEN
                          INDXSET = IS
                       END IF
                       IS = IS + 1
                    END DO
                    IS1 = IS1 + 1
                    
                    IF(INDXSET.EQ.0) WRITE(*,*) '  Can not find the intersest point for x-mesh line'

                    IF(IRRPT(I,J,K).EQ.0) THEN
                       IRRPT(I,J,K)   = 1
                       NIRR = NIRR + 1
                       INDIRR(I,J,K)  = NIRR
                       IRRXYZ(1,NIRR) = I
                       IRRXYZ(2,NIRR) = J
                       IRRXYZ(3,NIRR) = K
                    END IF
                    IIRR = INDIRR(I,J,K)
                    MCX(IIRR) = 1
                    DXL(IIRR) = PTSLINE(INDXSET+1) - XC(I)
                    DXR(IIRR) = XC(I+1) - PTSLINE(INDXSET+1)
                    CLOCAL(:,1,IIRR) = DPTSLINE(:,INDXSET)

                    IF(IRRPT(I+1,J,K).EQ.0) THEN
                       IRRPT(I+1,J,K) = 1
                       NIRR = NIRR + 1
                       INDIRR(I+1,J,K) = NIRR
                       IRRXYZ(1,NIRR)  = I+1
                       IRRXYZ(2,NIRR)  = J
                       IRRXYZ(3,NIRR)  = K
                    END IF
                    IIRR = INDIRR(I+1,J,K)
                    MCX(IIRR) =-1
                 ELSEIF(LPTIO(I,J,K).EQ.0) THEN
                    IS = IS1
                    INDXSET = 0	 
                    DO WHILE(IS.LE.NSET.AND.INDXSET.EQ.0)
                       IF(ABS(PTSLINE(IS+1)-XC(I)).LE.1.d-6) THEN
                          INDXSET = IS
                       END IF
                       IS = IS + 1
                    END DO
                    IS1 = IS1 + 1
                    
                    IF(INDXSET.EQ.0) WRITE(*,*) '  Can not find the intersest point for y-mesh line'

                    IF(IRRPT(I,J,K).EQ.0) THEN
                       NIRR = NIRR + 1
                       INDIRR(I,J,K)  = NIRR
                       IRRXYZ(1,NIRR) = I
                       IRRXYZ(2,NIRR) = J
                       IRRXYZ(3,NIRR) = K
                    END IF
                    IIRR = INDIRR(I,J,K)
                    IRRPT(I,J,K) =-1
                    CLOCAL(:,1,IIRR) = DPTSLINE(:,INDXSET)
                 END IF
              END DO
           ELSE
               LPTIO(:,J,K) =-1
           END IF
   
        END DO
     END DO

! fix the irregularity
     
     NIRROLD = NIRR
     DO NI=1,NIRROLD
        IX = IRRXYZ(1,NI)
        IY = IRRXYZ(2,NI)
        IZ = IRRXYZ(3,NI)

        IF(LPTIO(IX,IY,IZ).EQ.0) THEN 
           DO I=IX-1,IX+1
              DO J=IY-1,IY+1
                 DO K=IZ-1,IZ+1          
                    IF(IRRPT(I,J,K).EQ.0) THEN
                       NP  = LPTIO(I,J,K)
                       NP1 = LPTIO(I-1,J,K)*NP
                       NP2 = LPTIO(I+1,J,K)*NP
                       NP3 = LPTIO(I,J-1,K)*NP
                       NP4 = LPTIO(I,J+1,K)*NP
                       NP5 = LPTIO(I,J,K-1)*NP
                       NP6 = LPTIO(I,J,K+1)*NP              
           
                       NMIN = MIN(NP1,NP2,NP3,NP4,NP5,NP6)
                            
                       IF(NMIN.EQ.0.AND.NP.EQ.-1) THEN    
                          IRRPT(I,J,K)  = 1
                          NIRR = NIRR + 1
                          IRRXYZ(1,NIRR) = I
                          IRRXYZ(2,NIRR) = J
                          IRRXYZ(3,NIRR) = K
                          INDIRR(I,J,K)  = NIRR
                       END IF
                    END IF

                 END DO
              END DO
           END DO   
        END IF
     END DO

! mark the regular points

     NREG = 0
     DO I=1,NX
        DO J=1,NY
           DO K=1,NZ
              IF(IRRPT(I,J,K).EQ.0) THEN
                 NREG = NREG + 1    
                 INDREG(I,J,K) = NREG               
              END IF
           END DO
        END DO
     END DO

     write(*,'(a,i6)') '   Number of irregular points = ', nirr    

!---------------S.Y.---------------S.Y.---------------S.Y.---------------S.Y.---------------
! allocation for the grid points

! the representations of fictitious value and the associated 3 interface jumps

!     ALLOCATE(FICVA(20,NIRR), JPFIC(4,NIRR))
!     FICVZ = 0.D0
!     JPFIC = 0.D0

! the number and the position of grid points in the representation of fictitious value.
! the discretization state of irregular point
       
!     ALLOCATE(NXFIC(16,NIRR), NYFIC(16,NIRR), NZFIC(16,NIRR), NPFIC(NIRR), MDIS(9,NIRR))
!     NXFIC = 0
!     NYFIC = 0
!     NZFIC = 0
!     NPFIC = 0
!     MDIS = 0
     
!     ALLOCATE(AREG(7,NREG), AIRR(80,NIRR), NPAV(NIRR), IAV(80,NIRR), JAV(80,NIRR), KAV(80,NIRR))
!     AREG = 0.D0
!     AIRR = 0.D0
!     NPAV = 0
!     IAV  = 0
!     JAV  = 0
!     KAV  = 0
!---------------S.Y.---------------S.Y.---------------S.Y.---------------S.Y.---------------
END SUBROUTINE SETIRRPT2


SUBROUTINE PASSORTANGSID(THETA1,THETA2,IFACE1,IFACE2,SURFNORM,MXYZ)
     USE MOLECULE
     USE COMDATA
     IMPLICIT REAL*8(A-H,O-Z)
     REAL*8  :: SURFNORM(8,NFACE), THETA1, THETA2, VECLINE(3), VF1(3), VF2(3)
     INTEGER :: IFACE1, IFACE2, MXYZ
 
     SELECT CASE(MXYZ)
     CASE(1) 
        VECLINE(1) = 0.D0
        VECLINE(2) = 0.D0
        VECLINE(3) = 1.D0
     CASE(2)
        VECLINE(1) = 0.D0
        VECLINE(2) = 1.D0
        VECLINE(3) = 0.D0
     CASE(3)
        VECLINE(1) = 1.D0
        VECLINE(2) = 0.D0
        VECLINE(3) = 0.D0
     END SELECT

     VF1 = SURFNORM(6:8,IFACE1)
     VF2 = SURFNORM(6:8,IFACE2)

     THETA1 = SUM(VF1*VECLINE)
     THETA2 = SUM(VF2*VECLINE)
     
END SUBROUTINE PASSORTANGSID


SUBROUTINE PASSORTANGSPT(MULTP,NPOSPT,NFACESPT,SURFNORM,MXYZ, AP,BP,CP,DP,DPTB,PLANENORM,VPOS1,VPOS2,VPOS3,INDXF)
     USE MOLECULE
     USE COMDATA
     IMPLICIT REAL*8(A-H,O-Z)

     !COMMON/PLANE/AP,BP,CP,DP,DPTB,PLANENORM,VPOS1,VPOS2,VPOS3,INDXF
     
     REAL*8  :: AP, BP, CP, DP, DPTB, PLANENORM(3), VPOS1(3), VPOS2(3), VPOS3(3)
     REAL*8  :: VECLINE(3), SURFNORM(8,NFACE), DPTPLANE(30,2), PTTEST(3,2)
     INTEGER :: MULTP, NPOSPT, MXYZ, NFACESPT(31,NSPT), NV(3)
          
     SELECT CASE(MXYZ)
     CASE(1) 
        VECLINE(1) = 0.D0
        VECLINE(2) = 0.D0
        VECLINE(3) = 1.D0
     CASE(2)
        VECLINE(1) = 0.D0
        VECLINE(2) = 1.D0
        VECLINE(3) = 0.D0
     CASE(3)
        VECLINE(1) = 1.D0
        VECLINE(2) = 0.D0
        VECLINE(3) = 0.D0
     END SELECT
     
     NT =-NPOSPT
     PTTEST(:,1) = SPTPOS(:,NT) + 2.D-2*VECLINE
     PTTEST(:,2) = SPTPOS(:,NT) - 2.D-2*VECLINE
      
     NN = NFACESPT(1,NT)    
     
     DO I=1,NN 
        INDXF = NFACESPT(I+1,NT)
        AP    = SURFNORM(1,INDXF)
        BP    = SURFNORM(2,INDXF)
        CP    = SURFNORM(3,INDXF)
        DP    = SURFNORM(4,INDXF)
        DPTB  = SURFNORM(5,INDXF)
        PLANENORM = SURFNORM(6:8,INDXF)

        NV(1:3) = NVERT(:,INDXF)
        VPOS1(1:3) = SPTPOS(:,NV(1))
        VPOS2(1:3) = SPTPOS(:,NV(2))
        VPOS3(1:3) = SPTPOS(:,NV(3))

        DO NP=1,2 
           XX = PTTEST(1,NP) 
           YY = PTTEST(2,NP)
           ZZ = PTTEST(3,NP)
           !PRINT *,REAL(XX),REAL(YY),REAL(ZZ),NP,MXYZ
		   DPTPLANE(I,NP) = DPT(XX,YY,ZZ,NT,MXYZ, AP,BP,CP,DP,DPTB,PLANENORM,VPOS1,VPOS2,VPOS3,INDXF)
        END DO
     END DO
     
     DMIN1 = 1.D10
     DMIN2 = 1.D10
     DPT1  = 1.D10
     DPT2  = 1.D10

     DO I=1,NN

		IF(ABS(DPTPLANE(I,1)).LT.DMIN1) THEN
           DMIN1 = ABS(DPTPLANE(I,1))
           DPT1  = DPTPLANE(I,1)
        END IF

        IF(ABS(DPTPLANE(I,2)).LT.DMIN2) THEN
           DMIN2 = ABS(DPTPLANE(I,2))
           DPT2  = DPTPLANE(I,2)
        END IF
        
     END DO    

     IF(DPT1.GT.1.D9) THEN
        DPT1 =-SUM(VECLINE*SPTNRM(:,NT))
     END IF
	
     IF(DPT2.GT.1.D9) THEN
        DPT2 = SUM(VECLINE*SPTNRM(:,NT))
     END IF

     IF(DPT1*DPT2.LT.0.D0) THEN
        MULTP = 0
     ELSE
        MULTP = 1
     END IF     
 
END SUBROUTINE PASSORTANGSPT

SUBROUTINE SORTPTS(PTSLINE,DPTSLINE,PTMAX)
     USE COMDATA
     IMPLICIT REAL*8(A-H,O-Z)
     REAL*8  :: PTSLINE(31),DPTSLINE(4,30), DTMP(4), PTMAX
         
     NSET = INT(PTSLINE(1))
     IF(NSET.LE.1) RETURN
    
     DO I=1,NSET-1
        DMIN = PTSLINE(I+1)       
        MC   = 0
        DO II=I+1,NSET
           IF(PTSLINE(II+1).LE.DMIN) THEN
              DMIN = PTSLINE(II+1)
              IMIN = II
              MC   = 1
           END IF
        END DO
        IF(MC.EQ.1) THEN
           CTMP            = PTSLINE(IMIN+1)
           PTSLINE(IMIN+1) = PTSLINE(I+1)
           PTSLINE(I+1)    = CTMP
           
           DTMP             = DPTSLINE(:,IMIN)
           DPTSLINE(:,IMIN) = DPTSLINE(:,I)
           DPTSLINE(:,I)    = DTMP                      
        END IF
     END DO
   
     DO I=NSET,1,-1
        IF(PTSLINE(I+1).GT.1.D9) THEN 
           PTSLINE(1) = PTSLINE(1) - 1
        ELSEIF(PTSLINE(I+1).LE.PTMAX) THEN
           GOTO 100       
        END IF
     END DO

100  RETURN

END SUBROUTINE SORTPTS


FUNCTION DPT(XX,YY,ZZ,NPOSPT,MXYZ, AP,BP,CP,DP,DPTB,PLANENORM,VPOS1,VPOS2,VPOS3,INDXF)
     USE MOLECULE
     USE COMDATA
     IMPLICIT REAL*8(A-H,O-Z)
     !COMMON/PLANE/AP,BP,CP,DP,DPTB,PLANENORM,VPOS1,VPOS2,VPOS3,INDXF

     REAL*8  :: VPOS1(3), VPOS2(3), VPOS3(3), PLANENORM(3), XX, YY, ZZ, DPT
     REAL*8  :: SIDE1(3), SIDE2(3), SIDE3(3), XYZP(3), V1(3), V2(3), PTPROJ(3)
     
     PI    = ACOS(-1.D0)
     DCEL9 = 9.D0*DCEL
    
     VMAG =-(AP*XX + BP*YY + CP*ZZ + DP)/(AP*PLANENORM(1) + BP*PLANENORM(2) + CP*PLANENORM(3))

	 XYZP(1) = XX + VMAG*PLANENORM(1)
     XYZP(2) = YY + VMAG*PLANENORM(2)
     XYZP(3) = ZZ + VMAG*PLANENORM(3)
     
     SIDE1 = VPOS1 - XYZP
     SIDE2 = VPOS2 - XYZP
     SIDE3 = VPOS3 - XYZP  

     SIDE1L = DSQRT(SUM(SIDE1**2))
     SIDE2L = DSQRT(SUM(SIDE2**2))
     SIDE3L = DSQRT(SUM(SIDE3**2))  
          
!     IF(SIDE1L+SIDE2L+SIDE3L.GT.DCEL9)  THEN
!        DPT = 1.D10
!        RETURN
!     END IF

! the projection of point conincident with a vertex

     IF(SIDE1L.LE.1.d-6) THEN
        DPT = SIGN(ABS(AP*XX + BP*YY + CP*ZZ + DP)*DPTB,VMAG)
        RETURN
     ELSEIF(SIDE2L.LE.1.d-6) THEN
        DPT = SIGN(ABS(AP*XX + BP*YY + CP*ZZ + DP)*DPTB,VMAG)
        RETURN
     ELSEIF(SIDE3L.LE.1.d-6) THEN
        DPT = SIGN(ABS(AP*XX + BP*YY + CP*ZZ + DP)*DPTB,VMAG)
        RETURN                      
     END IF

! the projection of point is inside the triangle

     THETA1 = SUM(SIDE1*SIDE2)/(SIDE1L*SIDE2L)
     THETA2 = SUM(SIDE2*SIDE3)/(SIDE2L*SIDE3L)
     THETA3 = SUM(SIDE3*SIDE1)/(SIDE3L*SIDE1L)

     IF(ABS(THETA1).GT.1.D0) THETA1 = THETA1/ABS(THETA1)
     IF(ABS(THETA2).GT.1.D0) THETA2 = THETA2/ABS(THETA2)
     IF(ABS(THETA3).GT.1.D0) THETA3 = THETA3/ABS(THETA3)
   
     THETA1 = DACOS(THETA1)
     THETA2 = DACOS(THETA2)
     THETA3 = DACOS(THETA3)
   
     IF(ABS(THETA1-PI).LE.1.D-2) THEN ! the projection of point is on side 1
		V1(1) = XX - VPOS2(1)
        V1(2) = YY - VPOS2(2)
        V1(3) = ZZ - VPOS2(3)

        V1NM = SQRT(SUM(V1**2))
        V1   = V1/V1NM

        V2   = VPOS1 - VPOS2
        V2NM = SQRT(SUM(V2**2))
        V2   = V2/V2NM

        CTHETA  = SUM(V1*V2)				!BUG OF DABAO
        PTPROJ  = VPOS2 + V2*V1NM*CTHETA
        STHETA  = SQRT(1.D0 - CTHETA**2)
    
        DPT   = STHETA*V1NM
        PROJ  = SUM((V1-PTPROJ)*SPTNRM(:,NPOSPT))
        DPT   = SIGN(DPT,PROJ)
     ELSEIF(ABS(THETA2-PI).LE.1.D-2) THEN   ! the projection of point is on side 2
 		V1(1)   = XX - VPOS3(1)
        V1(2)   = YY - VPOS3(2)
        V1(3)   = ZZ - VPOS3(3)

        V1NM = SQRT(SUM(V1**2))
        V1   = V1/V1NM

        V2   = VPOS2 - VPOS3
        V2NM = SQRT(SUM(V2**2))
        V2   = V2/V2NM

        CTHETA  = SUM(V1*V2)
        PTPROJ  = VPOS3 + V2*V1NM*CTHETA
        STHETA  = SQRT(1.D0 - CTHETA**2)
    
        DPT   = STHETA*V1NM
        PROJ  = SUM((V1-PTPROJ)*SPTNRM(:,NPOSPT))
        DPT   = SIGN(DPT,PROJ)
     ELSEIF(ABS(THETA3-PI).LE.1.D-2) THEN   ! the projection of point is on side 3 
 		V1(1)   = XX - VPOS1(1)
        V1(2)   = YY - VPOS1(2)
        V1(3)   = ZZ - VPOS1(3)

        V1NM = SQRT(SUM(V1**2))
        V1   = V1/V1NM

        V2   = VPOS3 - VPOS1
        V2NM = SQRT(SUM(V2**2))
        V2   = V2/V2NM

        CTHETA = SUM(V1*V2)
        PTPROJ = VPOS1 + V2*V1NM*CTHETA
        STHETA = SQRT(1.D0 - CTHETA**2)
    
        DPT   = STHETA*V1NM
        PROJ  = SUM((V1-PTPROJ)*SPTNRM(:,NPOSPT))
        DPT   = SIGN(DPT,PROJ)   
     ElSEIF(ABS(THETA1 + THETA2 + THETA3 - 2.D0*PI).LT.1.d-6) THEN
		DPT = SIGN(ABS(AP*XX + BP*YY + CP*ZZ + DP)*DPTB,VMAG)
     ELSE
        DPT = 1.D10
     END IF

END FUNCTION DPT

SUBROUTINE EQPLANE(A,B,C,D,V1,V2,V3)
     REAL*8 :: A, B, C, D, V1(3), V2(3), V3(3)
     
     A = V2(2)*V3(3) + V3(2)*V1(3) + V1(2)*V2(3) - V3(2)*V2(3) - V2(2)*V1(3) - V1(2)*V3(3)
     B = V1(1)*V3(3) + V2(1)*V1(3) + V3(1)*V2(3) - V3(1)*V1(3) - V2(1)*V3(3) - V1(1)*V2(3)
     C = V1(1)*V2(2) + V3(1)*V1(2) + V2(1)*V3(2) - V2(2)*V3(1) - V1(1)*V3(2) - V1(2)*V2(1)
     D =-( V1(1)*V2(2)*V3(3) + V2(1)*V3(2)*V1(3) + V3(1)*V1(2)*V2(3) - &
           V3(1)*V2(2)*V1(3) - V1(1)*V3(2)*V2(3) - V1(2)*V2(1)*V3(3) )
     
END SUBROUTINE EQPLANE
      

SUBROUTINE CIRCUMCENTER(XYZCIRC,CIRCM)
     REAL*8 :: XYZCIRC(3), CIRCM(3,4)

END SUBROUTINE CIRCUMCENTER

SUBROUTINE NMPLANE(VNORM,V1,V2,V3)
     IMPLICIT REAL*8(A-H,O-Z)
     REAL*8 :: VNORM(3), V1(3), V2(3), V3(3), S1(3), S2(3), VNORML
     
     S1 = V3 - V1
     S2 = V2 - V1

     VNORM(1) = S1(2)*S2(3) - S1(3)*S2(2)
     VNORM(2) = S1(3)*S2(1) - S1(1)*S2(3)
     VNORM(3) = S1(1)*S2(2) - S1(2)*S2(1)

     VNORML = SQRT(SUM(VNORM**2))
     if(vnorml<1.d-6)then
	     vnorm(1:3)=(/1,0,0/)
	 else
		 VNORM  = VNORM/VNORML
	 endif

END SUBROUTINE NMPLANE

FUNCTION DETERMINANT3(V1,V2,V3)
     IMPLICIT REAL*8(A-H,O-Z)

     REAL*8 :: DETERMINANT3, V1(3), V2(3), V3(3)
     
     DETERMINANT3 = V1(1)*V2(2)*V3(3) + V1(2)*V2(3)*V3(1) + V1(3)*V2(1)*V3(2) - &
                    V1(3)*V2(2)*V3(1) - V1(1)*V2(3)*V3(2) - V1(2)*V2(1)*V3(3)
     
END FUNCTION DETERMINANT3
