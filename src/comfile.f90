MODULE MOLECULE
    integer, parameter :: max_natm = 10000, max_nspt = 50000, max_nface=120000, mapt = 100
    real*8, parameter  :: radw = 1.4d0
    integer :: natm, nspt, nsatm, nface, nchr, nt
    
    real*8,  dimension(:),     allocatable :: atmrad, atmchr
    real*8,  dimension(:,:),   allocatable :: atmpos, sptpos, sptnrm, curv,qxyz, qdxyz, atmpos_sph, charget
    real*8,  dimension(:,:),   allocatable :: chrpos, chrpos_sph, atmmch !yang atmmch=（atmchr+3di_mom+6quad_mom)*natm
    real*8,  dimension(:,:),   allocatable :: di_mom  !yang 3*natm
    real*8,  dimension(:,:,:), allocatable :: quad_mom !yang 3*3*natm
    real*8,  dimension(:,:,:), allocatable :: atrans, chgmnx, corlocqt
   
    integer, dimension(:),     allocatable :: natmsf, nsfatm, natmaff, nsftype, npture, nsptno, mface
    integer, dimension(:,:),   allocatable :: neigh, nsptatm    
    integer, dimension(:,:),   allocatable :: nvert
    integer, dimension(:,:,:), allocatable :: iqxyz, loc_qt 

END MODULE MOLECULE

MODULE COMDATA   

! atom data file   
   character(100) :: fname,pathname,den

! grid numbers in x, y and z direction

   integer :: nx,ny,nz,ntot,max_nirr,nirr,nreg,nelt,nsrfdrn,iratio,lenpath,lenfname

   real*8, dimension(:), allocatable :: xc, yc, zc
   real*8 :: xleft, xright, yleft, yright, zleft, zright, dx, dy, dz, dcel
   real*8 :: rdx1, rdx2, rdy1, rdy2, rdz1, rdz2, cmain  
   real*8 :: cent2(-1:1)
   real*8,  dimension(:),     allocatable :: dxl, dxr, dyl, dyr, dzl, dzr  
                        
   real*8,  dimension(:,:,:), allocatable :: clocal
   integer, dimension(:),     allocatable :: mcx, mcy, mcz
   integer, dimension(:,:),   allocatable :: IRRXYZ     
   integer, dimension(:,:,:), allocatable :: lptio, irrpt, indirr, indreg
            
END MODULE COMDATA 


Module PBEDATA
INTEGER nirrall, indsrf, ims, nfc, icg, isf, ibd, inl, ipm
Integer nf, ls, itrbl, iph,itnew	! Sining's pactch

REAL*8	eps0, eps1, eps, rds, prds, ra, y0, err_flux, cappa, cappa_nl
real*8  xa,xb !for two balls
     
REAL*8, DIMENSION(:),			ALLOCATABLE :: x, y, z, uirr, ueirr, b, u, uexa, bftc, uapbs, uirrapbs, xt,yt,zt 
REAL*8, DIMENSION(:),			ALLOCATABLE :: cftc, srfftc, u0, bphi0, phi0, phi0e 
REAL*8, DIMENSION(:),			allocatable :: phi_s, phi_v, phi_rxn !added for mibpb3
REAL*8, DIMENSION(:,:),			ALLOCATABLE :: srfpts,ftc, onsrf, rrintf, pchftc, ugrad !yang
REAL*8, DIMENSION(:,:,:),		ALLOCATABLE :: ftpts,dist_angle,akpatch, phi0e_d, ugradgrad !yang

INTEGER,DIMENSION(:),			ALLOCATABLE :: idpatch
INTEGER,DIMENSION(:,:),			ALLOCATABLE :: irrdrc,isrfpts,isrf,iftwrg, idre, itps, irrintf
INTEGER,DIMENSION(:,:),			ALLOCATABLE :: irrphi0,itpsnew, idpatchreal
INTEGER,DIMENSION(:,:,:),		ALLOCATABLE :: io, irrpts, iftpts, ioe,  icpatch

END MODULE

module ModInput

integer :: idAcSl, idSurf, igen, iFourth

end module

module bicg
integer nmax,nsize
real*8, dimension(:), allocatable:: sa,sb,u00,utemp,b00
integer, dimension(:), allocatable:: ijka,ijkb
end module bicg



