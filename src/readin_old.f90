!Read protein structure from msms
subroutine readin
use pbedata
use molecule
use comdata
implicit double precision(a-h,o-z)
real*8 pos(3),vector(3)
integer nind(5)
! temp. var. for read in
CHARACTER(100) :: FHEAD
character(10) :: c1,c3,c4,c5,c6,c8
integer i2,i7,i5, max_ires, max_iatm
integer i,j,iflag,nremark,MEOF, nremark1, nremark2, natm1, natm2, nchr2, nface2, nspt1, nspt2
real*8 xyzqr(5)
! real*8 xyzqr(14)


!Obtain path
pathname='test_proteins/'
lenpath = len(pathname)
do while (pathname(lenpath:lenpath) .eq. ' ')
    lenpath = lenpath - 1
enddo  

!Obtain filename
lenfname = len(fname)
do while (fname(lenfname:lenfname) .eq. ' ')
    lenfname = lenfname - 1
enddo 

nremark=0
open(102,file=pathname(1:lenpath)//fname(1:lenfname)//".pqr")
do
    READ(102,*) fhead
    if (fhead(1:6)=='REMARK') then
        nremark=nremark+1
    else
        exit
    endif
enddo
print *,'lines of remarks = ', nremark
close(102)

open(102,file=pathname(1:lenpath)//fname(1:lenfname)//".pqr")
open(103,file=pathname(1:lenpath)//fname(1:lenfname)//".xyzr")

do i=1,nremark
    read(102,*) fhead
enddo

natm=0
do
    read(102,*,IOSTAT = MEOF) c1,i2,c3,c4,i5,xyzqr 

    if ((c1(1:3) .ne. 'END') .and. (MEOF .eq. 0)) then
        write(103,*) real(xyzqr(1:3)),real(xyzqr(5))
!         write(103,*) real(xyzqr(1:3)),real(xyzqr(14))
        natm=natm+1
    endif
    IF(MEOF .LT. 0) EXIT
enddo

close(102)
close(103)
print *,'number of atoms = ', natm

!Read atom coordinates and partial charges
 nchr=natm
 allocate(atmpos(3,natm),atmrad(natm),atmchr(nchr),chrpos(3,nchr),STAT=ierr)
IF (ierr .NE. 0) THEN
    WRITE(6,*) 'Error allocating atmpos, atmrad, atmchr, chrpos!'
    STOP
END IF

open(102,file=pathname(1:lenpath)//fname(1:lenfname)//".pqr")

do i=1,nremark
    read(102,*) fhead
enddo

! allocate(di_mom(3,natm),quad_mom(3,3,natm)) !yang
do i=1,natm
    read(102,*,IOSTAT = MEOF) c1,i2,c3,c4,i5,xyzqr 
    atmpos(1:3,i)=xyzqr(1:3)
    atmrad(i)=xyzqr(5)
    chrpos(1:3,i)=xyzqr(1:3)
    atmchr(i)=xyzqr(4)
    
!     atmpos(1:3,i)=xyzqr(1:3)
!     atmrad(i)=xyzqr(14)
!     chrpos(1:3,i)=xyzqr(1:3)
!     atmchr(i)=xyzqr(4)
!     di_mom(:,i)=xyzqr(5:7)
!     quad_mom(:,1,i)=(/xyzqr(8), xyzqr(9), xyzqr(11)/)
!     quad_mom(:,2,i)=(/xyzqr(9), xyzqr(10), xyzqr(12)/)
!     quad_mom(:,3,i)=(/xyzqr(11), xyzqr(12), xyzqr(13)/)

!     print *,i
!     print *,atmchr(i)
!     print *,di_mom(:,i)
!     print *,real(quad_mom(:,:,i))
enddo
close(102)

rslt=system('msms -if '//pathname(1:lenpath)//fname(1:lenfname)//".xyzr"//' -prob 1.4 -de ' &
//den(1:5)//' -of '//pathname(1:lenpath)//fname(1:lenfname)//' > msms.out')    

!rslt=system('rm '//pathname(1:lenpath)//fname(1:lenfname)//".vert")

end

! read the surface points
!OPEN(102,FILE=pathname(1:lenpath)//FNAME(1:lenfname)//".vert")

!READ(102,*) FHEAD
!READ(102,*) FHEAD
!READ(102,*) NSPT, ppp, qqq, rrr

    !ALLOCATE(SPTPOS(3,NSPT), SPTNRM(3,NSPT), NATMAFF(NSPT), NSFTYPE(NSPT), STAT= ierr)
    !IF (ierr .NE. 0) THEN
    !    WRITE(6,*) 'Error allocating SPTPOS, SPTNRM, NATMAFF, NSFTYPE!'
    !STOP
    !END IF

    !SPTPOS=0.D0; SPTNRM=0.D0; NATMAFF=0; NSFTYPE=0;	
      
    !DO I=1,NSPT
    !    READ(102,*) POS(1:3), VECTOR(1:3), KK, NAFF, NAFFT 
         
    !    SPTPOS(:,I) = POS;   SPTNRM(:,I) = VECTOR
    !    NATMAFF(I)  = NAFF;  NSFTYPE(I)  = NAFFT
    !END DO
    !CLOSE(102)

! read the surface triangulization

    !OPEN(103,FILE=pathname(1:lenpath)//FNAME(1:lenfname)//".face")

    !READ(103,*) FHEAD
    !READ(103,*) FHEAD
    !READ(103,*) NFACE, PPP, QQQ, RRR

    !ALLOCATE(NVERT(3,NFACE), MFACE(NFACE), STAT=ierr)
    !IF (ierr .NE. 0) THEN
    !    WRITE(6,*) 'Error allocating NVERT, MFACE'
    !STOP
    !END IF

    !NVERT=0; MFACE=0
      
    !DO I=1,NFACE 
    !    READ(103,*) NIND(1:5) 
    !    NVERT(1:3,I) = NIND(1:3);  MFACE(I) = NIND(4)
    !END DO
    !CLOSE(103)
    !call surface_area(s_area) ! the post-MSMS code
    !print *,'surface area=', real(s_area)

    !rslt=system('rm '//pathname(1:lenpath)//fname(1:lenfname)//".vert")
    !rslt=system('rm '//pathname(1:lenpath)//fname(1:lenfname)//".face")


!#############################################################################################
!The face file contains three header lines followed by one triangle per line. 
!The first header line provides a comment and the file name of the sphere set. 
!The second header line holds comments about the content of the third line. 
!The third header line provides the number of triangles, the number of spheres in the set, 
!the triangulation density and the probe sphere radius. 

!The first three numbers are (1 based) vertex indices. 

!The next field can be: 
!1 for a triangle in a toric reen trant face, 
!2 for a triangle in a spheric reentrant face and 
!3 for a triangle in a contact face. 

!The last # on the line is the (1 based) face number in the analytical description of the solvent excluded surface. 
!These values are written in the following format ``%6d %6d %6d %2d %6d''.

!The vertex file contains three header lines (similar to the header in the .face file) 
!followed by one vertex per line and provides the coordinates (x,y,z) and the normals (nx,ny,nz) 
!followed by the number of the face (in the analytical description of the solvent excluded surface) 
!to which the vertex belongs. The vertices of the analytical surface have a value 0 in that field 
!and the vertices lying on edges of this surface have negative values. 
!The next field holds the (1 based) index of the closest sphere. 
!The next field is 
!1 for vertices which belong to toric reentrant faces (including ver tices of the analytical surface), 
!2 for vertices inside reentrant faces and 
!3 for vertices inside contact faces. 
!These values are written in the following format ``%9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %7d %7d %2d''.
