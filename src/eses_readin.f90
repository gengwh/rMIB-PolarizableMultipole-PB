!Read protein structure from eses
subroutine eses_readin
    use pbedata
    use molecule
    use comdata
    implicit double precision(a-h,o-z)

    ! temp. var. for read in
CHARACTER(100) :: FHEAD
character(10) :: c1,c3,c4,c5,c6,c8
integer i7,i5, max_ires, max_iatm
integer i,j,iflag,nremark,MEOF, nremark1, nremark2, natm1, natm2, nchr2, nface2, nspt1, nspt2
real*8 xyzqr(5)
! real*8 xyzqr(14)

    
    !For ESES
    real*8, dimension(3) :: INTERSECTION, NORMAL
    real*8, dimension(4) :: clocal_temp
    integer npts,inx,iny,inz,outx,outy,outz,rc,ix,iy,iz,pointio
    integer max_num_irregular, intersectionCount
    CHARACTER(7) :: dcel_STRING
    REAL*8 extension
    CHARACTER(8) :: extension_STRING
    CHARACTER(39):: box_STRING
    CHARACTER(100)::intersection_STRING,grid_STRING
    CHARACTER(9):: surfaceLocation_STRING
    LOGICAL:: box_exists, intersection_exists, grid_exists, always_calculate
    real*8 :: phi, theta
    real*8,dimension(:,:,:), allocatable:: irrpt_count
    integer, dimension(3) :: idx
    integer :: i2,j2,k2, old_nirr

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
do i=1,nremark
    read(102,*) fhead
enddo

natm=0
do
    read(102,*,IOSTAT = MEOF) c1,i2,c3,c4,i5,xyzqr 
    if ((c1(1:3) .ne. 'END') .and. (MEOF .eq. 0)) then
        !write(103,*) real(xyzqr(1:3)),real(xyzqr(5))
        natm=natm+1
    endif
    IF(MEOF .LT. 0) EXIT
enddo
close(102)
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
!     quad_mom(1,:,i)=(/xyzqr(8), xyzqr(9), xyzqr(11)/)
!     quad_mom(2,:,i)=(/xyzqr(9), xyzqr(10), xyzqr(12)/)
!     quad_mom(3,:,i)=(/xyzqr(11), xyzqr(12), xyzqr(13)/)

enddo
close(102)

    edge=2.d0*radw
    chrpos = atmpos
    
    !the locations of ESES surface outputs
    write(dcel_STRING,'(f7.5)') dcel
    surfaceLocation_STRING='surfaces/'
    box_STRING=surfaceLocation_STRING//fname(1:lenfname)//'-h'//dcel_STRING//'bounding_box.txt'
    intersection_STRING=surfaceLocation_STRING//fname(1:lenfname)//'-h'//dcel_STRING//'intersection_info.txt'
    grid_STRING=surfaceLocation_STRING//fname(1:lenfname)//'-h'//dcel_STRING//'grid_info.txt'


    write(*,*)"box_STRING=["//box_STRING//"]"
    !check if this ESES  surface has already been calculated
    INQUIRE(FILE=box_STRING,EXIST=box_exists)
    INQUIRE(FILE=intersection_STRING,EXIST=intersection_exists)
    INQUIRE(FILE=grid_STRING,EXIST=grid_exists)

    always_calculate = .TRUE.

    !if any part of this ESES surface has not been calculated, then calculate it
    IF (always_calculate .OR. (.NOT. box_exists) .OR.(.NOT.  intersection_exists) .OR. (.NOT. grid_exists)) THEN

        !extension = 3.5d0 !topology.f90 uses nlay*dcel1=max(7,int(1.d0/dcel1)) where dcel1=0.5d0 constant, so always 3.5d0
        extension = 3.0d0 !Siwen's rMIB using MSMS generates a [-5x5]^3 domain for 1ato (rad=2) !yang: 3.0 -> 1.0
        write(extension_STRING,'(f7.5)') extension

        !calculate surface
        write(dcel_STRING,'(f7.5)') dcel
        write(*,*)"calling ESES"    !yang: eses -> MS_Intersection; .xyzrq for eses -> .xyzr for MS_Intersection
        write(*,*)'MS_Intersection '//pathname(1:lenpath)//fname(1:lenfname)//'.pqr 1.4 '//dcel_STRING//' '//extension_STRING
        call system('MS_Intersection '//pathname(1:lenpath)//fname(1:lenfname)//'.pqr 1.4 '//dcel_STRING//' '//extension_STRING)

        !move outputs to surfaces folder
        write(*,*)"Finished calling ESES. Move ESES outputs to surfaces folder..."
        print*,"[cp bounding_box.txt "//box_STRING//"]"
        call system('mv \./bounding_box.txt '//box_STRING)
        call system('mv \./intersection_info.txt '//intersection_STRING)
        call system('mv \./grid_info.txt '//grid_STRING)
        write(*,*)"Finished moving ESES outputs to surfaces folder."
    ELSE
        !the whole surface exists, so do not waste time
        write(*,*)"ESES surface files already exist. Using those without recalculating."
    END IF
    
    write(*,*)"Reading in bounding box from "//box_STRING
    open(2,FILE=box_STRING)
    READ(2,*) xleft,yleft,zleft
    READ(2,*) xright,yright,zright
    READ(2,*) nx,ny,nz

    NTOT = nx*ny*nz
    npts = NTOT
    close(2)
    allocate(x(nx),y(ny),z(nz))
    
    print*, "xleft = ", xleft, "xright= ", xright, "nx= ",nx
    print*, "yleft = ", yleft, "yright= ", yright, "ny= ",ny
    print*, "zleft = ", zleft, "zright= ", zright, "nz= ",nz
    print*,"============"

    dx=dcel
    dy=dcel
    dz=dcel
    print*, dx, dy, dz
    do i=1,nx
        x(i)=xleft+(i-1)*dx
    enddo
    do i=1,ny
        y(i)=yleft+(i-1)*dy
    enddo
    do i=1,nz
        z(i)=zleft+(i-1)*dz
    enddo

    ! Read in grid inside/outside info
    allocate(io(nx,ny,nz))
    write(*,*)"Reading in grid info from "//grid_STRING
    open(2,FILE=grid_STRING)
    do i=1,npts
        READ(2,*)ix,iy,iz,pointio
        io(ix+1,iy+1,iz+1)=pointio
    end do
    close(2)
    io = io*(-1) !ESES is backwards

    !get number of intersections
    nirr = 0
    intersectionCount = 0

    write(*,*)"Reading in intersection info from "//intersection_STRING

    allocate(irrpt_count(nx,ny,nz),stat=ierr)
    if (ierr .ne. 0) then
        write(*,*) "error allocating irrpt_count"
        stop
    endif
    irrpt_count = 0

    open(2,FILE=intersection_STRING)
    do
        read(2,*,IOSTAT=rc) inx, iny, inz, outx, outy, outz, FHEAD
        if (rc /= 0 ) then
            EXIT
        end if
        irrpt_count(inx,iny,inz) = 1
        irrpt_count(outx,outy,outz) = 1
        intersectionCount = intersectionCount + 1
    end do
    close(2)

    do k = 1,nz
        do j = 1,ny
            do i = 1,nx
                nirr = nirr + irrpt_count(i,j,k)
            enddo
        enddo
    enddo
    deallocate(irrpt_count)

    print*,"preliminary count of irregular points: ",nirr
    max_num_irregular = nirr*2
    nirr = 0

    allocate(dxl(max_num_irregular),dxr(max_num_irregular),dyl(max_num_irregular),stat = ierr)   
    allocate(dyr(max_num_irregular),dzl(max_num_irregular),dzr(max_num_irregular),stat = ierr)
    allocate(mcx(max_num_irregular),mcy(max_num_irregular),mcz(max_num_irregular),stat = ierr)
    allocate(clocal(4,3,max_num_irregular))
    if (ierr .ne. 0) then
        write(*,*) "Error allocating exdl,exdr,dyl,dyr,dzl,dzr"
        stop
    endif

    mcx = 0; mcy = 0; mcz = 0

    do i=1,max_num_irregular
        dxl(i) = 0.d0
        dyl(i) = 0.d0
        dzl(i) = 0.d0
        dxr(i) = 0.d0
        dyr(i) = 0.d0
        dzr(i) = 0.d0
    enddo
    
    allocate(irrpt(nx,ny,nz),indirr(nx,ny,nz),irrxyz(3,max_num_irregular),stat = ierr)
    if (ierr .ne. 0) then   
        write(*,*) "Error allocating irrpt,indirr, irrxyz" 
        stop
    endif

    irrpt = 0
    indirr = 0
    irrxyz = 0
    open(2,FILE=intersection_STRING)
    do i=1,intersectionCount !nirr/2 lines in file, need to index
        READ(2,*) inx,iny,inz,outx,outy,outz, INTERSECTION, NORMAL
        !ESES uses 1-based indexing!
        inx = inx + 1
        iny = iny + 1
        inz = inz + 1
        outx = outx + 1
        outy = outy + 1
        outz = outz + 1
        !ESES uses 1-based indexing!

        

        ! set CLOCAL angles. The comment describing the CLOCAL format in msmsibw13.f90 disagrees with
        ! the rMIB usage of the values. The correct implimentation, per August 2021 correspondence with Weihua Geng
        ! and referring to "Treatment of charge singularities in implicit solvent models" by Geng, Yu, and Wei (2007),
        ! is clocal_temp(1:4) = (sin(phi), cos(phi), cos(theta), sin(theta)), where
        !   theta: angle between x-axis and the projection of the normal vector onto the x-y plane
        !   phi: angle between z-axis and normal vector

        ! cos(phi) = z
        clocal_temp(2) = NORMAL(3)  
        
        ! sin(phi) = opposite/hypotenuse = sqrt(x^2+y^2)/sqrt(x^2+y^2+z^2),
        ! but sqrt(x^2+y^2+z^2)  = magnitude of normal = 1 
        clocal_temp(1) = sqrt(NORMAL(1)**2 + NORMAL(2)**2)  
        
        IF (clocal_temp(1) .GE. 1d-6) THEN
            ! cos(theta) = x/sqrt(x^2+y^2) = x/sin(phi) in this special case
            clocal_temp(3) = NORMAL(1)/clocal_temp(1)
            ! sin(theta) = y/sqrt(x^2+y^2) = y/sin(phi) in this special case
            clocal_temp(4) = NORMAL(2)/clocal_temp(1)        
        ELSE
            ! sqrt(x^2+y^2) ~ 0, so use default values to avoid NaNs and Infs
            clocal_temp(3) = 1
            clocal_temp(4) = 0            
        ENDIF


        if (irrpt(inx,iny,inz) .eq. 0) then
            call add_irregular(inx,iny,inz)
        end if
        if (irrpt(outx,outy,outz) .eq. 0) then
            call add_irregular(outx,outy,outz)
        end if

!modified to fit msmsmibw13.f90
        !x irregularity
        if (inx < outx) then
            MCX(indirr(inx,iny,inz)) = 1
            MCX(indirr(outx,outy,outz))=  -1
            dxl(indirr(inx,iny,inz)) = INTERSECTION(1) - x(inx)
            dxr(indirr(inx,iny,inz)) = x(outx) - INTERSECTION(1)
            clocal(1:4,1,indirr(inx,iny,inz)) = clocal_temp
        else if (inx > outx) then   
            MCX(indirr(inx,iny,inz)) = -1
            MCX(indirr(outx,outy,outz))= 1
            dxl(indirr(outx,outy,outz)) = INTERSECTION(1) - x(outx)
            dxr(indirr(outx,outy,outz)) = x(inx) - INTERSECTION(1)
            clocal(1:4,1,indirr(outx,outy,outz)) = clocal_temp
        end if
    !y irregularity
        if (iny < outy) then   
            MCY(indirr(inx,iny,inz)) = 1
            MCY(indirr(outx,outy,outz))= -1
            dyl(indirr(inx,iny,inz)) = INTERSECTION(2) - y(iny)
            dyr(indirr(inx,iny,inz)) = y(outy) - INTERSECTION(2)
            clocal(1:4,2,indirr(inx,iny,inz)) = clocal_temp
        else if (iny > outy) then     
            MCY(indirr(inx,iny,inz)) = -1
            MCY(indirr(outx,outy,outz))= 1
            dyl(indirr(outx,outy,outz)) = INTERSECTION(2) - y(outy)
            dyr(indirr(outx,outy,outz)) = y(iny) - INTERSECTION(2)
            clocal(1:4,2,indirr(outx,outy,outz)) = clocal_temp
        end if
    !z irregularity
        if (inz < outz) then
            MCZ(indirr(inx,iny,inz)) = 1
            MCZ(indirr(outx,outy,outz))= -1
            dzl(indirr(inx,iny,inz)) = INTERSECTION(3) - z(inz)
            dzr(indirr(inx,iny,inz)) = z(outz) - INTERSECTION(3)
            clocal(1:4,3,indirr(inx,iny,inz)) = clocal_temp
        else if (inz > outz)  then
            MCZ(indirr(inx,iny,inz)) = -1
            MCZ(indirr(outx,outy,outz))= 1
            dzl(indirr(outx,outy,outz)) =  INTERSECTION(3) - z(outz)
            dzr(indirr(outx,outy,outz)) = z(inz) - INTERSECTION(3)  
            clocal(1:4,3,indirr(outx,outy,outz)) = clocal_temp 
        end if! inz neq outz
    enddo
    close(2)
    write(*,*)"Finished reading ",nirr,"irregular points. max_num_irregular=",max_num_irregular

    old_nirr = nirr
    
    ! MSMS criteria for on interface is if distance in 1 direction is < 1e-6
    do irr = 1, old_nirr
        if ((mcx(irr) .eq. 1 .and. abs(dxl(irr)) .lt. 1e-6) &
            .or. (mcy(irr) .eq.  1 .and. abs(dyl(irr)) .lt. 1e-6) &
            .or. (mcz(irr) .eq.  1 .and. abs(dzl(irr)) .lt. 1e-6)) then
             !interface points have irrpt=-1, io=0
            i1 = irrxyz(1,irr)
            j1 = irrxyz(2,irr)
            k1 = irrxyz(3,irr)
            irrpt(i1,j1,k1) = -1
            io(i1,j1,k1) = 0
           
            do ik=1,6  
            ! scan for inside neighbors of the interface point that should be marked irregular but aren't yet
                call idx6(i1,j1,k1,ik,idx)
                i2 = idx(1)
                j2 = idx(2)
                k2 = idx(3)
                if (io(i2,j2,k2) .eq. -1 .and. irrpt(i2,j2,k2) .eq. 0) then
                    call add_irregular(i2,j2,k2)
                endif

            enddo
            
        endif
    enddo


end subroutine eses_readin

subroutine add_irregular(i,j,k)
    use comdata
    irrpt(i,j,k) = 1
    nirr = nirr + 1
    indirr(i,j,k) = nirr
    irrxyz(:,nirr) = (/i,j,k/)
end subroutine add_irregular
