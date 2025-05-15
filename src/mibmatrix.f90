!---------------------------------------------------------------------------------------
Subroutine formmx
!the coodinates of (i,j,k) is stored in sa(ijka())
use comdata
use pbedata
use bicg
implicit double precision(A-H, O-Z)
integer*4 n, nsign,iaa(100)
dimension bt(1:6),betax(6),betay(6),betaz(6),aa(100)
integer	  ijkpt,ijk,idx(3)
real*8 cappa2 !squre of D-H parameter

sa=0.d0; ijka=0
itest=0
ijka(1)=nx*ny*nz+2				!matrix size plus 2
ijkpt=nx*ny*nz+1				!ijkpt trace the how many elements filled in the two biconjugate vectors
do k=1,nz
do j=1,ny
do i=1,nx
	ijk=(k-1)*nx*ny+(j-1)*nx+i

	ijka(ijk+1)=ijka(ijk) !cumulated non-zero elements so far
	
	if ((i==1) .or. (i==nx) .or. (j==1) .or. (j==ny) .or. (k==1) .or. (k==nz)) then
		sa(ijk)=1.d0
		!bftc(ijk)=upi(i,j,k)
		!bftc(ijk)=bdcond(x(i),y(j),z(k))
		!bftc(ijk)=uexa(id3d(i,j,k))
	else
		betax(1)=x(i);					betay(1)=y(j);					betaz(1)=(z(k)+z(k-1))/2.d0;
		betax(2)=x(i);					betay(2)=(y(j)+y(j-1))/2.d0;	betaz(2)=z(k);
		betax(3)=(x(i)+x(i-1))/2.d0;	betay(3)=y(j);					betaz(3)=z(k);
		betax(4)=(x(i)+x(i+1))/2.d0;	betay(4)=y(j);					betaz(4)=z(k);
		betax(5)=x(i);					betay(5)=(y(j)+y(j+1))/2.d0;	betaz(5)=z(k);
		betax(6)=x(i);					betay(6)=y(j);					betaz(6)=(z(k)+z(k+1))/2.d0;
		if (io(i,j,k)==1) then
			cappa2=cappap(x(i),y(j),z(k))
		else
			cappa2=cappan(x(i),y(j),z(k))
		endif
		
		if (irrpts(i,j,k)==0) then
			bt=beta(i,j,k)
			ijka(ijk+1)=ijka(ijk+1)+6
			ijka(ijkpt+1)	=ijk				-nx*ny
			ijka(ijkpt+2)	=ijk		-nx	
			ijka(ijkpt+3)	=ijk	-1
			ijka(ijkpt+4)	=ijk	+1							
			ijka(ijkpt+5)	=ijk		+nx
			ijka(ijkpt+6)	=ijk				+nx*ny
					
			sa(ijk)=-(1/dcel**2)*sum(bt)+cappa2
			sa(ijkpt+1:ijkpt+6)=(1/dcel**2)*bt(1:6)
			ijkpt=ijkpt+6
		else									! irregular points
			ixz=ijkpt+1
			if (io(i,j,k)==1) then
				do ll=1,6
					bt(ll)=betap(betax(ll),betay(ll),betaz(ll))
				enddo
				cappa2=cappap(x(i),y(j),z(k))
			elseif (io(i,j,k)==-1) then
				do ll=1,6
					bt(ll)=betan(betax(ll),betay(ll),betaz(ll))
				enddo
				cappa2=cappan(x(i),y(j),z(k))
			elseif (isrf(0,indirr(i,j,k))==1) then

				do ll=1,6
					bt(ll)=betap(betax(ll),betay(ll),betaz(ll))
				enddo
				cappa2=cappap(x(i),y(j),z(k))
			else 
				do ll=1,6
					bt(ll)=betan(betax(ll),betay(ll),betaz(ll))
				enddo
				cappa2=cappan(x(i),y(j),z(k))
			endif
			sa(ijk)=-(1/dcel**2)*sum(bt)+cappa2
			
			do l=1,6
				call idx6(i,j,k,l,idx)
				ix=idx(1);	iy=idx(2);	iz=idx(3);
				cof=(1/dcel**2)*bt(l)
				
				if (indirr(ix,iy,iz) .ne. 0 .and. io(i,j,k) .ne. 0) then	
					if (sum(abs(ftpts(indirr(ix,iy,iz),l,:))) >= 1.d-10) then	! Ficticious pt in this direction
						bftc(ijk)=bftc(ijk)-cof*ftc(indirr(ix,iy,iz),l)
						do m=1,nfc
							if (iftpts(indirr(ix,iy,iz),l,3*m-2)-i==0 .and. & 
							iftpts(indirr(ix,iy,iz),l,3*m-1)-j==0 &
							.and. iftpts(indirr(ix,iy,iz),l,3*m)-k==0) then
								sa(ijk)=sa(ijk)+cof*ftpts(indirr(ix,iy,iz),l,m)	! Uijk
							else
								if (abs(ftpts(indirr(ix,iy,iz),l,m))>=1.d-10) then
									if (abs(sa(ijkpt+1))<=1.d-10) then			! It has value already
										ijkpt=ijkpt+1
										sa(ijkpt)=cof*ftpts(indirr(ix,iy,iz),l,m)
										ijka(ijk+1)=ijka(ijk+1)+1
										ii=iftpts(indirr(ix,iy,iz),l,3*m-2)-i
										jj=iftpts(indirr(ix,iy,iz),l,3*m-1)-j
										kk=iftpts(indirr(ix,iy,iz),l,3*m)-k
										ijka(ijkpt)=ijk+kk*nx*ny+jj*nx+ii
									else										! First time it has value
										sa(ijkpt+1)=sa(ijkpt+1)+cof*ftpts(indirr(ix,iy,iz),l,m)
									endif
								endif
							endif
						enddo
					!elseif (iftwrg(indirr(ix,iy,iz),l)==1 .and. (iftpts(indirr(ix,iy,iz),l,1) .ne. -1))  then
					elseif (sum(abs(iftpts(indirr(ix,iy,iz),:,2)))>0 .and. iftwrg(indirr(ix,iy,iz),l)==1) then		!get fc pt from other drn
						itest=itest+1
						do ll=1,6
							if ((ll .ne. l) .and. sum(abs(ftpts(indirr(ix,iy,iz),ll,:)))>1.d-10) then
								bftc(ijk)=bftc(ijk)-cof*ftc(indirr(ix,iy,iz),ll)
								do m=1,nfc
									if (iftpts(indirr(ix,iy,iz),ll,3*m-2)-i==0 .and. & 
									iftpts(indirr(ix,iy,iz),ll,3*m-1)-j==0 &
									.and. iftpts(indirr(ix,iy,iz),ll,3*m)-k==0) then
										sa(ijk)=sa(ijk)+cof*ftpts(indirr(ix,iy,iz),ll,m)	! Uijk
									else
										if (abs(ftpts(indirr(ix,iy,iz),ll,m))>=1.d-10) then
											if (abs(sa(ijkpt+1))<=1.d-10) then			! It has value already
												ijkpt=ijkpt+1
												sa(ijkpt)=cof*ftpts(indirr(ix,iy,iz),ll,m)
												ijka(ijk+1)=ijka(ijk+1)+1
												ii=iftpts(indirr(ix,iy,iz),ll,3*m-2)-i
												jj=iftpts(indirr(ix,iy,iz),ll,3*m-1)-j
												kk=iftpts(indirr(ix,iy,iz),ll,3*m)-k
												ijka(ijkpt)=ijk+kk*nx*ny+jj*nx+ii
											else										! First time it has value
												sa(ijkpt+1)=sa(ijkpt+1)+cof*ftpts(indirr(ix,iy,iz),ll,m)
											endif
										endif
									endif
								enddo
							goto 500
							endif
						enddo
						print *,i,j,k,'missing'
						print *,ix,iy,iz
						500 continue	
!testing #######################################################################################
						!if (i==21 .and. j==41 .and. k==51 .and. l==5) then
							!print *,i,j,k,l
							!print *,ix,iy,iz
							!uintp=0.d0; utrue=0.d0
							!utrue=un(x(ix),y(iy),z(iz))
							!do m=1,16
							!	ii=iftpts(indirr(ix,iy,iz),l,3*m-2)
							!	jj=iftpts(indirr(ix,iy,iz),l,3*m-1)
							!	kk=iftpts(indirr(ix,iy,iz),l,3*m)
							!	print *,m,ftpts(indirr(ix,iy,iz),l,m),uexa(id3d(ii,jj,kk))
							!	uintp=uintp+ftpts(indirr(ix,iy,iz),l,m)*uexa(id3d(ii,jj,kk))
							!enddo
							!print *,ftc(indirr(ix,iy,iz),l),utrue
							!uintp=uintp+ftc(indirr(ix,iy,iz),l)
							!print *,'intp error', abs(uintp-utrue)
						!endif 
!###########################################################################################
					elseif (iftwrg(indirr(ix,iy,iz),l)==1 .and. sum(iftpts(indirr(ix,iy,iz),:,1))<0) then
						if (iftpts(indirr(ix,iy,iz),l,1)==-1) then !Adding the patch information
							ll=l
						else
							do ii=1,6
								if (iftpts(indirr(ix,iy,iz),ii,1)==-1) then
									ll=ii
									goto 15
								endif
							enddo

						15 continue	
						endif
						!print *,i,j,k,'llf'
						!print *,ix,iy,iz,'llf'
						do iii=1,iph
							jjj=idpatch(iii)
							call line2cube(ipp,jpp,kpp,jjj)
							if (ix==ipp .and. iy==jpp .and. iz==kpp) then
								ipch=iii
								goto 20
							endif
						enddo
						20 continue
						bftc(ijk)=bftc(ijk)-cof*pchftc(ipch,ll)

						!print *,ipch,ll,-cof,pchftc(ipch,ll)
						do m=1,50
							call line2cube(i0,j0,k0,icpatch(ipch,ll,m))
							if (i0-i==0 .and. j0-j==0 .and. k0-k==0) then
								sa(ijk)=sa(ijk)+cof*akpatch(ipch,ll,m)	! Uijk
							else
								if (abs(akpatch(ipch,ll,m))>=1.d-10 .and. (icpatch(ipch,ll,m) .ne. 0)) then
									if (abs(sa(ijkpt+1))<=1.d-10) then			! It has value already
										ijkpt=ijkpt+1
										sa(ijkpt)=cof*akpatch(ipch,ll,m)
										ijka(ijk+1)=ijka(ijk+1)+1
										ii=i0-i
										jj=j0-j
										kk=k0-k
										ijka(ijkpt)=ijk+kk*nx*ny+jj*nx+ii
									else										! First time it has value
										sa(ijkpt+1)=sa(ijkpt+1)+cof*akpatch(ipch,ll,m)
									endif
								endif
							endif
						enddo
					else						!Regular in this direction
						ijka(ijk+1)=ijka(ijk+1)+1
						ijkpt=ijkpt+1
						ii=ix-i;	jj=iy-j;	kk=iz-k
						ijka(ijkpt)=ijk+kk*nx*ny+jj*nx+ii
						sa(ijkpt)=cof
					endif
				elseif (isrf(l,indirr(i,j,k)) .ne. 0) then	!Adding the surface information

					isfc=isrf(l,indirr(i,j,k))
					bftc(ijk)=bftc(ijk)-cof*srfftc(isfc)
					! print *,i,j,k,ijk,bftc(ijk)
					do m=1,10
						if (isrfpts(isfc,3*m-2)-i==0 .and. isrfpts(isfc,3*m-1)-j==0 &
						.and. isrfpts(isfc,3*m)-k==0) then
							sa(ijk)=sa(ijk)+cof*srfpts(isfc,m)			! Uijk
						else
							if (abs(srfpts(isfc,m))>=1.d-10) then
								if (abs(sa(ijkpt+1))<=1.d-10) then		! It has value already
									ijkpt=ijkpt+1
									sa(ijkpt)=cof*srfpts(isfc,m)
									ijka(ijk+1)=ijka(ijk+1)+1
									ii=isrfpts(isfc,3*m-2)-i
									jj=isrfpts(isfc,3*m-1)-j
									kk=isrfpts(isfc,3*m)-k
									ijka(ijkpt)=ijk+kk*nx*ny+jj*nx+ii
								else									! First time it has value
									sa(ijkpt+1)=sa(ijkpt+1)+cof*srfpts(isfc,m)
								endif
							endif
						endif
					enddo	
				else													! No ficticious pt in this direction
					!if (io(i,j,k)==0) print *,i,j,k,l,'lsf_reg'
					ijka(ijk+1)=ijka(ijk+1)+1
					ijkpt=ijkpt+1
					ii=ix-i;	jj=iy-j;	kk=iz-k
					ijka(ijkpt)=ijk+kk*nx*ny+jj*nx+ii
					sa(ijkpt)=cof
				endif
			enddo
		jxz=ijkpt
		nrowi=ijka(ijk+1)-ijka(ijk)

		call sasort(nrowi,ijka(ixz:jxz),sa(ixz:jxz),nrdc)
		
		ijkpt=jxz-(nrowi-nrdc)
		ijka(ijk+1)=ijka(ijk)+nrdc
		!########################################## test
		!if (i==42 .and. j==16 .and. k==20) then
		!	do iii=ixz,ixz+nrdc-1
		!		print *,iii,sa(iii),ijka(iii)
		!		fleft=fleft+sa(iii)*uexa(ijka(iii))
		!	enddo
		!	fleft=fleft+sa(ijk)*uexa(ijk)
		!	fleft=fleft
		!	fright=bftc(ijk)
		!	print *,'pde error',i,j,k,fright-fleft
		!endif
		!##########################################
		endif
	endif
	!if (i==11 .and. j==9 .and. k==9) then
	!	print *,irrpts(i,j,k),io(i,j,k),bftc(ijk)
	!endif
enddo
enddo
enddo


ijka(nx*ny*nz+1)=ijkpt+1

nsize=ijka(nx*ny*nz+1)
!print *,'Solving Matrix Size=',nsize 

End

!------------------------------------------------------------------------------
subroutine fdsolver
use comdata
use pbedata
use bicg
implicit double precision(a-h,o-z)
real*8 d(0:2,0:2),xx(3),yy(3),zz(3)

sa=0.d0; 

d=0.d0
ww=x(2);	xx=x(1:3);
call weights(ww,xx,2,2,2,d)
xx=d(:,2)

d=0.d0
ww=y(2);	yy=y(1:3);
call weights(ww,yy,2,2,2,d)
yy=d(:,2)

d=0.d0
ww=z(2);	zz=z(1:3);
call weights(ww,zz,2,2,2,d)
zz=d(:,2)

sa=0.d0

ijka(1)=nx*ny*nz+2				!matrix size plus 2
ijkpt=nx*ny*nz+1				!ijkpt trace the how many elements filled in the two biconjugate vectors

do k=1,nz
	do j=1,ny
		do i=1,nx
			ijk=(k-1)*nx*ny+(j-1)*nx+i
			ijka(ijk+1)=ijka(ijk) !cumulated non-zero elements so far
			if ((i==1) .or. (i==nx) .or. (j==1) .or. (j==ny) .or. (k==1) .or. (k==nz)) then
				sa(ijk)=1.d0
				if (ibd==2 .or. ibd==3) then
					b(ijk)=bdcond(x(i),y(j),z(k))	! update boundary condition for eps=1
				else
					b(ijk)=up(x(i),y(j),z(k))                 
				endif
			else
				ijka(ijk+1)=ijka(ijk+1)+6
				ijka(ijkpt+1)	=ijk				-nx*ny
				ijka(ijkpt+2)	=ijk		-nx	
				ijka(ijkpt+3)	=ijk	-1
				ijka(ijkpt+4)	=ijk	+1							
				ijka(ijkpt+5)	=ijk		+nx
				ijka(ijkpt+6)	=ijk				+nx*ny
							
				sa(ijk)=yy(2)+xx(2)+zz(2)
				sa(ijkpt+1:ijkpt+6)=(/zz(1),yy(1),xx(1),xx(3),yy(3),zz(3)/)
				ijkpt=ijkpt+6
			endif
		enddo
	enddo
enddo
	ijka(nx*ny*nz+1)=ijkpt+1
End


!-------------------------------------------------------------------------------------------
SUBROUTINE indexx(n,narr,indx)
implicit double precision(a-h,o-z)
INTEGER n,indx(n),M,NSTACK
integer narr(n),na
PARAMETER (M=7,NSTACK=50)
INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
do j=1,n
	indx(j)=j
enddo
jstack=0
l=1
ir=n

1 if(ir-l.lt.M)then
	do j=l+1,ir
		indxt=indx(j)
		na=narr(indxt)
		do i=j-1,l,-1
			if(narr(indx(i)) .le. na)goto 2
			indx(i+1)=indx(i)
		enddo
		i=l-1
		2 indx(i+1)=indxt
	enddo

		if(jstack.eq.0)return
		ir=istack(jstack)
		l=istack(jstack-1)
		jstack=jstack-2
else
	k=(l+ir)/2
	itemp=indx(k)
	indx(k)=indx(l+1)
	indx(l+1)=itemp
	if(narr(indx(l)) .gt. narr(indx(ir)))then
		itemp=indx(l)
		indx(l)=indx(ir)
		indx(ir)=itemp
	endif

	if(narr(indx(l+1)) .gt. narr(indx(ir)))then
		itemp=indx(l+1)
		indx(l+1)=indx(ir)
		indx(ir)=itemp	
	endif

	if(narr(indx(l)) .gt. narr(indx(l+1)))then
		itemp=indx(l)
		indx(l)=indx(l+1)
		indx(l+1)=itemp
	endif

	i=l+1
	j=ir
	indxt=indx(l+1)
	na=narr(indxt)
	
	3 continue
	i=i+1
	if(narr(indx(i)) .lt. na)goto 3
		
	4 continue
	j=j-1
	if(narr(indx(j)) .gt. na)goto 4
	if(j.lt.i)goto 5
	itemp=indx(i)
	indx(i)=indx(j)
	indx(j)=itemp
	goto 3
		
	5 indx(l+1)=indx(j)
	indx(j)=indxt
	jstack=jstack+2
	if(jstack.gt.NSTACK) write(*,*) 'NSTACK too small in indexx'
	if(ir-i+1.ge.j-l)then
		istack(jstack)=ir
		istack(jstack-1)=i
		ir=j-1
	else
		istack(jstack)=j-1
		istack(jstack-1)=l
		l=i
	endif
endif
	goto 1
END


SUBROUTINE sort3(n,nra,ra)
implicit double precision(a-h,o-z)
INTEGER n,iwksp(n), nra(n), nwksp(n),indx(n)
INTEGER j
real*8 ra(n),wksp(n)
iwksp=0

call indexx(n,nra,iwksp)
do j=1,n 
	wksp(j)=ra(j)
enddo 

do j=1,n 
	ra(j)=wksp(iwksp(j))
enddo 


do j=1,n
	nwksp(j)=nra(j)
enddo

do j=1,n
	nra(j)=nwksp(iwksp(j))
enddo
End


SUBROUTINE rank(n,indx,irank)
implicit double precision(a-h,o-z)
INTEGER n,indx(n),irank(n)
!Given indx(1:n) as output from the routine indexx, this routine returns an array irank(1:n),
!the corresponding table of ranks.
INTEGER j
do j=1,n
	irank(indx(j))=j
enddo
return
END






!-------------------------------------------------------------------------------------------
subroutine pbeqmatrix
use comdata
use pbedata
use bicg
implicit double precision(a-h,o-z)
real*8 d(0:2,0:2),xx(3),yy(3),zz(3),fdcof(6),beta_apbs(6)
integer idx(3)

sa=0.d0; 

d=0.d0
ww=x(2);	xx=x(1:3);
call weights(ww,xx,2,2,2,d)
xx=d(:,2)

d=0.d0
ww=y(2);	yy=y(1:3);
call weights(ww,yy,2,2,2,d)
yy=d(:,2)

d=0.d0
ww=z(2);	zz=z(1:3);
call weights(ww,zz,2,2,2,d)
zz=d(:,2)

sa=0.d0

ijka(1)=nx*ny*nz+2				!matrix size plus 2
ijkpt=nx*ny*nz+1				!ijkpt trace the how many elements filled in the two biconjugate vectors

do k=1,nz
	do j=1,ny
		do i=1,nx
			ijk=(k-1)*nx*ny+(j-1)*nx+i
			ijka(ijk+1)=ijka(ijk) !cumulated non-zero elements so far
			if ((i==1) .or. (i==nx) .or. (j==1) .or. (j==ny) .or. (k==1) .or. (k==nz)) then
				sa(ijk)=1.d0
			else
				ijka(ijk+1)=ijka(ijk+1)+6
				ijka(ijkpt+1)	=ijk				-nx*ny
				ijka(ijkpt+2)	=ijk		-nx	
				ijka(ijkpt+3)	=ijk	-1
				ijka(ijkpt+4)	=ijk	+1							
				ijka(ijkpt+5)	=ijk		+nx
				ijka(ijkpt+6)	=ijk				+nx*ny
							
				if (irrpts(i,j,k)==0) then					! For regular points, just standard way
					if (io(i,j,k)<0) then
						beta_apbs=betan(x(i),y(j),z(k))		! Careful, for pbe, beta has only two values
					else
						beta_apbs=betap(x(i),y(j),z(k))
					endif
				else										! Irregular points, different surfaces varies
					do ik=1,6
						call idx6(i,j,k,ik,idx)
						xbt=(x(i)+x(idx(1)))/2.d0
						ybt=(y(j)+y(idx(2)))/2.d0
						zbt=(z(k)+z(idx(3)))/2.d0
																	! Irregular direction on irregular pt
						if (isf==0) then						! Analytic surface
							if (varphi(xbt,ybt,zbt)<1.d-10) then
								beta_apbs(ik)=betan(xbt,ybt,zbt)
							else
								beta_apbs(ik)=betap(xbt,ybt,zbt)
							endif
						else									! Msms surface and so on
							if (io(i,j,k)+io(idx(1),idx(2),idx(3)) .ne. 0) then 
								! Regular direction on irregular point
								! regular point		(+1)+(+1)=2 or (-1)+(-1)=-2
								! interface point	(+1)+0=1 or (-1)+0=-1
								if (io(idx(1),idx(2),idx(3))<=0) then
									beta_apbs(ik)=betan(xbt,ybt,zbt)
								else
									beta_apbs(ik)=betap(xbt,ybt,zbt)
								endif
							else		
								if (io(i,j,k)>0) then
									ikk=7-ik;	iin=idx(1);	jin=idx(2);	kin=idx(3)
								else
									ikk=ik;		iin=i;		jin=j;		kin=k
								endif

								if (i .ne. idx(1)) then ! Find inside point and beta location
									ptin=x(iin);		ptbt=xbt
								elseif (j .ne. idx(2)) then
									ptin=y(jin);		ptbt=ybt
								else
									ptin=z(kin);		ptbt=zbt
								endif

								irr=indirr(iin,jin,kin)
								ptsf=rrintf(ikk,irr)	! Find surface point

								if (abs(ptbt-ptin)<abs(ptin-ptsf)) then   
									beta_apbs(ik)=betan(xbt,ybt,zbt) ! ptbt is inbetween ptsf and ptin
								else
									beta_apbs(ik)=betap(xbt,ybt,zbt)
								endif
							endif
						endif		
					enddo
				endif
				
				fdcof=(/zz(1),yy(1),xx(1),xx(3),yy(3),zz(3)/)
				do ii=1,6
					sa(ijkpt+ii)=fdcof(ii)*beta_apbs(ii)
					sa(ijk)=sa(ijk)-sa(ijkpt+ii)
				enddo
				ijkpt=ijkpt+6
			endif
		enddo
	enddo
enddo
	ijka(nx*ny*nz+1)=ijkpt+1
End
