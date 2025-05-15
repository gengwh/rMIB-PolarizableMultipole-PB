# compile necessary components of mibpb3
 F90 = gfortran
 flag = -c #-O2 
#F90 = ifort
#flag = -c -heap-arrays -fast
rMIB.exe: comfile.o addin_3d_weihua.o bicggt_new.o chgdist.o fdweights.o \
       kirkwood.o lu.o main.o mibpb_para.o mibmatrix.o msmsmibw13.o para_func_oneb.o \
       topology.o para_interface.o nonlinearsolver.o readin.o eses_readin.o 
	$(F90) -o rMIB.exe *.o 
comfile.o :		comfile.f90
	$(F90) $(flag) comfile.f90
addin_3d_weihua.o:	addin_3d_weihua.f90
	$(F90) $(flag) addin_3d_weihua.f90
bicggt_new.o:		bicggt_new.f90
	$(F90) $(flag) bicggt_new.f90
chgdist.o :		chgdist.f90
	$(F90) $(flag) chgdist.f90
fdweights.o :		fdweights.f90
	$(F90) $(flag) fdweights.f90
kirkwood.o :		kirkwood.f90 
	$(F90) $(flag) kirkwood.f90
lu.o :			lu.f90
	$(F90) $(flag) lu.f90
main.o :	main.f90
	$(F90) $(flag) main.f90
mibpb_para.o :	mibpb_para.f90
	$(F90) $(flag) mibpb_para.f90
mibmatrix.o :		mibmatrix.f90
	$(F90) $(flag) mibmatrix.f90
msmsmibw13.o :		msmsmibw13.f90
	$(F90) $(flag) msmsmibw13.f90
para_func_oneb.o :	para_func_oneb.f90
	$(F90) $(flag) para_func_oneb.f90
para_interface.o : para_interface.f90 
	$(F90) $(flag) para_interface.f90
topology.o :		topology.f90
	$(F90) $(flag) topology.f90
vdw_surface.o:		vdw_surface.f90
	$(F90) $(flag) vdw_surface.f90
nonlinearsolver.o: nonlinearsolver.f90
	$(F90) $(flag) nonlinearsolver.f90
readin.o: readin.f90
	$(F90) $(flag) readin.f90
eses_readin.o: eses_readin.f90
	$(F90) $(flag) eses_readin.f90
clean:
	rm *.o *.mod *.dx *.txt *.obj *.out *.exe surfaces/*.txt surfaces/*box 
