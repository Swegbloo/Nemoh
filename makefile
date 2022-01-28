# makefile written by Adrien Combourieu from INNOSEA (adrien.combourieu@innosea.fr)

#COMPILATEUR  
 FC=ifort
#FC=gfortran
#OPTIONS  
# FFLAGS= -c 
FFLAGS=  -O3 -ffree-line-length-0 -fbounds-check -c -r8
FFLAGS2=  -O3 -ffree-line-length-0 -finit-local-zero -mcmodel=medium -fbounds-check
#FFLAGS2= -O3 -fopenmp -ffree-line-length-0 -finit-local-zero -mcmodel=medium -fbounds-check -fno-automatic
outputdir=./bin

#### same process for second order (QTF) ################
#SOURCES FORTRAN QTFpreProc(modules de postprocessing)
#ELEMENTARY_FNS is added by RK for preparing CREK 2nd term Green function
SRCQPRE=./Common/Identification.f90\
./Common/Environment.f90\
./Common/Mesh.f90\
./QTF/QTFpreProcessor/QTFCOM_VAR.f90\
./QTF/QTFpreProcessor/QTFRead.f90\
./QTF/QTFpreProcessor/QTFWrite.f90\
./QTF/QTFpreProcessor/QTFGeom.f90\
./QTF/QTFpreProcessor/ELEMENTARY_FNS.f90\
./QTF/QTFpreProcessor/QTFBaseFunctions.f90\
./QTF/QTFpreProcessor/QTFInit.f90\
./QTF/QTFpreProcessor/Main.f90\

#### same process for second order (QTF) ################
OBJQPRE=$(SRCQPRE:.f90=.o)

################
build: bin qtfpre qtfsol

bin:
	mkdir -p $(outputdir)

#
#Build QTFpreProc executable
qtfpre:	QTFpreProc
qtfsol: duokap hasbo hasfs hasfscalc asymp QTFSolver

#Rules to Build MAIN EXECUTABLE of QTFpreProc
QTFpreProc:	$(OBJQPRE) 
		$(FC) -o QTFpreProc $(OBJQPRE)


#Rules to Build QTFSolver
duokap:	./QTF/QTFSolver/duokap.f ./QTF/QTFSolver/para.o
		$(FC) -c ./QTF/QTFSolver/para.f90
		$(FC) $(FFLAGS2) ./QTF/QTFSolver/duokap.f -o duokap
		
hasbo: ./QTF/QTFSolver/hasbo.f ./QTF/QTFSolver/para.o
		$(FC) -c ./QTF/QTFSolver/para.f90
		$(FC) $(FFLAGS2) ./QTF/QTFSolver/hasbo.f -o hasbo
	
hasfs: ./QTF/QTFSolver/hasfs.f ./QTF/QTFSolver/subhasp.f ./QTF/QTFSolver/para.o 
		$(FC) -c ./QTF/QTFSolver/para.f90
		$(FC) $(FFLAGS2) ./QTF/QTFSolver/hasfs22.f ./QTF/QTFSolver/subhasp.f -o hasfs

hasfscalc: ./QTF/QTFSolver/hasfscalc.f90 ./QTF/QTFSolver/para.o
		$(FC) -c ./QTF/QTFSolver/para.f90
		$(FC) $(FFLAGS2) ./QTF/QTFSolver/hasfscalc22.f90 -o hasfscalc

asymp: ./QTF/QTFSolver/hasfsasymp.f90 ./QTF/QTFSolver/para.o
		$(FC) -c ./QTF/QTFSolver/para.f90
		$(FC) $(FFLAGS2) ./QTF/QTFSolver/hasfsasymp22.f90 -o asymp

QTFSolver: ./QTF/QTFSolver/Solver2.f90 ./QTF/QTFSolver/para.o
		$(FC) $(FFLAGS2) ./QTF/QTFSolver/Solver2.f90 -o QTFSolver


# rules for object files
.f.o:
	$(FC) $(FFLAGS) $<
%.o::	%.f90
	$(FC) $(FFLAGS) $< -o $@
	
	
# # #Copy to local bin directory
install : build
	mv QTFpreProc duokap hasbo hasfs hasfscalc asymp QTFSolver $(outputdir)

# Remove *.o and main executable
 clean:	
	rm -f *.mod *.o

 cleanall:
	rm ./Common/*.o
	rm ./QTF/QTFpreProcessor/*.o
	rm ./QTF/QTFSolver/*.o
	rm $(outputdir)/duokap
	rm $(outputdir)/hasbo
	rm $(outputdir)/hasfs
	rm $(outputdir)/hasfscalc
	rm $(outputdir)/asymp
	rm $(outputdir)/QTFSolver
	rm $(outputdir)/QTFpreProc
	rm *.o
	rm *.mod
