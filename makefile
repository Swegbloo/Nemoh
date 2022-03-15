# makefile written by Christophe Peyrard from EDF R&D
# edited by Adrien Combourieu from INNOSEA (adrien.combourieu@innosea.fr)
# edited by Ruddy Kurnia, ECN 2022
#

#COMPILATEUR
gtest=$(shell which gfortran 2> /dev/null | grep -o gfortran)
itest=$(shell which ifort 2> /dev/null | grep -o ifort)

ifeq ($(gtest), gfortran)
	FC=gfortran
	FFLAGS=-cpp -DGNUFORT -O2 -ffree-line-length-none -c -fdefault-real-8
# if with  -fdefault-real-8 change the source library for solver zPackgmres.f
# 	then change in Solver/Core/Com_Var.f90 ID_DP=1
# else cPackgmres.f
# 	then change in Solver/Core/Com_Var.f90 ID_DP=0
endif

ifeq ($(itest), ifort)
	FC=ifort
	FFLAGS=-O3 -c -cpp
# if with -r8 change in Solver/Core/Com_Var.f90 ID_DP=1
# else  change in Solver/Core/Com_Var.f90 ID_DP=0
endif

outputdir=./bin

#SOURCES FORTRAN Mesh(modules de maillage)
SRCM=./Common/Identification.f90\
./Mesh/calCol.f90\
./Mesh/coque.f90\
./Mesh/ExMaillage.f90\
./Mesh/hydre.f90\
./Mesh/Mailleur.f90\
./Mesh/mesh.f90\

# LISTE DES .o mesh
#TRANSFORME f90 en o  
OBJM=$(SRCM:.f90=.o)

#SOURCES FORTRAN Hydrostatic
SRCHS=./Common/Environment.f90\
./Common/Identification.f90\
./preProcessor/Mesh.f90\
./Mesh/hydre.f90\
./Mesh/coque.f90\
./Mesh/hydrostatic_cal.f90\

# LISTE DES .o hydrostatic
#TRANSFORME f90 en o  
OBJHS=$(SRCHS:.f90=.o)


#SOURCES FORTRAN preProc(modules de preprocessing)
SRCP=./Common/Identification.f90\
./Common/Environment.f90\
./preProcessor/Mesh.f90\
./preProcessor/BodyConditions.f90\
./preProcessor/Integration.f90\
./preProcessor/Main.f90\

# LISTE DES .o preProc
#TRANSFORME f90 en o  
OBJP=$(SRCP:.f90=.o)

#SOURCES FORTRAN Solver(modules de preprocessing)
SRCS=./Solver/Core/FIC_COM.f90\
./Solver/Core/COM_VAR.f90\
./Common/Environment.f90\
./Common/Identification.f90\
./Common/Mesh.f90\
./Solver/Core/Bodyconditions.f90\
./Solver/Core/ELEMENTARY_FNS.f90\
./Solver/Core/PREPARE_MESH.f90\
./Solver/Core/INITIALIZATION.f90\
./Solver/Core/OUTPUT.f90\
./Solver/Core/cPackgmres.f\
./Solver/Core/zPackgmres.f\
./Solver/Core/blas_rot.f\
./Solver/Core/M_SOLVER.f90\
./Solver/Core/ALLOCATE_DATA.f90\
./Solver/Core/COMPUTE_GREEN_INFD.f90\
./Solver/Core/SOLVE_BEM_INFD_DIRECT.f90\
./Solver/Core/COMPUTE_GREEN_FD.f90\
./Solver/Core/SOLVE_BEM_FD_DIRECT.f90\
./Solver/Core/SOLVE_BEM_FD_GMRES.f90\
./Solver/Core/SOLVE_BEM_INFD_GMRES.f90\
./Solver/Core/SOLVE_BEM.f90\
./Solver/Core/COMPUTE_KOCHIN.f90\
./Solver/Core/COMPUTE_GREEN_FREESURFACE.f90\
./Solver/Core/COMPUTE_POTENTIAL_DOMAIN.f90\
./Solver/NEMOH.f90\
./Solver/Core/DEALLOCATE_DATA.f90\



# LISTE DES .o solver
#TRANSFORME f90 en o  
OBJS=$(SRCS:.f90=.o)

#SOURCES FORTRAN postProc(modules de postprocessing)
SRCO=./Common/Identification.f90\
./Common/Environment.f90\
./Common/Results.f90\
./Common/Mesh.f90\
./postProcessor/Compute_RAOs.f90\
./postProcessor/IRF.f90\
./postProcessor/Plot_WaveElevation.f90\
./postProcessor/Main.f90\

# LISTE DES .o postProc
#TRANSFORME f90 en o  
OBJO=$(SRCO:.f90=.o)


build: clean bin msh hyd pre solver post  

bin:
	mkdir -p $(outputdir)

bin:
	mkdir -p $(outputdir)
#
#Build Mesh executable
msh:	meshexe
#Rules to Build MAIN EXECUTABLE  (dependances et regle d'execution)
meshexe:	$(OBJM) 
		$(FC)  -o meshExe $(OBJM)
#

#Build Hydrostatic cal executable
hyd:	hydstexe
#Rules to Build MAIN EXECUTABLE  (dependances et regle d'execution)
hydstexe:	$(OBJHS) 
		$(FC) -o hydroStaticExe $(OBJHS)
#

#Build preProc executable
pre:	preProc
#Rules to Build MAIN EXECUTABLE  (dependances et regle d'execution)
preProc:	$(OBJP) 
		$(FC) -o preProcExe $(OBJP)

#
#Rules to Build MAIN EXECUTABLE  (dependances et regle d'execution)
solver:	        $(OBJS)
		$(FC) -llapack -lblas -o solverExe $(OBJS)

#
#Build postProc executable
post:	postProc
#Rules to Build MAIN EXECUTABLE  (dependances et regle d'execution)
postProc:	$(OBJO) 
		$(FC) -o postProcExe $(OBJO)

		
# Rules for .f compilation
.f.o:
	$(FC) $(FFLAGS) $<
%.o:	%.f90
	$(FC) $(FFLAGS) $< -o $@
	

# # #Copy to local bin directory
install: build
	mv meshExe hydroStaticExe preProcExe solverExe postProcExe ./bin/

# Remove *.o and *.mod
clean:
	rm -f *.o */*.o */*/*.o *.mod
	
#remove executable as well
cleanall:
	rm -r ./bin
	rm ./preProcessor/*.o
	rm ./postProcessor/*.o
	rm ./Solver/*.o
	rm ./Solver/Core/*.o
	rm ./Mesh/*.o 
	rm ./Common/*.o
	rm -f *.o *.mod
