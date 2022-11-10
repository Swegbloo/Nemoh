# Makefile written by Christophe Peyrard from EDF R&D
# Extended to OS X by Yi-Hsiang Hu & Eliot Quin from NREL
# Clean up by Matthieu Ancellin (only tested with gfortran on Linux for the moment)

# Compiler
gtest=$(shell which gfortran 2> /dev/null | grep -o gfortran)
itest=$(shell which ifort 2> /dev/null | grep -o ifort)
# M.A.: this test does not work by me...

MOD_DIR=/tmp/

ifeq ($(gtest), gfortran)
	FC=gfortran
	FFLAGS=  -c                                     # No linker (yet)
	FFLAGS+= -g                                     # Add extra informations for debugging
	FFLAGS+= -O2                                    # Optimization level
	FFLAGS+= -J$(MOD_DIR)                           # Where to put .mod files
	FFLAGS+= -cpp -DGNUFORT -ffree-line-length-none # Run preprocessor
	#FFLAGS+=-fdefault-real-8			# forcing variables to be double precision
	#if with -fdefault-real-8  change in Common/Constants.f90 ID_DP=1 else ID_DP=0
endif

ifeq ($(itest), ifort)
	FC=ifort
	FFLAGS=-c -cpp
	FFLAGS+= -O2                                    # Optimization level
	FFLAGS+= -g                                     # Add extra informations for debugging
	FFLAGS+= -module $(MOD_DIR)                           # Where to put .mod files
	# FFLAGS+= -r8					# forcing variables to be double precision
	# if with -r8  change in Common/Constants.f90 ID_DP=1 else ID_DP=0
endif

# Output directory
outputdir=./bin

# Default rule: build all
.PHONY: all clean_all remake
all:		mesh preProc solver postProc hydroCal

clean_all:	clean_mesh clean_preProc clean_solver clean_postProc clean_hydroCal
			@rm -f $(MOD_DIR)/*.mod

remake:		clean_all all

# Rule to compile f90 file
%.o:	%.f90
		@$(FC) $(FFLAGS) $< -o $@

# Rule to compile f file
%.o:	%.f
		@$(FC) $(FFLAGS) $< -o $@

##################
#  Meshing tool  #
##################

# Sources (relative to DIRM)
SRCM=./Common/Identification.f90\
./Mesh/calCol.f90\
./Mesh/coque.f90\
./Mesh/ExMaillage.f90\
./Mesh/hydre.f90\
./Mesh/Mailleur.f90\
./Mesh/mesh.f90

OBJM=$(SRCM:.f90=.o)

# Rules to build
mesh:		$(OBJM)
			@test -d $(outputdir) || mkdir $(outputdir)
			@$(FC) -o $(outputdir)/mesh $(OBJM)
			@echo "Meshing tool compilation successful!"

clean_mesh:
			@rm -f $(OBJM)
			@rm -f $(outputdir)/mesh
##################
#  HYDROSTATIC  #
##################

# Sources (relative to DIRM)
SRCHS=./Common/Identification.f90\
./Common/Elementary_functions.f90\
./Common/Environment.f90\
./Common/Mesh.f90\
./Mesh/coque.f90\
./Mesh/hydre.f90\
./Mesh/hydrostatic_cal.f90

OBJHS=$(SRCHS:.f90=.o)

# Rules to build
hydroCal:	$(OBJHS)
			@test -d $(outputdir) || mkdir $(outputdir)
			@$(FC) -o $(outputdir)/hydrosCal $(OBJHS)
			@echo "hydrosCal compilation successful!"

clean_hydroCal:
			@rm -f $(OBJHS)
			@rm -f $(outputdir)/hydrosCal
###################
#  Pre-processor  #
###################

# Sources
SRCP=./Common/Constants.f90\
./Common/Elementary_functions.f90\
./Common/Identification.f90\
./Common/Environment.f90\
./Common/MNemohCal.f90\
./preProcessor/Mesh.f90\
./preProcessor/BodyConditions.f90\
./preProcessor/Integration.f90\
./preProcessor/Main.f90

OBJP=$(SRCP:.f90=.o)

# Rules to build
preProc:	$(OBJP)
			@test -d $(outputdir) || mkdir $(outputdir)
			@$(FC) -o $(outputdir)/preProc $(OBJP)
			@echo "Preprocessor compilation successful!"

clean_preProc:
			@rm -f $(OBJP)
			@rm -f $(outputdir)/preProc

############
#  Solver  #
############
SRCSEXT=./Solver/Core/cPackgmres.f\
./Solver/Core/zPackgmres.f\
./Solver/Core/blas_rot.f

OBJSEXT=$(SRCSEXT:.f=.o)

SRCS=./Common/Constants.f90\
./Common/Logfile.f90\
./Common/Elementary_functions.f90\
./Common/Bodyconditions.f90\
./Common/Environment.f90\
./Common/Mesh.f90\
./Common/Face.f90\
./Solver/Core/OUTPUT.f90\
./Solver/Core/M_SOLVER.f90\
./Solver/Core/GREEN_1.f90\
./Solver/Core/INITIALIZE_GREEN.f90\
./Solver/Core/GREEN_2.f90\
./Solver/Core/SOLVE_BEM_DIRECT.f90\
./Solver/Core/KOCHIN.f90\
./Solver/Core/FREESURFACE.f90\
./Solver/Core/FORCES.f90\
./Solver/NEMOH.f90

OBJS=$(SRCS:.f90=.o)

# Rules to build
solver:		$(OBJSEXT) $(OBJS)
			@test -d $(outputdir) || mkdir $(outputdir)
			@$(FC) -llapack -lblas -o $(outputdir)/solver $(OBJSEXT) $(OBJS)
			@echo "Solver compilation succesful!"

clean_solver:
			@rm -f $(OBJSEXT) $(OBJS)
			@rm -f $(outputdir)/solver

####################
#  Post-processor  #
####################
SRCOEXT=./Solver/Core/cPackgmres.f\
./Solver/Core/zPackgmres.f\
./Solver/Core/blas_rot.f

OBJOEXT=$(SRCOEXT:.f=.o)

# Sources
SRCO=./Common/Constants.f90\
./Common/Elementary_functions.f90\
./Common/Identification.f90\
./Common/Environment.f90\
./Common/Results.f90\
./Common/Mesh.f90\
./Common/MNemohCal.f90\
./Solver/Core/M_SOLVER.f90\
./postProcessor/MPP_ReadInputFiles.f90\
./postProcessor/MPP_Compute_RAOs.f90\
./postProcessor/IRF.f90\
./postProcessor/Plot_WaveElevation.f90\
./postProcessor/Main.f90

OBJO=$(SRCO:.f90=.o)

# Rules to build
postProc:	$(OBJOEXT) $(OBJO)
			@test -d $(outputdir) || mkdir $(outputdir)
			@$(FC) -llapack -lblas -o $(outputdir)/postProc $(OBJOEXT) $(OBJO)
			@echo "Postprocessor compilation succesful!"

clean_postProc:
			@rm -f $(OBJO) $(OBJOEXT)
			@rm -f $(outputdir)/postProc

################
#  Test cases  #
################
# Verification directory
verdir=./
.PHONY: run_cylinder clean_cylinder
run_cylinder: preProc solver postProc
	$(MAKE) -C $(verdir)/Verification/Cylinder/ run

clean_cylinder:
	$(MAKE) -C $(verdir)/Verification/Cylinder/ clean

.PHONY: run_nonsymmetrical clean_nonsymmetrical
run_nonsymmetrical: preProc solver postProc
	$(MAKE) -C $(verdir)/Verification/NonSymmetrical/ run

clean_nonsymmetrical:
	$(MAKE) -C $(verdir)/Verification/NonSymmetrical/ clean

.PHONY: run_2Bodies clean_2Bodies
run_2Bodies: preProc solver postProc
	$(MAKE) -C $(verdir)/Verification/2Bodies/ run

clean_2Bodies:
	$(MAKE) -C $(verdir)/Verification/2Bodies/ clean

.PHONY: run_Postprocessing clean_Postprocessing
run_Postprocessing: preProc solver postProc
	$(MAKE) -C $(verdir)/Verification/Postprocessing/ run

clean_Postprocessing:
	$(MAKE) -C $(verdir)/Verification/Postprocessing/ clean


.PHONY: run_quicktest clean_quicktest
run_quicktest: preProc solver postProc
	@echo ""
	@echo "Sphere"
	@$(MAKE) --silent -C $(verdir)/Verification/QuickTests/1_Sphere/                     test
	@echo ""
	@echo "Sphere using y-symmetry"
	@$(MAKE) --silent -C $(verdir)/Verification/QuickTests/2_SymmetricSphere/            test
	@echo ""
	@echo "Sphere in finite depth"
	@$(MAKE) --silent -C $(verdir)/Verification/QuickTests/3_FiniteDepthSphere/          test
	@echo ""
	@echo "Sphere in finite depth using y-symmetry"
	@$(MAKE) --silent -C $(verdir)/Verification/QuickTests/4_SymmetricFiniteDepthSphere/ test
	@echo ""
	@echo "Alien sphere"
	@$(MAKE) --silent -C $(verdir)/Verification/QuickTests/5_AlienSphere/                test

clean_quicktest:
	@$(MAKE) --silent -C $(verdir)/Verification/QuickTests/1_Sphere/                     clean
	@$(MAKE) --silent -C $(verdir)/Verification/QuickTests/2_SymmetricSphere/            clean
	@$(MAKE) --silent -C $(verdir)/Verification/QuickTests/3_FiniteDepthSphere/          clean
	@$(MAKE) --silent -C $(verdir)/Verification/QuickTests/4_SymmetricFiniteDepthSphere/ clean
	@$(MAKE) --silent -C $(verdir)/Verification/QuickTests/5_AlienSphere/                clean
