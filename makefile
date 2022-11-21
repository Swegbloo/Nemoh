# Makefile written by Christophe Peyrard from EDF R&D
# Extended to OS X by Yi-Hsiang Hu & Eliot Quin from NREL
# Clean up by Matthieu Ancellin (only tested with gfortran on Linux for the moment)

# Compiler
gtest=$(shell which gfortran 2> /dev/null | grep -o gfortran)
itest=$(shell which ifort 2> /dev/null | grep -o ifort)
# M.A.: this test does not work by me...

MOD_DIR=/tmp/

# External libraries location
LIBDIR=/home/libre-service/Documents/Codes/util/

ifeq ($(gtest), gfortran)
	FC=gfortran
	FFLAGS=  -c                                     # No linker (yet)
	FFLAGS+= -g                                     # Add extra informations for debugging
	FFLAGS+= -O2                                    # Optimization level
	FFLAGS+= -J$(MOD_DIR)                           # Where to put .mod files
	FFLAGS+= -cpp -DGNUFORT -ffree-line-length-none # Run preprocessor
	#FFLAGS+=-fdefault-real-8			# forcing variables to be double precision
	#if with -fdefault-real-8  change in Common/Constants.f90 ID_DP=1 else ID_DP=0
	FFLAGS+= -L$(LIBDIR)
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
			@$(FC) -o $(outputdir)/solver $(OBJSEXT) $(OBJS) -llapack -lblas
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
			@$(FC) -o $(outputdir)/postProc $(OBJOEXT) $(OBJO) -llapack -lblas
			@echo "Postprocessor compilation succesful!"

clean_postProc:
			@rm -f $(OBJO) $(OBJOEXT)
			@rm -f $(outputdir)/postProc

################
#  Test cases  #
################
# Verification directory
verdir=.
.PHONY: run_1_cylinder clean_1_cylinder
run_1_cylinder:
	$(MAKE) -C $(verdir)/Verification/1_Cylinder/ run

clean_1_cylinder:
	$(MAKE) -C $(verdir)/Verification/1_Cylinder/ clean

.PHONY: run_2_2Bodies clean_2_2Bodies
run_2_2Bodies:
	$(MAKE) -C $(verdir)/Verification/2_2Bodies/ run

clean_2_2Bodies:
	$(MAKE) -C $(verdir)/Verification/2_2Bodies/ clean

.PHONY: run_3_nonsymmetrical clean_3_nonsymmetrical
run_3_nonsymmetrical:
	$(MAKE) -C $(verdir)/Verification/3_NonSymmetrical/ run

clean_3_nonsymmetrical:
	$(MAKE) -C $(verdir)/Verification/3_NonSymmetrical/ clean

.PHONY: run_4_Postprocessing clean_4_Postprocessing
run_4_Postprocessing:
	$(MAKE) -C $(verdir)/Verification/4_Postprocessing/ run

clean_4_Postprocessing:
	$(MAKE) -C $(verdir)/Verification/4_Postprocessing/ clean


.PHONY: run_5_quicktest clean_5_quicktest
run_5_quicktest:
	@echo ""
	@echo "Sphere"
	@$(MAKE) --silent -C $(verdir)/Verification/5_QuickTests/1_Sphere/                     test
	@echo ""
	@echo "Sphere using y-symmetry"
	@$(MAKE) --silent -C $(verdir)/Verification/5_QuickTests/2_SymmetricSphere/            test
	@echo ""
	@echo "Sphere in finite depth"
	@$(MAKE) --silent -C $(verdir)/Verification/5_QuickTests/3_FiniteDepthSphere/          test
	@echo ""
	@echo "Sphere in finite depth using y-symmetry"
	@$(MAKE) --silent -C $(verdir)/Verification/5_QuickTests/4_SymmetricFiniteDepthSphere/ test
	@echo ""
	@echo "Alien sphere"
	@$(MAKE) --silent -C $(verdir)/Verification/5_QuickTests/5_AlienSphere/                test

clean_5_quicktest:
	@$(MAKE) --silent -C $(verdir)/Verification/5_QuickTests/1_Sphere/                     clean
	@$(MAKE) --silent -C $(verdir)/Verification/5_QuickTests/2_SymmetricSphere/            clean
	@$(MAKE) --silent -C $(verdir)/Verification/5_QuickTests/3_FiniteDepthSphere/          clean
	@$(MAKE) --silent -C $(verdir)/Verification/5_QuickTests/4_SymmetricFiniteDepthSphere/ clean
	@$(MAKE) --silent -C $(verdir)/Verification/5_QuickTests/5_AlienSphere/                clean

.PHONY: run_6_box_coarsemesh clean_6_box_coarsemesh
run_6_box_coarsemesh:
	$(MAKE) -C $(verdir)/Verification/6_box_coarsemesh/ run

clean_6_box_coarsemesh:
	$(MAKE) -C $(verdir)/Verification/6_box_coarsemesh/ clean

clean_all_tests:	clean_1_cylinder clean_2_2Bodies clean_3_nonsymmetrical clean_4_Postprocessing clean_5_quicktest clean_6_box_coarsemesh
