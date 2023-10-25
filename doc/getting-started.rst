
###############
Getting started
###############

This section describes the installation procedure and how to run the codes. It also includes a description of input-output files.


****************
Executable files
****************

Executable/binary files are provided in the ``Release`` section of the `GitLab project <https://gitlab.com/lheea/Nemoh>`__. They can be used directly without the need to do the compilation procedure described in the next subsection. Windows and Linux executable files are provided. In the case of a compilation on your computer, executable files will be located in the ``bin/`` folder.
The following executable files are available:

-  NEMOH1: ``mesh``, ``preProc``, ``hydrosCal``, ``solver``, ``postProc``,

-  NEMOH2: ``QTFpreProc``, ``QTFsolver``, ``QTFpostProc``.

Note that a Matlab wrapper is provided to use those executable files in a Matlab environment. More details are provided in :ref:`getting-started:supporting matlab files`.

************
Installation
************

This procedure is intended for a developer to perform the compilation after changing the source code. Users may skip this step.

As a first preliminary step, it is necessary to install the following external libraries that are used by NEMOH:

-  BLAS, https://netlib.org/blas/

-  LAPACK, https://netlib.org/lapack/

Windows platform
================

An updated manual will provide the details of the compilation on a Windows platform.

Linux platform
==============

Before compiling, the following tools need to be available:

-  A Fortran compiler. The code has been tested using:

   -  gfortran, https://gcc.gnu.org/wiki/GFortran,

   -  intel fortran compiler, `ifort <https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html#gs.jik1s6>`__.

-  CMake, https://cmake.org/, a cross-platform tool for building and testing the software package.

Compile all Nemoh executables using CMake (from the root of the repository):

.. code:: bash

   cmake -S. -Bbuild
       cmake --build build

The resulting executables are in the ``bin/`` directory. To compile only one of the executables, use the ``–target`` option of CMake. The available targets are:

-  for NEMOH: ``mesh``, ``preProc``, ``hydrosCal``, ``solver``, ``postProc``

-  for NEMOH QTF:``QTFpreProc``, ``QTFsolver``, ``QTFpostProc``

The choice of the compiler is left to CMake, but can be overridden by setting the ``CMAKE_Fortran_COMPILER`` at the configuration step, e.g.:

.. code:: bash

   cmake -S. -Bbuild -DCMAKE_Fortran_COMPILER=gfortran

After building, the tests can be run from the ``build/`` directory:

.. code:: bash

   ctest -V -j <N_concurrent>

Where ``<N_ concurrent>`` is the number of simultaneous workers (processes). The tests can be restricted using their labels and the ``-L`` option of ctest:

.. code:: bash

   ctest -V -j <N_concurrent> -L <label>

Where label is one of the following:

-  ``NEMOH1``: only the non-QTF test cases

-  ``PREPROC``: only the pre-processing operations

-  ``SOLVER``: only the solving operations (depend on the pre-processing tests)

-  ``POSTPROC``: only the post-processing operations (depend on the pre-processing and solving tests)

-  ``NEMOH2``: only the QTF test cases

-  ``QTF``: only the computation of the QTF (depend on the prior non-QTF Nemoh computation)

Tests with unsatisfied requirements will fail.

*******
Running
*******

The binary files of NEMOH1 and NEMOH2 have to be executed following the order provided in :ref:`getting-started:executable files`.

The following steps are for executing the binary files in the command window.

-  Suppose a project directory, *e.g.* ``/NEMOH/projdir/``, that contains all the input files and that is in the same location as the binary directory, i.e ``NEMOH/bin/``.

-  The program can be run depending on your current working directory in the command line. For example, the following commands are possible, with ``binfile`` being an executable file i.e. ``preProc``, etc., as

   -  if you are in the ``projdir``:

      .. code:: bash

         ./../bin/binfile

   -  if you are in the bin:

      .. code:: bash

         ./binfile ./../projdir

Before executing the binary files, the input files are needed. These are described in the next subsection.
Finally, a set of test cases are provided. The results are presented in :ref:`test-cases:test cases` and we provide here the methodology to run those tests. The first possibility is to run each of the cases in the command line, following the above instructions. This applies to both Windows and Linux environments.

To simplify the procedure for Linux platforms, ``Makefile`` is provided in the ``TestCases/`` directory. It is then possible to run the NEMOH1 test cases by executing the following commands in a Terminal (each line being a test case):

.. code:: bash

   make run_1_cylinder
   make run_2_2Bodies
   make run_3_nonsymmetrical
   make run_4_postprocessing
   make run_5_quicktest
   make run_6_box_coarsemesh
   make run_7_Solvers_Check_OC3
   make run_8a_Cylinder_irregfreq

For the QTF test cases, the following commands can be used:

.. code:: bash

   make run_8b_QTF_Cylinder
   make run_9_oc4_semisub
   make run_10a_softwind
   make run_10b_softwind_FS
   make run_11_QTF_OC3_Hywind.

Commands to clean the test cases are also available to clean all the output files. They can apply either to a specific tests case, *e.g.*

.. code:: bash

   make clean_1_cylinder

Or to remove a range of test cases

.. code:: bash

   make clean_all_testsNEMOH1
   make clean_all_testsNEMOH2
   make clean_all_tests

The description and the benchmark results of those test cases are described in :ref:`test-cases:test cases`.


***********************
Supporting Matlab files
***********************

Following Matlab directories, containing a set of functions, are provided in ``matlabRoutines/``,

-  ``NemohWrapper``: This is for running NEMOH executables in MATLAB environment.

-  ``GMSHconverter``: There are two codes, first, for converting body mesh file output from GMSH to NEMOH, DIODORE and HYDROSTAR formats and second, for converting free-surface mesh file output from GMSH to NEMOH and HYDROSTAR formats.

-  ``postproc_testcases``: There are two main codes for plotting results from NEMOH and HYDROSTAR. First, for plotting hydrodynamic coefficients results and second for plotting QTF results. This code can be executed after all data in one specific test cases are obtained.
