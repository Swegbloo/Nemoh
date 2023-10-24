Preamble
========

NEMOH is a Boundary Element Methods (BEM) code dedicated to the computation of wave loads on offshore structures (added mass, radiation damping, diffraction forces, etc.). It has been developed by researchers at Ecole Centrale de Nantes for 30 years. Since its first release in January 2014, NEMOH, as the world’s first open-source BEM code, has been widely used by many researchers, and engineers. Typical use is the estimation of the dynamic response of floating structures or performance assessment of wave energy converters.

Now together with this manual, we gladly release NEMOH v3.0 which includes all recent developments in the first-order module and a new extension module for computing quadratic transfer functions (QTFs). This release will be, to our knowledge, the sole and only open-source BEM software that provides the second-order module.

We sincerely hope that the use of this software, just as the design of it has been, will be fascinating, and challenging for academicians and practitioners in the ocean engineering community. We hope to receive comments and suggestions for further improvement and extension of the software that can be profitable for the community.

Introduction
============

This document is the Manual of NEMOH Software in version v3.0 Released on 2nd December 2022. It serves as a guide for using, installing and running the software.

| NEMOH in its original version is an open-source potential flow boundary element solver for computing first-order hydrodynamic coefficients in the frequency domain. NEMOH v3.0 includes the treatment of irregular frequencies, has an extended module to post-process the first-order hydrodynamic results, and compute complete Quadratic Transfer Functions (QTFs).
| The new developments/extensions in this NEMOH v3.0 are summarised as follows.

-  GNU General Public v3 license is used instead of Apache.

-  Newly written source codes for better readability and coding practice.

-  Green function uses finer interpolation points.

-  The influence coefficient is constructed by applying Gauss-quadrature integration over a panel.

-  Irregular frequency removal method is available.

-  Two new options of the linear system solvers, LU decomposition and GMRES iterative solver, are available for enhanced computational efficiency.

-  | Full difference- and sum-frequencies QTF module.

The present User Manual is organized as follows. Section `2 <#Sec:Descrip_NEMOH>`__ describes briefly the mathematical background and the capabilities of the software; it is advised to read this Section before continuing to the rest of the manual. Section `3 <#Sec:Getstarted>`__ describes the installation procedure, how to run the codes, and includes a description of the input and output files. Section `4 <#Sec:MatlabFiles>`__ describes supporting Matlab files such as a NEMOH wrapper, a mesh converter and the routines to post-process the test cases. Finally, Section `5 <#Sec:Testcase>`__ describes briefly test cases that will show the capabilities of the code and its use. For new NEMOH users, they can be used as the first elements for setting up a new configuration.

.. _`Sec:Descrip_NEMOH`:

Description of NEMOH v3.0
=========================

Modelling and Numerical aspects
-------------------------------

This section provides background information on NEMOH v3.0 with a focus on basic scientific ideas.

| NEMOH v3.0 contains two main modules. First, NEMOH1 solves linear diffraction and radiation problems of wave-structure interaction using a 3-D boundary element method in the frequency domain. Second, NEMOH2, an extended module for computing difference- and sum- frequencies Quadratic Transfer Functions (QTFs) for fixed or floating structures.
| The following subsections describe the underlying modelling and numerical approaches used in NEMOH.

Notations
~~~~~~~~~

As sketched in :numref:`fig:sketch`, we consider fluid domain in the Cartesian coordinate :math:`\boldsymbol x=(\vec{x},z)` with :math:`\vec{x}=(x,y)` the horizontal coordinates perpendicular to the :math:`z` axis in the opposite direction of gravity :math:`\boldsymbol g`. Free-surface boundary :math:`S_F` is defined by the free surface elevation at time :math:`t`, denoted as :math:`\eta(\vec{x},t)` with respect to the mean water level at :math:`z=0`. The fluid velocity potential is denoted as :math:`\Phi(\boldsymbol x,t)` with :math:`\boldsymbol x` in fluid domain :math:`V_{\Omega}`.

.. figure:: figures/Sketch.png
   :name: fig:sketch

   Sketch definition of the system

The floating body has 6 degrees of freedom (DoF), :math:`\boldsymbol\xi=(\boldsymbol{X},\boldsymbol{\theta})` where the positions, :math:`\boldsymbol{X}=(X,Y,Z)` and the orientations, :math:`\boldsymbol{\theta}=(\theta_1,\theta_2,\theta_3)` are determined at the center of gravity (COG). Displacements of points at the hull are specified by a body of vector :math:`\boldsymbol r` with respect to the COG as :math:`\boldsymbol{\mathcal{X}}=\boldsymbol{X}+R(\boldsymbol{r})`. :math:`R` is a rotation operator where :math:`R(\boldsymbol r)\approx \boldsymbol\theta \times \boldsymbol r`. The velocity of the points at the hull is expressed as :math:`\dot{\boldsymbol{\mathcal{X}}}`.

On the body hull :math:`S_B`, the wetted part is defined as a function :math:`z=\zeta(\vec{x},t)`. The normalized normal vector is defined as directed toward the fluid domain, :math:`\boldsymbol n=-\boldsymbol N/|\boldsymbol N|` with :math:`\boldsymbol N=\left(-\nabla_2\zeta,1 \right)` where :math:`\nabla_2` is the two-dimensional gradient in :math:`\vec{x}`. Then the six-dimensional generalized normal vector is defined as :math:`\boldsymbol\nu=(\boldsymbol n,\boldsymbol r \times \boldsymbol n)^T`, with :math:`( )^T` the matrix transpose operator.

Modelling
~~~~~~~~~

NEMOH1, the first-order solver, is based on the following modelling principles:

-  Potential flow theory for an inviscid and incompressible fluid in irrotational motion. The set of equations solved looks for unknowns satisfying free-surface conditions, impermeable bottom condition, diffraction and radiation conditions on the body hull and radiation wave condition in the far field.

-  The harmonic fluid potential is defined as

   .. math::

      \begin{aligned}
      \label{Eq:PhiHarm}
      \Phi(\boldsymbol x,t)=Re\left\lbrace\Phi^{(1)}(\boldsymbol x)e^{-i\omega t}\right\rbrace.
      \end{aligned}

   The total potential, :math:`\Phi`, is the sum of the incident potential, the diffraction potential and the radiation potential.
   The incident potential is defined as, with :math:`k` and :math:`\omega` related with the dispersion relation, :math:`\vec{k}=k(\cos \beta,\sin \beta)`, :math:`\beta` is wave direction and :math:`a` is a unit wave amplitude,

   .. math::

      \begin{aligned}
      \label{Eq:PhiI}
      \Phi_{I}^{(1)}(\boldsymbol x)=-i\frac{a g}{\omega}\frac{\cosh(k(D+z))}{\cosh(kD)} e^{i\vec{k}\cdot \vec{x}}.
      \end{aligned}

   The radiation potential is defined as :math:`\Phi_R(\boldsymbol x,t)=Re\left\lbrace \dot{\boldsymbol\xi}^{(1)}(t) \cdot \boldsymbol\psi(x)\right\rbrace` where :math:`\boldsymbol\psi(\boldsymbol x)` is the normalized vector radiation potential.

-  The three-dimensional linear potential flow problem around arbitrary body condition is reformulated in the Boundary Integral Equation (BIE) and transformed into the two-dimensional problem of the source distribution, :math:`\sigma`, on the body surface, :math:`S_B`, using Green’s second identity and the appropriate Green function, :math:`G(\boldsymbol x,\boldsymbol x')`.

-  The Green function is based on Delhommeau’s formulation and is available for finite and infinite water-depth, see :raw-latex:`\cite{Delhommeau}`.

-  The source distribution depends on the considered boundary condition problem. For each frequency and wave direction, the diffraction source distribution, :math:`\sigma_D(\boldsymbol x)`, depends on the position of the panels while the radiation source distribution, :math:`\sigma_{R_j}(\boldsymbol x)`, depends on the position of the panels and the considered degree of freedom :math:`j`.

-  Then, the BIE for :math:`\boldsymbol x \in S_B`, is expressed as, with flow points :math:`\boldsymbol x` and source points :math:`\boldsymbol x'`,

   .. math::

      \begin{aligned}
      \frac{1}{2}\sigma_{D,R_j}(\boldsymbol x)-\frac{1}{4\pi}\int_{S_B} \partial_n G(\boldsymbol x, \boldsymbol x') \sigma_{D,R_j}(\boldsymbol x') dS'=\mathcal{N}_{D,R_j}(\boldsymbol x). \label{Eq:BIE_source_distribution}
      \end{aligned}

   where :math:`\mathcal{N}(\boldsymbol x)` is the body normal condition. The diffraction normal condition is defined as :math:`\mathcal{N}_D (\boldsymbol x)=-\partial_{n} \Phi_I^{(1)}(\boldsymbol x)`, the normalized radiation condition, :math:`\mathcal{N}_R (\boldsymbol x)=\partial_{n} \Phi_{R_j}(\boldsymbol x)`, with :math:`\Phi_{R_j}(\boldsymbol x)` is the vector component-:math:`j` of the normalized radiation potential :math:`\boldsymbol\psi(\boldsymbol x)`, explicitly :math:`\boldsymbol\psi=(\Phi_{R_1},\Phi_{R_2},\cdots,\Phi_{R_{Ndof}})`.

-  The diffraction potential, :math:`\Phi^{(1)}_{D}`, the normalized radiation potential vector component-:math:`j`, :math:`\Phi_{R_j}` and the corresponding velocities are then computed as follows, for the flow points in the fluid domain :math:`\boldsymbol x \in S_B \cup V_{\Omega_F}`,

   .. math::

      \begin{aligned}
       \label{Eq:BIE_Sol_Pot_Sb}
      \Phi^{(1)}_{D,R_j}(\boldsymbol x)=&-\frac{1}{4\pi}\int_{S_B} G(\boldsymbol x, \boldsymbol x') \sigma_{D,R_j}(\boldsymbol x') dS'\\
      \partial_{\boldsymbol x} \Phi^{(1)}_{D,R_j}(\boldsymbol x)=&\frac{1}{2}\sigma_{D,R_j}(\boldsymbol x)\boldsymbol{n}\delta_{\boldsymbol x \boldsymbol x'}-\frac{1}{4\pi}\int_{S_B} \partial_{\boldsymbol{x}} G(\boldsymbol x, \boldsymbol x') \sigma_{D,R_j}(\boldsymbol x') dS'
      \end{aligned}

   where the Kronecker delta :math:`\delta_{\boldsymbol x \boldsymbol x'}=1` for :math:`\boldsymbol x = \boldsymbol x'`, and :math:`\delta_{\boldsymbol x \boldsymbol x'}=0` otherwise.

-  The hydrodynamic coefficients are then computed as follows, the excitation force is defined as

   .. math::

      \begin{aligned}
      \boldsymbol F_{exc}^{(1)}&=\rho \iint_{S_{B}} -i\omega\left[ \Phi_I^{(1)}+ \Phi_D^{(1)}\right]\boldsymbol\nu dS.
      \end{aligned}

   The added mass matrix and damping coefficient matrix components are computed as

   .. math::

      \begin{aligned}
      M^a_{ij}= -\rho \iint_{S_{B}} \nu_{i} Re \left\lbrace\psi_{R_j} \right\rbrace dS\\
      B_{ij}= -\rho \omega \iint_{S_{B}} \nu_{i} Im \left\lbrace\psi_{R_j} \right\rbrace dS.
      \end{aligned}

-  In post-processing, the radiation damping impulse response matrix function (:math:`\boldsymbol{IRF}(t)`), the infinite frequency added mass matrix (:math:`[\boldsymbol M^a](\infty)`), and the excitation force impulse response vector function (:math:`\boldsymbol{IRF}_{ex}(t)`) are provided. They are computed as,

   .. math::

      \begin{aligned}
      \boldsymbol{IRF}(t)&\approx\frac{2}{\pi}\int_0^{\omega_{max}}[\boldsymbol B](\omega)\cos(\omega t)d\omega, \\
      [\boldsymbol M^a](\infty)&\approx  \frac{1}{N_{\omega}}\sum_{i=1}^{N_{\omega}}[\boldsymbol M^a](\omega_i)+\frac{1}{\omega_i}\int_0^{t_{max}}\boldsymbol{IRF}(t)\sin(\omega_i t)dt\\
      \boldsymbol{IRF}_{exc}(t)&\approx\frac{1}{2\pi}\int_{-\omega_{max}}^{\omega_{max}}\boldsymbol F_{exc}(\omega)e^{-i\omega t}d\omega.
      \end{aligned}

   where :math:`\boldsymbol F_{exc}(-\omega)=\boldsymbol F^*_{exc}(\omega)`. Note that :math:`\omega_{max}` is a user-specified input, for better accuracy of :math:`\boldsymbol{IRF}(t)` make sure that :math:`[\boldsymbol B ](\omega_{max})` has reached an asymptotic value.

-  Response Amplitude Operators (RAO) are obtained by solving the following equation of motion

   .. math::

      \begin{aligned}
      \label{Eq:RAO}
      \left[-[\boldsymbol M+\boldsymbol M^a(\omega)]\omega^2-i\omega[\boldsymbol B(\omega)+\boldsymbol B_{add}]+[\boldsymbol K_h+\boldsymbol K_M]\right]\mathcal{\boldsymbol\xi}(\omega)=\boldsymbol F_{exc}(\omega)
      \end{aligned}

   where :math:`[\boldsymbol B_{add}]` and :math:`[\boldsymbol K_M]` are user-specified additional damping and stiffness matrices.

|  
| NEMOH2, the second-order QTF module, is based on the following principles

-  The second-order loads are composed of the quadratic part and the potential part, the detailed formulation is given in :raw-latex:`\cite{Kurnia22_JH,Kurnia22}`.

-  The quadratic part is based on the near-field method :raw-latex:`\cite{CHEN88}`.

-  The potential part is based on the
   indirect method :raw-latex:`\cite{CHEN88,MOLIN79}`.

Numerical Methods
~~~~~~~~~~~~~~~~~

NEMOH1 uses the following numerical approach:

-  The BIE, Eq. `[Eq:BIE_source_distribution] <#Eq:BIE_source_distribution>`__, is discretised using the constant panel method with quadrilateral mesh. This leads to a linear system with the influence coefficients matrix. The mesh is user-specified with the normal direction towards fluid.

-  Numerical implementation of the Green function is described in :raw-latex:`\cite{Babarit15}`.

-  Free-surface Green function integrands are pre-calculated with the discretized :math:`\omega^2r/g\in [0,100]` with 676 points in a constant scale and :math:`\omega^2(z+z')/g \in [-251,-1.6\, 10^{-6}]` with 130 points in logarithmic scale. A polynomial surface interpolation with the :math:`5^{th}` order Lagrange formula is used for interpolating any values in the specified interval.

-  The specified points for the interpolation of the Green function are finer than in the previous release. However, an option to switch the two different tabulated Green function data is available in the source file ``\Solver\Core\INITIALIZE_GREEN.f90`` with the parameter FLAG_IGREEN=1 or 2, 2 being the default.

-  Influence coefficients, the integration of :math:`\partial_n G(\boldsymbol x, \boldsymbol x')` over a body panel, is computed using Gauss-quadrature integration with a user-input number of Gauss-quadrature points.

-  The source distributions on body panels are then obtained after solving the corresponding linear system.

-  The linear system is solved using a user-choice solver among the available ones, which are Gauss elimination, LU-decomposition (default) and GMRES-iterative solvers.

-  The GMRES solver code :raw-latex:`\cite{GMRES}` from `CERFACS <https://www.cerfacs.fr/algor/Softs/GMRES/index.html>`__ is embedded in NEMOH solver module. For using the GMRES solver, the user has to obtain a license in https://www.cerfacs.fr/algor/Softs/GMRES/license.html.

-  For free-surface piercing bodies problem, the irregular frequencies removal (IRR) method is applied by the user providing lid panels at :math:`z=0`. Then, the extended boundary integral equation will be solved :raw-latex:`\cite{Babarit15,Malenica98}`. As in :raw-latex:`\cite{Malenica98}`, the IRR may be influenced by the input parameter :math:`\epsilon` in ``input_solver.txt`` that shifts the lid panels from :math:`z=0` to :math:`z=-\epsilon d_B` where :math:`d_B` is a maximum horizontal distance of points on the body. :math:`d_B` is computed by the software.

-  RAO in Eq. `[Eq:RAO] <#Eq:RAO>`__ is obtained by applying the inverse matrix using LU-decomposition.

-  The software can solve multi-bodies problems, as well as multi-directional waves.

|  
| NEMOH2 uses the following numerical approach

-  The QTF module can be run only after the first order-hydrodynamic coefficients are computed in NEMOH1.

-  In the potential part, the computation of the free-surface integral is an option:

   -  For the difference-frequency QTFs, it is in general acceptable not to compute the free-surface integral terms.

   -  For the sum-frequency QTFs, it is necessary to compute the free-surface integrals.

-  Important notice: the computation with the free-surface integral still has an issue if the lid body panels exist (cf. IRR method). For now, the user is suggested not to specify the lid body panels in the mesh file input for NEMOH1 computation if he wants to compute the full QTFs with the free surface integral.

-  For the free-surface integral, a quadrilateral free-surface mesh has to be specified.

-  The computation can be done for bi-directional or uni-directional wave for the specified multiple wave direction.

-  QTF computations have not been tested yet for the multi-bodies problem.

|  
| NEMOH related publications to be referred are :raw-latex:`\cite{Babarit15}` for the first order NEMOH and :raw-latex:`\cite{Philippe15,Kurnia22_JH,Kurnia22}` for the QTF module. A publication related with this release is in preparation as in :raw-latex:`\cite{Kurnia23}`.

Units
-----

NEMOH expects all quantities to be expressed in S.I. units: :math:`m, kg, s, rad` (meter, kilogram, seconds, radian, respectively). But some of the phase outputs may be expressed in :math:`deg` or :math:`^{\circ}`, in this case it will be indicated in the file header.

The force unit is [:math:`N`], the moment unit is [:math:`Nm`], added Mass [:math:`kg`], damping coefficient [:math:`kg/s`]. As the force output is normalized with the unit wave amplitude :math:`a` :math:`[m]`, then the normalized force unit is [:math:`N/m`] and the normalized moment is [:math:`N`].

Response amplitude operator for translation motion has unit [:math:`m/m`] and for rotation it is [:math:`deg/m`].

The force quadratic transfer function (QTF) has unit [:math:`N/m^2`] and for the moment QTF it is [:math:`N/m`]. The QTF output is normalized by :math:`\rho g` where the fluid density :math:`\rho,\ [kg/m^3],` and the gravitation constant :math:`g,\ [m/s^2]`.

Software features and capabilities
----------------------------------

.. _`fig:flowchart`:
.. figure:: figures/FlowChart.png

   Global flowchart of NEMOH software

:numref:`fig:flowchart` shows a global overview of the software. There are three main programs: a mesh preprocessor, NEMOH1 and NEMOH2. The program features and capabilities are described as follows.

Mesh Preprocessor
~~~~~~~~~~~~~~~~~

NEMOH mesh preprocessor, the executable file **``mesh``**, is for generating the NEMOH mesh file with a given geometry input file and an input ``Mesh.cal`` file. This **``mesh``** is not a meshing code but allows the user to refine an existing mesh and to calculate properties such as displacement, buoyancy center, and hydrostatic stiffness. It also makes estimates of masses and inertia matrix. The concept with this program is to write by hand a coarse description of the body under consideration in a ``GeomInput`` file and to have **``mesh``** make the refined mesh for NEMOH calculations.

NEMOH1: 1st-order solver
~~~~~~~~~~~~~~~~~~~~~~~~

NEMOH1 solves the first-order potential flow problem. There are four modules: **``preProc``**, **``hydrosCal``**, **``solver``** and **``postProc``**, described as follows.

-  **``preProc``**: processes the input mesh file and generates the body condition for each calculation case (diffraction and radiation). The outputs are used as input for **``solver``**.

-  **``hydrosCal``**: computes hydrostatic parameters, i.e. stiffness matrix and inertia matrix. The output file will be used in the **``postProc``** for computing the RAOs. If the input mesh is generated by the NEMOH mesh preprocessor, **``mesh``**, the hydrostatic parameters are already computed and then it is not necessary to execute this program.

-  **``solver``**: solves the boundary value problems for each problem, diffraction and radiation, defined in the file ``Normalvelocities.dat``, provided by the **``preProc``**.

   -  The influence coefficients matrix is constructed with the infinite/finite depth Green function.

   -  If a finite depth is specified, then the finite depth green function is applied only for :math:`\frac{\omega^2}{g}D<20`, otherwise infinite depth case is applied.

   -  The integration of the Green function on a panel for the influence coefficients is obtained by the Gauss-quadrature integration. The number of Gauss quadrature points is a user input.

   -  The minimum distance, :math:`\epsilon`, between the flow and source points for the influence coefficient computation is user-specified.

   -  The source distributions are then obtained by solving the linear system. There are three options for the solver: Gauss elimination, LU-decomposition and GMRES. If the GMRES solver :raw-latex:`\cite{GMRES}` is used and the target tolerance is not achieved after the maximum number of iterations, the problem is automatically solved by LU-decomposition. License for using GMRES has to be obtained in https://www.cerfacs.fr/algor/Softs/GMRES/license.html.

-  **``postProc``**: post-processes the **``solver``**\ ’s output files. The results are the excitation forces, added mass and damping coefficients. Optionally, the program computes

   -  the radiation damping impulse response function, the infinite frequency added mass and the excitation force impulse response function,

   -  the Kochin coefficient,

   -  the free-surface elevation,

   -  the motion response amplitude operator (RAO). For the RAO computation, additional stiffness matrix :math:`[\boldsymbol K_m]` and additional damping :math:`[\boldsymbol B_{add}]` can be user-specified in the ``Mechanics/`` folder.

NEMOH2: 2nd-order QTF module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

NEMOH2 computes the second-order wave loads that are expressed as Quadratic Transfer Function (QTF). It is suggested to verify the first-order results before running the QTF module. There are three modules in this program: **``QTFpreProc``**, **``QTFsolver``** and **``QTFpostProc``**, described as follows

-  **``QTFpreProc``**: computes the perturbed potential, the total potential, the normalized radiation potential and the corresponding velocities on the body panels, the water-line and the free-surface panels.

   -  The computation on free-surface panels requires possibly long computational time. Then, it is suggested not to compute the free-surface integral for the first execution of NEMOH2. This is controlled by the flag HASFS, which is available in the input file ``Nemoh.cal``.

   -  In general, the free-surface integral may be negligible for the difference-frequency QTFs computation.

   -  The potential on the waterline is rather sensitive with the :math:`\epsilon` value. For default, :math:`\epsilon=0.001`, it can be adjusted in ``input_solver.txt``. The :math:`\epsilon` can be set differently for NEMOH1 and NEMOH2. Further investigation into this is needed.

   -  In case the body lid panels exist, the influence coefficients are affected and give a somewhat larger error for higher frequencies on the free-surface potentials and velocities. This also needs to be investigated.

   -  For now, in the case of full-QTFs computation, the user is suggested not to specify the lid body panels in a mesh file input for NEMOH1 computation.

-  **``QTFsolver``**: computes the quadratic part and the potential part of the second order loads. The free-surface integrals in the potential part QTF are optionnally computed (flag HASFS in ``Nemoh.cal``).

-  **``QTFpostProc``**: adds all the computed QTF parts and produces the total QTF. The option to sum only some parts of the QTF is available in ``Nemoh.cal``.

.. _`Sec:Getstarted`:

Getting-started
===============

This section describes the installation procedure and how to run the codes. It also includes a description of input-output files.

.. _`Sec:Execute`:

Executable files
----------------

| Executable/binary files are provided in the ``Release`` section of the `GitLab project <https://gitlab.com/lheea/Nemoh>`__. They can be used directly without the need to do the compilation procedure described in the next subsection. Windows and Linux executable files are provided. In the case of a compilation on your computer, executable files will be located in the ``bin/`` folder.
| The following executable files are available:

-  NEMOH1: **``mesh``**, **``preProc``**, **``hydrosCal``**, **``solver``**, **``postProc``**,

-  NEMOH2: **``QTFpreProc``**, **``QTFsolver``**, **``QTFpostProc``**.

Note that a Matlab wrapper is provided to use those executable files in a Matlab environment. More details are provided in Sec. `4 <#Sec:MatlabFiles>`__.

Installation
------------

This procedure is intended for a developer to perform the compilation after changing the source code. Users may skip this step.

As a first preliminary step, it is necessary to install the following external libraries that are used by NEMOH:

-  BLAS, https://netlib.org/blas/

-  LAPACK, https://netlib.org/lapack/

Windows platform
~~~~~~~~~~~~~~~~

An updated manual will provide the details of the compilation on a Windows platform.

Linux platform
~~~~~~~~~~~~~~

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

Running
-------

The binary files of NEMOH1 and NEMOH2 have to be executed following the order provided in Sec. `3.1 <#Sec:Execute>`__.

The following steps are for executing the binary files in the command window.

-  Suppose a project directory, *e.g.* ``/NEMOH/projdir/``, that contains all the input files and that is in the same location as the binary directory, i.e ``NEMOH/bin/``.

-  The program can be run depending on your current working directory in the command line. For example, the following commands are possible, with **``binfile``** being an executable file i.e. **``preProc``**, etc., as

   -  if you are in the ``projdir``:

      .. code:: bash

         ./../bin/binfile

   -  if you are in the bin:

      .. code:: bash

         ./binfile ./../projdir

| Before executing the binary files, the input files are needed. These are described in the next subsection.
| Finally, a set of test cases are provided. The results are presented in Sec. `5 <#Sec:Testcase>`__ and we provide here the methodology to run those tests. The first possibility is to run each of the cases in the command line, following the above instructions. This applies to both Windows and Linux environments.

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

The description and the benchmark results of those test cases are described in Sec. `5 <#Sec:Testcase>`__.

Input/Output
------------

Following is the list of the user’s input files and the output files for each executable file:

-  **``mesh``**

   -  User’s input files: ``projdir/Mesh.cal``, ``projdir/geomInput``,

   -  Output files: ``projdir/meshfile``, ``projdir/mesh/[*.dat, *.tec]``,

-  **``preProc``**

   -  User’s input files: ``projdir/Nemoh.cal``, ``projdir/meshfile``,

   -  | Output files: ``projdir/Normalvelocities.dat``,
      | ``projdir/results/[FKForce.dat, FKForce.tec, index.dat]``,
      | ``projdir/mesh/[L10.dat, L12.dat]``,

-  **``hydrosCal``**

   -  User’s input files: ``projdir/Nemoh.cal``, ``projdir/mesh.cal``,

   -  Output files: ``projdir/mesh/[*.dat, *.tec]``, ``projdir/Mechanics/``,

-  **``solver``**

   -  User’s input files: ``projdir/Nemoh.cal``, ``projdir/input_solver.txt``

   -  | Output files: ``projdir/results/Forces.dat``,
      | ``projdir/results/sources`` (if QTF will be computed, indicated in ``Nemoh.cal``)

-  **``postProc``**

   -  | User’s input files: ``projdir/Nemoh.cal``,
      | ``projdir/Mechanics/[Km.dat,Badd.dat]``

   -  | Output files: ``projdir/results/[ExcitationForce.tec,``
      | ``DiffractionForce.tec,RadiationCoefficients.tec,``\ :math:`\cdots`\ ``]``,
      | ``projdir/Motion/RAO.dat``,

-  **``QTFpreProc``**

   -  | User’s input files: ``projdir/Nemoh.cal``, ``projdir/FSmeshfile`` (If the free-surface integral, HASFS flag, is computed),
      | ``projdir/Mechanics/[Km.dat,Badd.dat]``

   -  Output files: ``projdir/QTFPreprocOut/*.bin``,

-  **QTFsolver**

   -  User’s input files: ``projdir/Nemoh.cal``

   -  Output files: ``projdir/results/QTF/*.dat``,

-  **QTFpostproc**

   -  User’s input files: ``projdir/Nemoh.cal``

   -  Output files: ``projdir/results/QTF/[OUT_QTFM_N.dat,OUT_QTFP_N.dat]``

|  
| As summary, following files are needed for the input, some depends on user-specified choice in ``Nemoh.cal``:

-  ``Nemoh.cal`` contains all NEMOH computation parameters

-  ``Mesh.cal`` contains information of ``geomInput`` file. It is an input for **``mesh``** and **``hydrosCal``**.

-  a ``meshfile``, input for **``preProc``**, or ``geomInput`` file, input for **``mesh``**

-  ``input_solver.txt`` contains **``solver``** parameters

-  ``Km.dat`` and ``Badd.dat``, are the additional stiffness and damping matrices. These optional input are for **``postProc``**/**``QTFpreProc``**

-  ``FSmeshfile`` contains the free-surface mesh if contrib=3 in ``Nemoh.cal``.

|  
| NEMOH produces the following main output files, some depending on user-specified choice in ``Nemoh.cal``,

-  hydrostatic files: inertia and stiffness matrices,

-  hydrodynamics coefficients: Froude-Krylov force, excitation force, added-mass, damping coefficient,

-  Kochin function, free-surface elevation,

-  Response Amplitude Operator (RAO),

-  Total difference- and sum-frequencies QTFs.

|  
| Detail descriptions of the input/output files are discussed in the next subsections.

User’s input files
~~~~~~~~~~~~~~~~~~

.. _`fig:NemohCal`:
.. figure:: figures/NemohCal.png

   ``Nemoh.cal`` input file

``Nemoh.cal``: contains all computation parameters with the format as in :numref:`fig:NemohCal`. The following parameters have to be specified:

-  *Environment*: fluid density, gravity constant, water-depth and wave reference point. Specify :math:`0.` for infinite water depth case.

-  *Description of floating bodies*: number of bodies, name of the ``meshfile``, number of points and number of panels, number of degrees of freedom, motion description, number of resulting generalized forces and its description.

   -  ``meshfile`` has to be provided in the folder ``projdir/``

   -  Number of points and panels correspond with the data in the ``meshfile``

   -  For each motion and resulting generalized force, 7 parameters have to be specified in a row. The first parameter, 1 for translation motion/force, 2 for rotation motion/force. The second to fourth parameters is a unit vector of motion/force, 1 0 0 for surge/roll, 0 1 0 for sway/roll, 0 0 1 for heave/yaw. The fifth to seventh parameters are the reference point coordinate :math:`(x,y,z)`.

   -  In the multibody cases, all the parameters have to be specified in rows for each body.

-  *Load case to be solved*: contains wave frequency and wave direction parameter.

   -  Four wave frequency parameters have to be specified in a row. First, the frequency type, 1 for radial-frequency [rad/s], 2 for frequency [Hz] and 3 for period [s]. The second to fourth parameters are the number of frequencies/periods, and their minimum and maximum values.

   -  In the case of QTF computation, it is suggested that the computed first-order hydrodynamic coefficients in NEMOH1 cover all difference-frequencies and sum-frequencies intervals for the QTF computation. The suggested radial frequency interval is :math:`\omega \in [\Delta \omega, \omega_{max}]` with a step :math:`\Delta \omega`, :math:`\omega_{max}=N_{\omega}\Delta \omega` and :math:`N_{\omega}` is the total number of radial frequencies. The :math:`\omega_{max}` should be chosen as the maximum computed sum-frequencies, :math:`\omega_1+\omega_2`.

   -  The wave direction parameters are the number of directions, and the minimum and maximum angle directions [deg].

-  *Post processing* contains parameters for Impulse Response Functions, pressure, Kochin, free surface elevation, RAO and frequency type output.

   -  The IRFs computation requires 3 parameters; a switch [0 is not calculated, 1 calculated], a time-step and a maximum time.

   -  A switch pressure parameter [0 is not calculated, 1 calculated].

   -  The Kochin parameters are a number of directions (specify 0 if it is not calculated), minimum and maximum values of angle directions [deg].

   -  The free surface parameters are the number of points in :math:`x`-direction (specify 0 if it is not calculated) and :math:`y`-direction, the dimension of the domain in :math:`x` and :math:`y` directions [m].

   -  A switch for RAO computation [0 is not calculated and 1 is calculated]. If QTF will be computed, the RAO has to be computed and then the switch has to be set to 1.

   -  Frequency output option: 1 for the radial frequency [rad/s], 2 for the frequency [Hz] and 3 for the period [s].

-  *Quadratic Transfer Function (QTF)* contains following parameters

   -  A flag to decide if QTFs are computed (1) or not (0). If flag=1, then the NEMOH1 module produces the source-distribution file for each problem, which is saved in ``/projdir/results/source``.

   -  The frequency parameters are provided only in radial frequency [rad/s] under the form of: the number of frequencies, and minimum and maximum values. The values are not necessarily the same as the input in NEMOH 1 but should be within the range of NEMOH1 input, then an interpolation may be applied.

   -  A flag to choose between uni-directional (0) and bi-directional (1) computations of QTFs. If multi-directions are computed in NEMOH1, all the direction interactions will be computed for the bi-directional case. In the uni-directional case, only the same-direction interaction, :math:`\beta_1=\beta_2`, will be computed.

   -  Contribution parameter: 1 computes only the quadratic terms of QTFs (referenced as DUOK), 2 computes the quadratic and the body force contribution in the potential QTFs referenced as DUOK+HASBO), 3 includes the computation of the free-surface integrals in the finite domain and semi-infinite domain (this is referenced as DUOK+HASBO+HASFS+ASYMP).

   -  If Contrib=3, a free-surface mesh file name has to be specified and the file is placed in ``/projdir/mesh``. Type NA if it is not applicable, Contrib\ :math:`<3`.

   -  If Contrib=3, specify the free-surface parameters: an external radius, :math:`R_e` is the maximum radius as in the free-surface mesh, a number of discretized radius points in between the waterline :math:`R_b` and the external radius, :math:`[R_b,R_e]`, and the number of Bessel functions. The number of Bessel functions is used during the computation in the semi-infinite integral, specify 30 as a default value.

   -  Switch 1 for adding to the quadratic QTF (DUOK) the Hydrostatic terms of the quadratic first-order motion, :math:`-[\boldsymbol K] \tilde{\boldsymbol\xi}^{(2)}`, where, with :math:`z_G` is the vertical component of CoG,

      .. math::

         \begin{aligned}
         \tilde{\boldsymbol\xi}^{(2)}=[0,0,z_G(\theta_1^{(1)2}+\theta_2^{(1)2})/2,\theta_2^{(1)}\theta_3^{(1)}/2,-\theta_3^{(1)}\theta_1^{(1)}/2,0]^T.
         \end{aligned}

      Note that this term is optional and needed only in QTFsolver. In other software this term is not always included, *e.g.* HYDROSTAR :raw-latex:`\cite{HYDROSTAR}` does not included it.

   -  Frequency output option: 1 for the radial frequency [rad/s], 2 for the frequency [Hz] and 3 for the period [s].

   -  In **``QTFpostproc``**, QTFs total is calculated with summation of all the terms. Option to exclude/include the terms are available with the corresponding switch for DUOK, HASBO and HASFS+ASYMP terms. Switch 1 to include, 0 to exclude. If Contrib=2, then HASFS+ASYMP switch has to be 0.

 

| ``meshfile``: contains all the mesh information with a format as shown in Table `1 <#tab:meshfile>`__. Lid panels (:math:`z=0`) of the structure may be included in this file to activate the irregular frequencies removal method. This mesh file may be generated by NEMOH **``mesh``** preprocessor or by an external mesh generator.
  External mesh generators, *e.g.* the open-source software GMSH :raw-latex:`\cite{gmsh}`, may be used to generate mesh files but they must be adapted to the NEMOH format. A Matlab file for converting GMSH mesh file to the NEMOH format is provided in the dedicated repository. The Matlab file will be described in the next section.

.. container:: center

   .. container::
      :name: tab:meshfile

      .. table:: ``meshfile`` format

         +---+-------------+-------------+-------------+---------------+
         | 2 | 1           |             |             | First column  |
         |   |             |             |             | must be a 2.  |
         +---+-------------+-------------+-------------+---------------+
         |   |             |             |             | Second column |
         |   |             |             |             | is 1 if half  |
         |   |             |             |             | symmetric     |
         |   |             |             |             | body mesh,    |
         |   |             |             |             | about         |
         |   |             |             |             | (             |
         |   |             |             |             | :math:`xOz`), |
         +---+-------------+-------------+-------------+---------------+
         |   |             |             |             | specified, 0  |
         |   |             |             |             | otherwise.    |
         +---+-------------+-------------+-------------+---------------+
         | 1 | :math:`x_1` | :math:`y_1` | :math:`z_1` | Table of      |
         |   |             |             |             | nodes. First  |
         |   |             |             |             | column is     |
         |   |             |             |             | node ID.      |
         +---+-------------+-------------+-------------+---------------+
         |   |             |             |             | Other columns |
         |   |             |             |             | are the node  |
         |   |             |             |             | coordinates   |
         |   |             |             |             | :m            |
         |   |             |             |             | ath:`(x,y,z)` |
         +---+-------------+-------------+-------------+---------------+
         | : | :           | :           | :           | All nodes     |
         |   |             |             |             | coordinated   |
         |   |             |             |             | listed in the |
         |   |             |             |             | rows          |
         +---+-------------+-------------+-------------+---------------+
         | 0 | 0.          | 0.          | 0.          | Last line of  |
         |   |             |             |             | table of      |
         |   |             |             |             | nodes         |
         +---+-------------+-------------+-------------+---------------+
         | 1 | 2           | 3           | 4           | Table of      |
         |   |             |             |             | co            |
         |   |             |             |             | nnectivities. |
         |   |             |             |             | Number of     |
         |   |             |             |             | node IDs      |
         +---+-------------+-------------+-------------+---------------+
         | : | :           | :           | :           | C             |
         |   |             |             |             | onnectivities |
         |   |             |             |             | in each panel |
         |   |             |             |             | listed in the |
         |   |             |             |             | rows          |
         +---+-------------+-------------+-------------+---------------+
         | 0 | 0           | 0           | 0           | Last line of  |
         |   |             |             |             | table of      |
         |   |             |             |             | c             |
         |   |             |             |             | onnectivities |
         +---+-------------+-------------+-------------+---------------+

| ``geomInput``: contain coarse description of mesh, that are number of nodes, number of panels, table of nodes and table of connectivities. The input file has to follow the format as shown in Table `2 <#tab:geomInput>`__.

.. container:: center

   .. container::
      :name: tab:geomInput

      .. table:: ``geomInput`` file format

         +-------------+-------------+-------------+---+---------------+
         | 100         |             |             |   | Total number  |
         |             |             |             |   | of nodes.     |
         +=============+=============+=============+===+===============+
         | 25          |             |             |   | Total number  |
         |             |             |             |   | of panels.    |
         +-------------+-------------+-------------+---+---------------+
         | :math:`x_1` | :math:`y_1` | :math:`z_1` |   | Table of      |
         |             |             |             |   | nodes.        |
         +-------------+-------------+-------------+---+---------------+
         |             |             |             |   | The node      |
         |             |             |             |   | coordinates   |
         |             |             |             |   | :m            |
         |             |             |             |   | ath:`(x,y,z)` |
         +-------------+-------------+-------------+---+---------------+
         | :           | :           | :           |   | All nodes     |
         |             |             |             |   | coordinated   |
         |             |             |             |   | listed in the |
         |             |             |             |   | rows          |
         +-------------+-------------+-------------+---+---------------+
         | 1           | 2           | 3           | 4 | Table of      |
         |             |             |             |   | co            |
         |             |             |             |   | nnectivities. |
         |             |             |             |   | Number of     |
         |             |             |             |   | node IDs      |
         +-------------+-------------+-------------+---+---------------+
         | :           | :           | :           | : | C             |
         |             |             |             |   | onnectivities |
         |             |             |             |   | in each panel |
         |             |             |             |   | listed in the |
         |             |             |             |   | rows          |
         +-------------+-------------+-------------+---+---------------+

| ``Mesh.cal:`` contains mesh and environmental parameters with a format as in Table `3 <#tab:meshcal>`__. This file is used as input for **``mesh``** and **``hydroCal``**. All the parameters are used in **``mesh``**. Only center of gravity, water density, and gravity are used in **``hydroCal``**.

.. container:: center

   .. container::
      :name: tab:meshcal

      .. table:: ``Mesh.cal`` file format

         +----------------+----+----+---+---------------------------+
         | geomInput_name |    |    |   | Name of the geomInput     |
         |                |    |    |   | file.                     |
         +================+====+====+===+===========================+
         | 0              |    |    |   | 1 if a half symmetric     |
         |                |    |    |   | body mesh, about          |
         |                |    |    |   | (:math:`xOz`), specified. |
         +----------------+----+----+---+---------------------------+
         | 0.             | 0. |    |   | Possible translation      |
         |                |    |    |   | about x axis (first       |
         |                |    |    |   | column)                   |
         +----------------+----+----+---+---------------------------+
         |                |    |    |   | and y axis (second        |
         |                |    |    |   | column)                   |
         +----------------+----+----+---+---------------------------+
         | 0.             | 0. | -7 |   | Coordinates of gravity    |
         |                |    |    |   | centre                    |
         +----------------+----+----+---+---------------------------+
         | 500.           |    |    |   | Target for the number of  |
         |                |    |    |   | panels in refined mesh    |
         +----------------+----+----+---+---------------------------+
         | 2.             |    |    |   |                           |
         +----------------+----+----+---+---------------------------+
         | 0.             |    |    |   |                           |
         +----------------+----+----+---+---------------------------+
         | 1.             |    |    |   |                           |
         +----------------+----+----+---+---------------------------+
         | 1025           |    |    |   | water density             |
         |                |    |    |   | :math:`(kg/m^3)`          |
         +----------------+----+----+---+---------------------------+
         | 9.81           |    |    |   | gravity :math:`(m/s^2)`   |
         +----------------+----+----+---+---------------------------+

``input_solver.txt`` contains solver parameters with format as in Table `4 <#tab:input_solver>`__. The parameters are described as follows.

-  Number of Gauss Quadrature points, :math:`N^2`, is used for the surface integration in the influence coefficients. User specifies an integer value of :math:`N\in [1,4]`, default :math:`N=2`.

-  Minimum z of flow and source points is defined with a factor :math:`\epsilon_{zmin}` multiplied by the maximal horizontal distance between two point of the mesh, default :math:`\epsilon_{zmin}=0.001`.

-  Three linear-system solvers are available; 1 Gauss elimination, 2 LU Decomposition, 3 GMRES iterative solver.

-  | If GMRES solver is chosen then the three parameters, the restart parameter, the relative tolerance and the maximum number of iterations, have to be specified. If the tolerance is not achieved after the maximum iteration exceeded then LU decomposition solves the system directly.

.. container:: center

   .. container::
      :name: tab:input_solver

      .. table:: ``input_solver.txt`` file format

         +-------+------+------+---------------------------------------------+
         | 2     |      |      | Gauss quadrature (GQ) surface integration,  |
         |       |      |      | :math:`N^2` GQ Nodes,                       |
         +-------+------+------+---------------------------------------------+
         |       |      |      | specify N=[1,4]                             |
         +-------+------+------+---------------------------------------------+
         | 0.001 |      |      | eps_zmin for determine minimum z of flow    |
         |       |      |      | and                                         |
         +-------+------+------+---------------------------------------------+
         |       |      |      | source points of panel.                     |
         +-------+------+------+---------------------------------------------+
         | 1     |      |      | Solver option: 0 GAUSS ELIM., 1 LU DECOMP., |
         |       |      |      | 2 GMRES                                     |
         +-------+------+------+---------------------------------------------+
         | 10    | 1e-5 | 1000 | GMRES parameters: restart parameter, Rel    |
         |       |      |      | Tol, max iter                               |
         +-------+------+------+---------------------------------------------+

| ``Km.dat`` and ``Badd.dat`` are additional stiffness matrix and damping coefficient matrix. The files contains the matrix components with size :math:`(Nbody\cdot Nradiation)\times (Nbody\cdot Nradiation)`.
| ``FSmeshfile`` contains all the free-surface mesh information with a format as shown in Table `5 <#tab:FSmeshfile>`__. Quadrilateral panels discretized free-surface area in between the body waterline, :math:`R_B`, and the exterior radius :math:`R_e`. Waterline on :math:`R_B` and :math:`R_e` has to discretized by line segments.

.. container:: center

   .. container::
      :name: tab:FSmeshfile

      .. table:: ``FSmeshfile`` format (Free surface mesh file)

         +------+-------------+-------------+-------------+--------------+
         | 1    | 5000        | 4900        | 400         | This row     |
         |      |             |             |             | contais the  |
         |      |             |             |             | free-surface |
         |      |             |             |             | computation  |
         |      |             |             |             | parameters.  |
         +------+-------------+-------------+-------------+--------------+
         |      |             |             |             | First column |
         |      |             |             |             | is 1 if half |
         |      |             |             |             | symmetric    |
         |      |             |             |             | free surface |
         |      |             |             |             | mesh         |
         +------+-------------+-------------+-------------+--------------+
         |      |             |             |             | specified, 0 |
         |      |             |             |             | otherwise.   |
         +------+-------------+-------------+-------------+--------------+
         |      |             |             |             | Column 2-4   |
         |      |             |             |             | are Number   |
         |      |             |             |             | of points,   |
         |      |             |             |             | Number of    |
         |      |             |             |             | panels,      |
         +------+-------------+-------------+-------------+--------------+
         |      |             |             |             | Number of    |
         |      |             |             |             | segmented    |
         |      |             |             |             | waterline,   |
         |      |             |             |             | r            |
         |      |             |             |             | espectively. |
         +------+-------------+-------------+-------------+--------------+
         | 1    | :math:`x_1` | :math:`y_1` | :math:`z_1` | Table of     |
         |      |             |             |             | nodes. First |
         |      |             |             |             | column is    |
         |      |             |             |             | node ID.     |
         +------+-------------+-------------+-------------+--------------+
         |      |             |             |             | Other        |
         |      |             |             |             | columns are  |
         |      |             |             |             | the node     |
         |      |             |             |             | coordinates  |
         |      |             |             |             | :ma          |
         |      |             |             |             | th:`(x,y,z)` |
         +------+-------------+-------------+-------------+--------------+
         | :    | :           | :           | :           | All nodes    |
         |      |             |             |             | coordinated  |
         |      |             |             |             | listed in    |
         |      |             |             |             | the rows     |
         +------+-------------+-------------+-------------+--------------+
         | 0    | 0.          | 0.          | 0.          | Last line of |
         |      |             |             |             | table of     |
         |      |             |             |             | nodes        |
         +------+-------------+-------------+-------------+--------------+
         | 1    | 2           | 3           | 4           | Table of     |
         |      |             |             |             | co           |
         |      |             |             |             | nnectivities |
         |      |             |             |             | in a panel.  |
         +------+-------------+-------------+-------------+--------------+
         |      |             |             |             | Number of    |
         |      |             |             |             | node IDs     |
         +------+-------------+-------------+-------------+--------------+
         | :    | :           | :           | :           | Co           |
         |      |             |             |             | nnectivities |
         |      |             |             |             | in each      |
         |      |             |             |             | panel listed |
         |      |             |             |             | in the rows  |
         +------+-------------+-------------+-------------+--------------+
         | 4901 | 4902        |             |             | Table of     |
         |      |             |             |             | co           |
         |      |             |             |             | nnectivities |
         |      |             |             |             | in a         |
         |      |             |             |             | segmented    |
         |      |             |             |             | waterline.   |
         +------+-------------+-------------+-------------+--------------+
         |      |             |             |             | Number of    |
         |      |             |             |             | node IDs     |
         +------+-------------+-------------+-------------+--------------+
         | :    | :           |             |             | Co           |
         |      |             |             |             | nnectivities |
         |      |             |             |             | in each line |
         |      |             |             |             | listed in    |
         |      |             |             |             | the rows     |
         +------+-------------+-------------+-------------+--------------+
         | 0    | 0           | 0           | 0           | Last line of |
         |      |             |             |             | table of     |
         |      |             |             |             | co           |
         |      |             |             |             | nnectivities |
         +------+-------------+-------------+-------------+--------------+

Output files
~~~~~~~~~~~~

Hydrostatic output files such as inertia and stiffness matrices are produced by **``mesh``**, if ``geomInput`` is prescribed, or by **``hydroCal``**, if ``meshfile`` is prescribed. The files contain the matrix components with size :math:`(Nbody\cdot Nradiation)\times (Nbody\cdot Nradiation)`.

The following hydrodynamic coefficients are produced in Tecplot format, which can be opened by the Tecplot program or by a simple text-editor program,

-  ``FKForce.tec``, ``DiffractionForce.tec`` and ``ExcitationForce.tec`` are the output files of the Froude-Krylov, the diffraction and the excitation forces respectively. The output file format is given in Table `6 <#tab:WaveForce>`__. The file contains the absolute value and the phase [deg] of the force for each ’frequency’ :math:`f`. The force is given for each specified force axis (i.e. surge, heave, pitch) for each body. The ’frequency’ is given based on the chosen type, [rad/s, Hz, s], of the post-processing parameter in ``Nemoh.cal``, except the Froude-Krylov force, which is only in the radial frequency [rad/s].

-  ``RadiationCoefficients.tec`` is the output file for added mass and damping coefficients with format as in Table `7 <#tab:addedmass_damping_coeffs>`__. The radiation coefficients are given for each :math:`DoF`, each force axis and for each frequency. The frequency is given based on the chosen ’frequency’ type, [rad/s, Hz, s], of the post-processing parameter in ``Nemoh.cal``.

| The hydrodynamic coefficients are also produced in the *.dat* files, i.e. *CA.dat* for the damping coefficients, *CM.dat* for the added mass coefficients, *Fe.dat* for the excitation force and *FKForce.dat* for the excitation force. The frequency type of the output files is only radial frequency [rad/s]. These output files are used as input files for the QTF module.

.. container:: center

   .. container::
      :name: tab:WaveForce

      .. table:: Output file format of Froude-Krylov, diffraction and excitation forces

         +---------+---------+---------+---------+---------+---------+---------+
         | :mat    | :mat    | :m      | :math:` | :math:` | :mat    | :m      |
         | h:`f_1` | h:`|F_1 | ath:`\a | \cdots` | \cdots` | h:`|F_{ | ath:`\a |
         |         | (f_1)|` | ngle F_ |         |         | Ninteg} | ngle F_ |
         |         |         | 1(f_1)` |         |         | (f_1)|` | {Ninteg |
         |         |         |         |         |         |         | }(f_1)` |
         +=========+=========+=========+=========+=========+=========+=========+
         | :mat    | :mat    | :m      | :math:` | :math:` | :mat    | :m      |
         | h:`f_2` | h:`|F_1 | ath:`\a | \cdots` | \cdots` | h:`|F_{ | ath:`\a |
         |         | (f_2)|` | ngle F_ |         |         | Ninteg} | ngle F_ |
         |         |         | 1(f_2)` |         |         | (f_2)|` | {Ninteg |
         |         |         |         |         |         |         | }(f_2)` |
         +---------+---------+---------+---------+---------+---------+---------+
         | :math:` | :math:` | :math:` | :math:` | :math:` | :math:` | :math:` |
         | \vdots` | \vdots` | \vdots` | \vdots` | \vdots` | \vdots` | \vdots` |
         +---------+---------+---------+---------+---------+---------+---------+
         | :math:` | :math:` | :math   | :math:` | :math:` | :math:` | :math   |
         | f_{Nf}` | |F_1(f_ | :`\angl | \cdots` | \cdots` | |F_{Nin | :`\angl |
         |         | {Nf})|` | e F_1(f |         |         | teg}(f_ | e F_{Ni |
         |         |         | _{Nf})` |         |         | {Nf})|` | nteg}(f |
         |         |         |         |         |         |         | _{Nf})` |
         +---------+---------+---------+---------+---------+---------+---------+

.. container:: center

   .. container::
      :name: tab:addedmass_damping_coeffs

      .. table:: Output file format of the radiation coefficients

         +---------+---------+---------+---------+---------+---------+---------+
         | :mat    | :math:` | :math   | :math:` | :math:` | :math   | :ma     |
         | h:`f_1` | M^a_{11 | :`B_{11 | \cdots` | \cdots` | :`M^a_{ | th:`B_{ |
         |         | }(f_1)` | }(f_1)` |         |         | 1Ninteg | 1Ninteg |
         |         |         |         |         |         | }(f_1)` | }(f_1)` |
         +=========+=========+=========+=========+=========+=========+=========+
         | :mat    | :math:` | :math   | :math:` | :math:` | :math   | :ma     |
         | h:`f_2` | M^a_{11 | :`B_{11 | \cdots` | \cdots` | :`M^a_{ | th:`B_{ |
         |         | }(f_2)` | }(f_2)` |         |         | 1Ninteg | 1Ninteg |
         |         |         |         |         |         | }(f_2)` | }(f_2)` |
         +---------+---------+---------+---------+---------+---------+---------+
         | :math:` | :math:` | :math:` | :math:` | :math:` | :math:` | :math:` |
         | \vdots` | \vdots` | \vdots` | \vdots` | \vdots` | \vdots` | \vdots` |
         +---------+---------+---------+---------+---------+---------+---------+
         | :       | :mat    | :m      | :math:` | :math:` | :m      | :math:` |
         | math:`f | h:`M^a_ | ath:`B_ | \cdots` | \cdots` | ath:`M^ | B_{1Nin |
         | _{N_f}` | {11}(f_ | {11}(f_ |         |         | a_{1Nin | teg}(f_ |
         |         | {N_f})` | {N_f})` |         |         | teg}(f_ | {N_f})` |
         |         |         |         |         |         | {N_f})` |         |
         +---------+---------+---------+---------+---------+---------+---------+
         |         |         |         |         |         |         |         |
         +---------+---------+---------+---------+---------+---------+---------+
         | :mat    | :math:` | :math   | :math:` | :math:` | :math   | :ma     |
         | h:`f_1` | M^a_{21 | :`B_{21 | \cdots` | \cdots` | :`M^a_{ | th:`B_{ |
         |         | }(f_1)` | }(f_1)` |         |         | 2Ninteg | 2Ninteg |
         |         |         |         |         |         | }(f_1)` | }(f_1)` |
         +---------+---------+---------+---------+---------+---------+---------+
         | :math:` | :math:` | :math:` | :math:` | :math:` | :math:` | :math:` |
         | \vdots` | \vdots` | \vdots` | \vdots` | \vdots` | \vdots` | \vdots` |
         +---------+---------+---------+---------+---------+---------+---------+
         | :       | :mat    | :m      | :math:` | :math:` | :m      | :math:` |
         | math:`f | h:`M^a_ | ath:`B_ | \cdots` | \cdots` | ath:`M^ | B_{2Nin |
         | _{N_f}` | {21}(f_ | {21}(f_ |         |         | a_{2Nin | teg}(f_ |
         |         | {N_f})` | {N_f})` |         |         | teg}(f_ | {N_f})` |
         |         |         |         |         |         | {N_f})` |         |
         +---------+---------+---------+---------+---------+---------+---------+
         | :math:` | :math:` | :math:` | :math:` | :math:` | :math:` | :math:` |
         | \vdots` | \vdots` | \vdots` | \vdots` | \vdots` | \vdots` | \vdots` |
         +---------+---------+---------+---------+---------+---------+---------+
         | :       | :ma     | :       | :math:` | :math:` | :       | :math:  |
         | math:`f | th:`M^a | math:`B | \cdots` | \cdots` | math:`M | `B_{N_{ |
         | _{N_f}` | _{N_{Do | _{N_{Do |         |         | ^a_{N_{ | DoF}Nin |
         |         | F}1}(f_ | F}1}(f_ |         |         | DoF}Nin | teg}(f_ |
         |         | {N_f})` | {N_f})` |         |         | teg}(f_ | {N_f})` |
         |         |         |         |         |         | {N_f})` |         |
         +---------+---------+---------+---------+---------+---------+---------+

| ``RAO.dat`` is the output file of the response amplitude operator with the file format as in Table. `8 <#tab:RAO>`__. The output file gives the absolute value and the phase of RAO for each degree of freedom and each frequency. The frequency is given based on the chosen ’frequency’ type, [rad/s, Hz, s], of the post-processing parameter in ``Nemoh.cal``. Only radial frequency output file will be produced in the case of the QTF computed.

.. container:: center

   .. container::
      :name: tab:RAO

      .. table:: Output file format of ``RAO.dat``

         +---------+---------+---------+---------+---------+---------+---------+
         | :mat    | :math:  | :math:` | :math:  | :mat    | :math:` | :mat    |
         | h:`f_1` | `|\xi_1 | \cdots` | `|\xi_6 | h:`\ang | \cdots` | h:`\ang |
         |         | (f_1)|` |         | (f_1)|` | le \xi_ |         | le \xi_ |
         |         |         |         |         | 1(f_1)` |         | 6(f_1)` |
         +=========+=========+=========+=========+=========+=========+=========+
         | :math:` | :math:` | :math:` | :math:` | :math:` | :math:` | :math:` |
         | \vdots` | \vdots` | \vdots` | \vdots` | \vdots` | \vdots` | \vdots` |
         +---------+---------+---------+---------+---------+---------+---------+
         | :       | :ma     | :math:` | :ma     | :       | :math:` | :       |
         | math:`f | th:`|\x | \cdots` | th:`|\x | math:`\ | \cdots` | math:`\ |
         | _{N_f}` | i_1(f_{ |         | i_6(f_{ | angle \ |         | angle \ |
         |         | N_f})|` |         | N_f})|` | xi_1(f_ |         | xi_6(f_ |
         |         |         |         |         | {N_f})` |         | {N_f})` |
         +---------+---------+---------+---------+---------+---------+---------+

| ``IRF.tec`` and ``IRF_excForce.tec`` are the impulse response functions for the radiation damping and the excitation force, respectively. The radiation damping IRF has the file format as in Table `9 <#tab:IRF>`__ and the excitation force IRF as in Table `10 <#tab:IRFExcF>`__.

.. container:: center

   .. container::
      :name: tab:IRF

      .. table:: Output file format of ``IRF.tec``

         +---------+---------+---------+---------+---------+---------+---------+
         | :mat    | :ma     | :math:` | :math:` | :math:` | :       | :math   |
         | h:`t_1` | th:`M^a | IRF_{11 | \cdots` | \cdots` | math:`M | :`IRF_{ |
         |         | _{11}(\ | }(t_1)` |         |         | ^a_{1Ni | 1Ninteg |
         |         | infty)` |         |         |         | nteg}(\ | }(t_1)` |
         |         |         |         |         |         | infty)` |         |
         +=========+=========+=========+=========+=========+=========+=========+
         | :mat    | :ma     | :math:` | :math:` | :math:` | :       | :math   |
         | h:`t_2` | th:`M^a | IRF_{11 | \cdots` | \cdots` | math:`M | :`IRF_{ |
         |         | _{11}(\ | }(t_2)` |         |         | ^a_{1Ni | 1Ninteg |
         |         | infty)` |         |         |         | nteg}(\ | }(t_2)` |
         |         |         |         |         |         | infty)` |         |
         +---------+---------+---------+---------+---------+---------+---------+
         | :math:` | :math:` | :math:` | :math:` | :math:` | :math:` | :math:` |
         | \vdots` | \vdots` | \vdots` | \vdots` | \vdots` | \vdots` | \vdots` |
         +---------+---------+---------+---------+---------+---------+---------+
         |         |         |         |         |         |         |         |
         +---------+---------+---------+---------+---------+---------+---------+
         | :mat    | :ma     | :math:` | :math:` | :math:` | :       | :math   |
         | h:`t_1` | th:`M^a | IRF_{21 | \cdots` | \cdots` | math:`M | :`IRF_{ |
         |         | _{21}(\ | }(t_1)` |         |         | ^a_{2Ni | 2Ninteg |
         |         | infty)` |         |         |         | nteg}(\ | }(t_1)` |
         |         |         |         |         |         | infty)` |         |
         +---------+---------+---------+---------+---------+---------+---------+
         | :math:` | :math:` | :math:` | :math:` | :math:` | :math:` | :math:` |
         | \vdots` | \vdots` | \vdots` | \vdots` | \vdots` | \vdots` | \vdots` |
         +---------+---------+---------+---------+---------+---------+---------+
         | :mat    | :m      | :math:  | :math:` | :math:` | :math:` | :mat    |
         | h:`t_N` | ath:`M^ | `IRF_{N | \cdots` | \cdots` | M^a_{N_ | h:`IRF_ |
         |         | a_{N_{D | _{DoF}1 |         |         | {DoF}Ni | {N_{DoF |
         |         | oF}1}(\ | }(t_N)` |         |         | nteg}(\ | }Ninteg |
         |         | infty)` |         |         |         | infty)` | }(t_N)` |
         +---------+---------+---------+---------+---------+---------+---------+

.. container:: center

   .. container::
      :name: tab:IRFExcF

      .. table:: Output file format of ``IRF_excForce.tec``

         +----------------+----------------+----------------+----------------+
         | :math:`t_1`    | :math:         | :math:`\cdots` | :math:`IRF_    |
         |                | `IRF_{1}(t_1)` |                | {Ninteg}(t_1)` |
         +================+================+================+================+
         | :math:`t_2`    | :math:         | :math:`\cdots` | :math:`IRF_    |
         |                | `IRF_{1}(t_2)` |                | {Ninteg}(t_2)` |
         +----------------+----------------+----------------+----------------+
         | :math:`\vdots` | :math:`\vdots` | :math:`\vdots` | :math:`\vdots` |
         +----------------+----------------+----------------+----------------+
         | :math:`t_N`    | :math:         | :math:`\cdots` | :math:`IRF_    |
         |                | `IRF_{1}(t_N)` |                | {Ninteg}(t_N)` |
         +----------------+----------------+----------------+----------------+

``pressure.00XXX.dat``, ``kochin.00XXX.dat`` and ``freesurface.00XXX.dat`` are output files of pressure, Kochin and free surface, respectively, for a specific problem-XXX. The problem number is defined as in order of the diffraction problem (:math:`Nbeta`), the radiation problem (:math:`Ndof`) and for each frequency. So problem-001 is the, first frequency and first wave direction, diffraction problem. Suppose :math:`Nbeta=1`, then problem-002 is the first frequency radiation problem DoF 1. If :math:`Ndof=6` then problem-008 is the second frequency diffraction problem.

-  ``pressure.00XXX.dat`` is a pressure output file for the problem-XXX. In each file, the absolute value of pressure, :math:`|P|`, (Pa) and the phase, :math:`\angle P`, (rad) are given for each panel. The format of the output file is given in Table `11 <#tab:pressure>`__.

   .. container:: center

      .. container::
         :name: tab:pressure

         .. table:: Output file format of ``pressure.00XXX.dat``

            +-------------+-------------+-------------+-------------+-------------+
            | :math:`x_1` | :math:`y_1` | :math:`z_1` | :math:      | :math:`\ang |
            |             |             |             | `|P(\boldsy | le P(\bolds |
            |             |             |             | mbol x_1)|` | ymbol x_1)` |
            +=============+=============+=============+=============+=============+
            | ⋮           | ⋮           | ⋮           | ⋮           | ⋮           |
            +-------------+-------------+-------------+-------------+-------------+
            | :math:`     | :math:`     | :math:`     | :m          | :math:`     |
            | x_{Npanel}` | y_{Npanel}` | z_{Npanel}` | ath:`|P(\bo | \angle P(\b |
            |             |             |             | ldsymbol x_ | oldsymbol x |
            |             |             |             | {Npanel})|` | _{Npanel})` |
            +-------------+-------------+-------------+-------------+-------------+

-  ``kochin.00XXX.dat`` is an output file of the Kochin function on a prescribed direction for the problem-XXX. In each file, depending on the diffraction/radiation problem, the computed absolute value of the Kochin, :math:`|\mathcal{H}|`, and the phase, :math:`\angle \mathcal{H}`, (rad) are saved for each direction, :math:`\vartheta`. The format of the output file is given in Table `12 <#tab:kochin>`__.

   .. container:: center

      .. container::
         :name: tab:kochin

         .. table:: Output file format of *kochin.00XXX.dat*

            +----------------------+----------------------+----------------------+
            | :math:`\vartheta_1`  | :math:`|\mathc       | :math:`\angle \math  |
            |                      | al{H}(\vartheta_1)|` | cal{H}(\vartheta_1)` |
            +======================+======================+======================+
            | ⋮                    | ⋮                    | ⋮                    |
            +----------------------+----------------------+----------------------+
            | :math:`\va           | :math                | :math:`\an           |
            | rtheta_{N\vartheta}` | :`|\mathcal{H}(\vart | gle \mathcal{H}(\var |
            |                      | heta_{N\vartheta})|` | theta_{N\vartheta})` |
            +----------------------+----------------------+----------------------+

-  ``freesurface.00XXX.dat`` is an output file of the free-surface elevation on a prescribed free-surface domain for the problem-XXX. In each file, depending on the diffraction/radiation problem, the computed absolute value of the free-surface elevation, :math:`|\eta|`, and the phase, :math:`\angle \eta`, (rad) are saved for each free-surface panel position. The format of the output file is given in Table `13 <#tab:freesurface>`__.

   .. container:: center

      .. container::
         :name: tab:freesurface

         .. table:: Output file format of ``freesurface.00XXX.dat``

            +----------+----------+----------+----------+----------+----------+
            | :ma      | :ma      | :        | :math:   | :mat     | :mat     |
            | th:`x_1` | th:`y_1` | math:`|\ | `\angle  | h:`Re[ \ | h:`Im[ \ |
            |          |          | eta(\vec | \eta(\ve | eta(\vec | eta(\vec |
            |          |          | {x}_1)|` | c{x}_1)` | {x}_1)]` | {x}_1)]` |
            +==========+==========+==========+==========+==========+==========+
            | ⋮        | ⋮        | ⋮        | ⋮        | ⋮        | ⋮        |
            +----------+----------+----------+----------+----------+----------+
            | :m       | :m       | :math:`| | :math    | :ma      | :ma      |
            | ath:`x_{ | ath:`y_{ | \eta(\ve | :`\angle | th:`Re[  | th:`Im[  |
            | Npanel}` | Npanel}` | c{x}_{Np |  \eta(\v | \eta(\ve | \eta(\ve |
            |          |          | anel})|` | ec{x}_{N | c{x}_{Np | c{x}_{Np |
            |          |          |          | panel})` | anel})]` | anel})]` |
            +----------+----------+----------+----------+----------+----------+

``OUT_QTFM_N.dat`` and ``OUT_QTFP_N.dat`` are the output files of difference- and sum-frequencies QTF. The QTF results are either the total QTF or parts of the QTF terms that depend on the user choice QTF post-processing parameters in ``Nemoh.cal``. The QTF values are given in the absolute value with the phase in deg and real-imaginary parts. The QTF values are normalized by :math:`\rho g`. The ’frequency’ type, [rad/s, Hz, s], depends on the user choice in the ``Nemoh.cal``. The format of the output file is given in Table `14 <#tab:QTF>`__. Only the lower triangular part of the QTF matrix is saved in the file. The full difference-frequency QTF matrix can be constructed with the lower triangular part of the matrix and the upper triangular part which is in conjugate-symmetric with the lower part. The upper triangular part of the sum-frequency QTF is symmetric with the lower triangular part. A Matlab file for reading this output file is provided in ``matlabRoutines/`` and will be described in the next section.

.. container:: center

   .. container::
      :name: tab:QTF

      .. table:: Output file format of ``OUT_QTFM_N.dat`` and ``OUT_QTFP_N.dat``

         +-------+-------+-------+-------+-------+-------+-------+-------+-------+
         | :math | :math | :mat  | :mat  | :ma   | :math | :ma   | :m    | :m    |
         | :`f_{ | :`f_{ | h:`\b | h:`\b | th:`D | :`|QT | th:`\ | ath:` | ath:` |
         | 1_1}` | 2_1}` | eta_{ | eta_{ | oF_1` | F|/\r | angle | Re[QT | Im[QT |
         |       |       | 1_1}` | 2_1}` |       | ho g` |  QTF` | F]/\r | F]/\r |
         |       |       |       |       |       |       |       | ho g` | ho g` |
         +=======+=======+=======+=======+=======+=======+=======+=======+=======+
         | :math | :math | :mat  | :mat  | :ma   | :math | :ma   | :m    | :m    |
         | :`f_{ | :`f_{ | h:`\b | h:`\b | th:`D | :`|QT | th:`\ | ath:` | ath:` |
         | 1_2}` | 2_1}` | eta_{ | eta_{ | oF_1` | F|/\r | angle | Re[QT | Im[QT |
         |       |       | 1_1}` | 2_1}` |       | ho g` |  QTF` | F]/\r | F]/\r |
         |       |       |       |       |       |       |       | ho g` | ho g` |
         +-------+-------+-------+-------+-------+-------+-------+-------+-------+
         | :mat  | :mat  | :mat  | :mat  | :mat  | :mat  | :mat  | :mat  | :mat  |
         | h:`\v | h:`\v | h:`\v | h:`\v | h:`\v | h:`\v | h:`\v | h:`\v | h:`\v |
         | dots` | dots` | dots` | dots` | dots` | dots` | dots` | dots` | dots` |
         +-------+-------+-------+-------+-------+-------+-------+-------+-------+
         | :ma   | :ma   | :math | :math | :ma   | :math | :ma   | :m    | :m    |
         | th:`f | th:`f | :`\be | :`\be | th:`D | :`|QT | th:`\ | ath:` | ath:` |
         | _{1_{ | _{2_{ | ta_{1 | ta_{2 | oF_{N | F|/\r | angle | Re[QT | Im[QT |
         | Nf}}` | Nf}}` | _{Nbe | _{Nbe | Dof}` | ho g` |  QTF` | F]/\r | F]/\r |
         |       |       | ta}}` | ta}}` |       |       |       | ho g` | ho g` |
         +-------+-------+-------+-------+-------+-------+-------+-------+-------+

.. _`Sec:MatlabFiles`:

Supporting Matlab files
=======================

Following Matlab directories, containing a set of functions, are provided in ``matlabRoutines/``,

-  ``NemohWrapper``: This is for running NEMOH executables in MATLAB environment.

-  ``GMSHconverter``: There are two codes, first, for converting body mesh file output from GMSH to NEMOH, DIODORE and HYDROSTAR formats and second, for converting free-surface mesh file output from GMSH to NEMOH and HYDROSTAR formats.

-  ``postproc_testcases``: There are two main codes for plotting results from NEMOH and HYDROSTAR. First, for plotting hydrodynamic coefficients results and second for plotting QTF results. This code can be executed after all data in one specific test cases are obtained.

.. _`Sec:Testcase`:

Test-cases
==========

The following test cases are provided for verification with the original Aquaplus software (which is the ancestor of NEMOH) and/or HYDROSTAR commercial software :raw-latex:`\cite{HYDROSTAR}`. Note that Tecplot’s layout files ``.lay`` are provided in the relevant test case folder for plotting in Tecplot.

-  **1_Cylinder**: half-symmetric body mesh, deep water case, wave direction :math:`0^{\circ}`. The results are shown in :numref:`fig:Cylinder`.

.. _`fig:Cylinder`:
.. figure:: figures/Ver_Cylinder.svg

   Comparison of the first order results between NEMOH and AQUAPLUS for the test case **1_Cylinder**

-  **2_2Bodies**: half-symmetric body mesh, two different bodies, water depth :math:`20` m, wave direction :math:`45^{\circ}`. The results are shown in :numref:`fig:2Bodies`.

.. _`fig:2Bodies`:
.. figure:: figures/Ver_2Bodies.svg

   Comparison of the first order results between NEMOH and AQUAPLUS for the test case **2_2Bodies**

-  **3_Nonsymmetrical**: full non-symmetrical body mesh, deep-water, wave direction :math:`0^{\circ}`. Comparison of NEMOH results against Aquaplus are shown in :numref:`fig:NonSymmetrical_1` and :numref:`fig:NonSymmetrical_2`, a slight difference are observed in the results. Added mass and damping coefficients comparison between NEMOH and HYDROSTAR are shown in :numref:`fig:NonSymmetrical_mass` and :numref:`fig:NonSymmetrical_damp`, and for the excitation force is in :numref:`fig:NonSymmetrical_excforce`. Good agreement between NEMOH and HYDROSTAR is achieved.

.. _`fig:NonSymmetrical_1`:
.. figure:: figures/Ver_NonSymmetrical_1.svg

   Comparison of the first order results between NEMOH and AQUAPLUS for the test case **3_Nonsymmetrical**

.. _`fig:NonSymmetrical_2`:
.. figure:: figures/Ver_NonSymmetrical_2.svg

   Comparison of the first order results between NEMOH and AQUAPLUS for the test case **3_Nonsymmetrical**

.. _`fig:NonSymmetrical_mass`:
.. figure:: figures/Ver_NonSymmetrical_addedmass.svg

   Comparison of added mass coefficients between NEMOH, red dashed-line, and HYDROSTAR, blue solid-line, for the test case **3_Nonsymmetrical**

.. _`fig:NonSymmetrical_damp`:
.. figure:: figures/Ver_NonSymmetrical_dampcoef.svg

   Comparison of damping coefficients between NEMOH, red dashed-line, and HYDROSTAR, blue solid-line, for the test case **3_Nonsymmetrical**

.. _`fig:NonSymmetrical_excforce`:
.. figure:: figures/Ver_NonSymmetrical_excitationforce.svg

   Comparison of excitation force between NEMOH, red dashed-line, and HYDROSTAR, blue solid-line, for the test case **3_Nonsymmetrical**

-  **4_Postprocessing**: half-symmetric body mesh, water depth :math:`10` m, wave direction :math:`0^{\circ}`. This test case shows a comparison of the free-surface elevation and the Kochin function. The results are shown in :numref:`fig:PostProcessing`. The phase difference, :math:`\pm \pi/2`, of wave elevation between NEMOH and AQUAPLUS is due to different conventions of the incident potential.

.. _`fig:PostProcessing`:
.. figure:: figures/Ver_PostProcessing.svg

   Comparison of the diffracted wave elevation, the diffraction Kochin function between NEMOH and AQUAPLUS, test case **4_Postprocessing**

-  **5_QuickTest** shows a quantitative comparison of force and free-surface for the first-frequency diffraction problem. The comparison results are shown in the command window for all the test cases inside the directory ``5_QuickTest``.

-  **6_box_coarsemesh** is showing the procedure for running the code starting with the executable **``mesh``** with a coarse description mesh file, ``meshbox``. No reference data is given in this test case.

-  **7_Solvers_Check_OC3** is testing the performance of the three difference linear solvers, Gauss elimination, LU decomposition and GMRES. Reference logfiles reporting the computational time of the solvers are provided.

-  **8a_Cylinder_irregfreq** shows the results with and without irregular frequencies removal (IRR) method. The results are verified against HYDROSTAR with IRR and shown in :numref:`fig:Cylinder_IRR_addedmass` and :numref:`fig:Cylinder_IRR_dampcoef` for the added mass and damping coefficients and in :numref:`fig:Cylinder_IRR_excforce` for the excitation forces. The mesh used was obtained using GMSH :raw-latex:`\cite{gmsh}` and is shown in :numref:`fig:meshesCylinder`.

.. _`fig:meshesCylinder`:
.. figure:: figures/Cylinder/mesh.svg

   Body boundary mesh for the Cylinder used for test case **8a_Cylinder_irregfreq** and **8b_QTF_Cylinder**.

.. _`fig:Cylinder_IRR_addedmass`:
.. figure:: figures/Cylinder/addedmass.svg

   Comparison of added masscoefficients between NEMOH without irregular frequencies removal (IRR), green dash-dotted line, NEMOH with IRR, red dashed-line and HYDROSTAR with IRR, blue solid-line, for the test-case **8a_Cylinder_irregfreq**

.. _`fig:Cylinder_IRR_dampcoef`:
.. figure:: figures/Cylinder/dampcoef.svg

   Comparison of damping coefficients between NEMOH without irregular frequencies removal (IRR), green dash-dotted line, NEMOH with IRR, red dashed-line and HYDROSTAR with IRR, blue solid-line, for the test-case **8a_Cylinder_irregfreq**

.. _`fig:Cylinder_IRR_excforce`:
.. figure:: figures/Cylinder/excForce.svg

   Comparison of excitation force between NEMOH without irregular frequencies removal (IRR), green dash-dotted line, NEMOH with IRR, red dashed-line and HYDROSTAR with IRR, blue solid-line, for the test-case 8a_Cylinder_irregfreq

The following test cases are provided for the QTF verification with HYDROSTAR software :raw-latex:`\cite{HYDROSTAR}`.

-  **8b_QTF_Cylinder**: full body mesh with lid panels, CoG :math:`(0,0,0)`, deep water, wave direction :math:`0^{\circ}`, the difference-frequency QTF DUOK+HASBO. The results are shown in the density plot, :numref:`fig:QTFM_Cylinder_surge`, :numref:`fig:QTFM_Cylinder_heave` and :numref:`fig:QTFM_Cylinder_pitch`, and in the off-diagonal line plot, :numref:`fig:QTFM_diag_Cylinder_surge`, :numref:`fig:QTFM_diag_Cylinder_heave` and :numref:`fig:QTFM_diag_Cylinder_pitch`. The mesh used was obtained using GMSH :raw-latex:`\cite{gmsh}` and is shown in :numref:`fig:meshesCylinder`.

.. _`fig:QTFM_Cylinder_surge`:
.. figure:: figures/Cylinder/QTFsurge.svg

   Density plots of the normalized surge difference frequency QTF magnitude (without the free-surface integrals) for the floating Cylinder (test case **8b_QTF_Cylinder**. HYDROSTAR results are on the left column, NEMOH results are on the middle column and the difference on the right column.

.. _`fig:QTFM_Cylinder_heave`:
.. figure:: figures/Cylinder/QTFheave.svg

   Density plots of the normalized heave difference frequency QTF magnitude (without the free-surface integrals) for the floating Cylinder (test case **8b_QTF_Cylinder**. HYDROSTAR results are on the left column, NEMOH results are on the middle column and the difference on the right column.

.. _`fig:QTFM_Cylinder_pitch`:
.. figure:: figures/Cylinder/QTFpitch.svg

   Density plots of the normalized pitch difference frequency QTF magnitude (without the free-surface integrals) for the floating Cylinder (test case **8b_QTF_Cylinder**. HYDROSTAR results are on the left column, NEMOH results are on the middle column and the difference on the right column.

.. _`fig:QTFM_diag_Cylinder_surge`:
.. figure:: figures/Cylinder/QTFsurge_diag.svg

   Comparison of the surge off-diagonal difference frequency QTF for the Cylinder (test case **8b_QTF_Cylinder**) between HYDROSTAR, real part (blue, solid-line), imaginary part (blue, dashed-dot line) and NEMOH, real part (red, dashed-line), imaginary part (red, dotted-line).

.. _`fig:QTFM_diag_Cylinder_heave`:
.. figure:: figures/Cylinder/QTFheave_diag.svg

   Comparison of the heave off-diagonal difference frequency QTF for the Cylinder (test case **8b_QTF_Cylinder**) between HYDROSTAR, real part (blue, solid-line), imaginary part (blue, dashed-dot line) and NEMOH, real part (red, dashed-line), imaginary part (red, dotted-line).

.. _`fig:QTFM_diag_Cylinder_pitch`:
.. figure:: figures/Cylinder/QTFpitch_diag.svg

   Comparison of the pitch off-diagonal difference frequency QTF for the Cylinder (test case **8b_QTF_Cylinder**) between HYDROSTAR, real part (blue, solid-line), imaginary part (blue, dashed-dot line) and NEMOH, real part (red, dashed-line), imaginary part (red, dotted-line).

-  **9_QTF_OC4_Semisubmersible**: full body mesh with lid panels, CoG :math:`(0,0,0)`, water depth 200 m, wave direction :math:`0^{\circ}` and :math:`30^{\circ}`, bi-directional QTF, the difference-frequency QTF DUOK+HASBO. The results are shown in the density plot, :numref:`fig:QTFM_OC4_surge`, :numref:`fig:QTFM_OC4_heave` and :numref:`fig:QTFM_OC4_pitch`, and in the off-diagonal line plot, :numref:`fig:QTFM_diag_OC4_surge`, :numref:`fig:QTFM_diag_OC4_heave` and :numref:`fig:QTFM_diag_OC4_pitch`, of the bi-directional QTF :math:`(\beta_1,\beta_2)=(0^{\circ},30^{\circ})`. The mesh used was obtained using GMSH :raw-latex:`\cite{gmsh}` and is shown in :numref:`fig:meshesOC4`.

.. _`fig:meshesOC4`:
.. figure:: figures/OC4/bodymesh.svg

   Body boundary mesh for for the OC4-platform used for test case **9_QTF_OC4_Semisubmersible**.

.. _`fig:QTFM_OC4_surge`:
.. figure:: figures/OC4/QTFM_Surge_beta030.svg

   Density plots of the normalized bi-directional, :math:`(\beta_1,\beta_2)=(0^{\circ},30^{\circ})`, surge difference frequency QTF magnitude (without the free-surface integrals) for the floating OC4-semisubmersible platform (test case **9_QTF_OC4_Semisubmersible**). HYDROSTAR results are on the left column, NEMOH results are on the middle column and the difference on the right column.

.. _`fig:QTFM_OC4_heave`:
.. figure:: figures/OC4/QTFM_Heave_beta030.svg

   Density plots of the normalized bi-directional, :math:`(\beta_1,\beta_2)=(0^{\circ},30^{\circ})`, heave difference frequency QTF magnitude (without the free-surface integrals) for the floating OC4-semisubmersible platform (test case **9_QTF_OC4_Semisubmersible**). HYDROSTAR results are on the left column, NEMOH results are on the middle column and the difference on the right column.

.. _`fig:QTFM_OC4_pitch`:
.. figure:: figures/OC4/QTFM_Pitch_beta030.svg

   Density plots of the normalized bi-directional, :math:`(\beta_1,\beta_2)=(0^{\circ},30^{\circ})`, pitch difference frequency QTF magnitude (without the free-surface integrals) for the floating OC4-semisubmersible platform (test case **9_QTF_OC4_Semisubmersible**). HYDROSTAR results are on the left column, NEMOH results are on the middle column and the difference on the right column.

.. _`fig:QTFM_diag_OC4_surge`:
.. figure:: figures/OC4/QTFM_Surge_beta030_diag.svg

   Comparison of the off-diagonal bi-directional, :math:`(\beta_1,\beta_2)=(0^{\circ},30^{\circ})`, surge difference frequency QTF for the OC4-semisubmersible platform (test case **9_QTF_OC4_Semisubmersible**) between HYDROSTAR, real part (blue, solid-line), imaginary part (blue, dashed-dot line) and NEMOH, real part (red, dashed-line), imaginary part (red, dotted-line).

.. _`fig:QTFM_diag_OC4_heave`:
.. figure:: figures/OC4/QTFM_Heave_beta030_diag.svg

   Comparison of the off-diagonal bi-directional, :math:`(\beta_1,\beta_2)=(0^{\circ},30^{\circ})`, heave difference frequency QTF for the OC4-semisubmersible platform (test case **9_QTF_OC4_Semisubmersible**) between HYDROSTAR, real part (blue, solid-line), imaginary part (blue, dashed-dot line) and NEMOH, real part (red, dashed-line), imaginary part (red, dotted-line).

.. _`fig:QTFM_diag_OC4_pitch`:
.. figure:: figures/OC4/QTFM_Pitch_beta030_diag.svg

   Comparison of the off-diagonal bi-directional, :math:`(\beta_1,\beta_2)=(0^{\circ},30^{\circ})`, pitch difference frequency QTF for the OC4-semisubmersible platform (test case **9_QTF_OC4_Semisubmersible**) between HYDROSTAR, real part (blue, solid-line), imaginary part (blue, dashed-dot line) and NEMOH, real part (red, dashed-line), imaginary part (red, dotted-line).

-  **10a_QTF_SOFTWIND**: half symmetric body mesh with lid panels, CoG :math:`(0,0,-71.56)`, water depth 200 m, wave direction :math:`0^{\circ}` and :math:`30^{\circ}`, bi-directional QTF, the difference-frequency QTF DUOK+HASBO. The results are shown in the density plot, :numref:`fig:QTFM_SOFTWIND_surge`, :numref:`fig:QTFM_SOFTWIND_heave` and :numref:`fig:QTFM_SOFTWIND_pitch`, and in the off-diagonal line plot, :numref:`fig:QTFM_diag_softwind_surge`, :numref:`fig:QTFM_diag_softwind_heave` and :numref:`fig:QTFM_diag_softwind_pitch`, of the bi-directional QTF :math:`(\beta_1,\beta_2)=(0^{\circ},30^{\circ})`. The mesh used was obtained using GMSH :raw-latex:`\cite{gmsh}` and is shown in :numref:`fig:meshesSoftwind_body`.

.. _`fig:meshesSoftwind_body`:
.. figure:: figures/Softwind/bodymesh.svg

   Body boundary mesh for the SOFTWIND platform, used in test cases **10a_QTF_SOFTWIND** and  **10b_QTF_SOFTWIND_FS**

.. _`fig:meshesSoftwind_FS`:
.. figure:: figures/Softwind/FSmesh.svg

   Free surface mesh for the SOFTWIND platform, used in test case **10b_QTF_SOFTWIND_FS**

.. _`fig:QTFM_SOFTWIND_surge`:
.. figure:: figures/Softwind/QTFM_Surge_beta030.svg

   Density plots of the normalized bi-directional, :math:`(\beta_1,\beta_2)=(0^{\circ},30^{\circ})`, surge difference frequency QTF magnitude (without the free-surface integrals) for the floating SOFTWIND platform (test case **10a_QTF_SOFTWIND**). HYDROSTAR results are on the left column, NEMOH results are on the middle column and the difference on the right column.

.. _`fig:QTFM_SOFTWIND_heave`:
.. figure:: figures/Softwind/QTFM_Heave_beta030.svg

   Density plots of the normalized bi-directional, :math:`(\beta_1,\beta_2)=(0^{\circ},30^{\circ})`, heave difference frequency QTF magnitude (without the free-surface integrals) for the floating SOFTWIND platform (test case **10a_QTF_SOFTWIND**). HYDROSTAR results are on the left column, NEMOH results are on the middle column and the difference on the right column.

.. _`fig:QTFM_SOFTWIND_pitch`:
.. figure:: figures/Softwind/QTFM_Pitch_beta030.svg

   Density plots of the normalized bi-directional, :math:`(\beta_1,\beta_2)=(0^{\circ},30^{\circ})`, pitch difference frequency QTF magnitude (without the free-surface integrals) for the floating SOFTWIND platform (test case **10a_QTF_SOFTWIND**). HYDROSTAR results are on the left column, NEMOH results are on the middle column and the difference on the right column.

.. _`fig:QTFM_diag_softwind_surge`:
.. figure:: figures/Softwind/QTFM_Surge_beta030_diag.svg

   Comparison of the off-diagonal bi-directional, :math:`(\beta_1,\beta_2)=(0^{\circ},30^{\circ})`, surge difference frequency QTF for the SOFTWIND platform (test case **10a_QTF_SOFTWIND**) between HYDROSTAR, real part (blue, solid-line), imaginary part (blue, dashed-dot line) and NEMOH, real part (red, dashed-line), imaginary part (red, dotted-line).

.. _`fig:QTFM_diag_softwind_heave`:
.. figure:: figures/Softwind/QTFM_Heave_beta030_diag.svg

   Comparison of the off-diagonal bi-directional, :math:`(\beta_1,\beta_2)=(0^{\circ},30^{\circ})`, heave difference frequency QTF for the SOFTWIND platform (test case **10a_QTF_SOFTWIND**) between HYDROSTAR, real part (blue, solid-line), imaginary part (blue, dashed-dot line) and NEMOH, real part (red, dashed-line), imaginary part (red, dotted-line).

.. _`fig:QTFM_diag_softwind_pitch`:
.. figure:: figures/Softwind/QTFM_Pitch_beta030_diag.svg

   Comparison of the off-diagonal bi-directional, :math:`(\beta_1,\beta_2)=(0^{\circ},30^{\circ})`, pitch difference frequency QTF for the SOFTWIND platform (test case **10a_QTF_SOFTWIND**) between HYDROSTAR, real part (blue, solid-line), imaginary part (blue, dashed-dot line) and NEMOH, real part (red, dashed-line), imaginary part (red, dotted-line).

-  **10b_QTF_SOFTWIND_FS**: half symmetric body mesh without lid panels, half symmetric free-surface mesh, CoG :math:`(0,0,-71.56)`, water depth 200 m, wave direction :math:`0^{\circ}`, the sum-frequency total QTF DUOK+HASBO+HASFS+ASYMP. The results are shown in the density plot, :numref:`fig:QTFP_SOFTWIND_surge`, :numref:`fig:QTFP_SOFTWIND_heave` and :numref:`fig:QTFP_SOFTWIND_pitch` and in the off-diagonal line plot, :numref:`fig:QTFP_SOFTWIND_DIAG_surge`, :numref:`fig:QTFP_SOFTWIND_DIAG_heave` and :numref:`fig:QTFP_SOFTWIND_DIAG_pitch`. The mesh used was obtained using GMSH :raw-latex:`\cite{gmsh}` and is shown in :numref:`fig:meshesSoftwind_body` and :numref:`fig:meshesSoftwind_FS`.

.. _`fig:QTFP_SOFTWIND_surge`:
.. figure:: figures/Softwind/QTFP_Surge_beta00.svg

   Density plots of the normalized surge sum-frequency full QTF magnitude (including the free-surface integrals) for the floating SOFTWIND platform (test case **10b_QTF_SOFTWIND_FS**). HYDROSTAR results are on the left column, NEMOH results are on the middle column and the difference in the right column.

.. _`fig:QTFP_SOFTWIND_heave`:
.. figure:: figures/Softwind/QTFP_Heave_beta00.svg

   Density plots of the normalized heave sum-frequency full QTF magnitude (including the free-surface integrals) for the floating SOFTWIND platform (test case **10b_QTF_SOFTWIND_FS**). HYDROSTAR results are on the left column, NEMOH results are on the middle column and the difference in the right column.

.. _`fig:QTFP_SOFTWIND_pitch`:
.. figure:: figures/Softwind/QTFP_Pitch_beta00.svg

   Density plots of the normalized pitch sum-frequency full QTF magnitude (including the free-surface integrals) for the floating SOFTWIND platform (test case **10b_QTF_SOFTWIND_FS**). HYDROSTAR results are on the left column, NEMOH results are on the middle column and the difference in the right column.

.. _`fig:QTFP_SOFTWIND_DIAG_surge`:
.. figure:: figures/Softwind/QTFP_Surge_beta00_diag.svg

   Comparison of the off-diagonal surge sum-frequency full QTF for SOFTWIND platform (test case **10b_QTF_SOFTWIND_FS**) between HYDROSTAR, real part (blue, solid-line), imaginary part (blue, dashed-dot line) and NEMOH, real part (red, dashed-line), imaginary part (red, dotted-line).

.. _`fig:QTFP_SOFTWIND_DIAG_heave`:
.. figure:: figures/Softwind/QTFP_Heave_beta00_diag.svg

   Comparison of the off-diagonal heave sum-frequency full QTF for SOFTWIND platform (test case **10b_QTF_SOFTWIND_FS**) between HYDROSTAR, real part (blue, solid-line), imaginary part (blue, dashed-dot line) and NEMOH, real part (red, dashed-line), imaginary part (red, dotted-line).

.. _`fig:QTFP_SOFTWIND_DIAG_pitch`:
.. figure:: figures/Softwind/QTFP_Pitch_beta00_diag.svg

   Comparison of the off-diagonal pitch sum-frequency full QTF for SOFTWIND platform (test case **10b_QTF_SOFTWIND_FS**) between HYDROSTAR, real part (blue, solid-line), imaginary part (blue, dashed-dot line) and NEMOH, real part (red, dashed-line), imaginary part (red, dotted-line).

-  **11_QTF_OC3_Hywind**: full body mesh with lid panels, CoG :math:`(0,0,0)`, water depth 320 m, wave direction :math:`0^{\circ}`, NEMOH1 uses GMRES solver, the difference-frequency QTF DUOK+HASBO. The results are shown in the density plot, :numref:`fig:QTFM_OC3_HYWIND_surge`, :numref:`fig:QTFM_OC3_HYWIND_heave` and :numref:`fig:QTFM_OC3_HYWIND_pitch`, and in the off-diagonal line plot, :numref:`fig:QTFM_diag_OC3_HYWIND_surge`, :numref:`fig:QTFM_diag_OC3_HYWIND_heave` and :numref:`fig:QTFM_diag_OC3_HYWIND_pitch`, of the difference-frequency QTF. The mesh used was obtained using GMSH :raw-latex:`\cite{gmsh}` and is shown in :numref:`fig:meshesHYWIND`.

.. _`fig:meshesHYWIND`:
.. figure:: figures/OC3_HYWIND/bodyMesh.svg

   Body boundary mesh for OC3-HYWIND platform, test case **11_QTF_OC3_Hywind**.

.. _`fig:QTFM_OC3_HYWIND_surge`:
.. figure:: figures/OC3_HYWIND/QTFM_Surge.svg

   Density plots of the normalized surge difference frequency QTF magnitude (without the free-surface integrals) for the floating OC3-HYWIND platform (test case **11_QTF_OC3_Hywind**). HYDROSTAR results are on the left column, NEMOH results are on the middle column and the difference on the right column.

.. _`fig:QTFM_OC3_HYWIND_heave`:
.. figure:: figures/OC3_HYWIND/QTFM_Heave.svg

   Density plots of the normalized heave difference frequency QTF magnitude (without the free-surface integrals) for the floating OC3-HYWIND platform (test case **11_QTF_OC3_Hywind**). HYDROSTAR results are on the left column, NEMOH results are on the middle column and the difference on the right column.

.. _`fig:QTFM_OC3_HYWIND_pitch`:
.. figure:: figures/OC3_HYWIND/QTFM_Pitch.svg

   Density plots of the normalized pitch difference frequency QTF magnitude (without the free-surface integrals) for the floating OC3-HYWIND platform (test case **11_QTF_OC3_Hywind**). HYDROSTAR results are on the left column, NEMOH results are on the middle column and the difference on the right column.

.. _`fig:QTFM_diag_OC3_HYWIND_surge`:
.. figure:: figures/OC3_HYWIND/QTFM_Surge_diag.svg

   Comparison of the off-diagonal surge difference frequency QTF for the OC3-HYWIND platform (test case **11_QTF_OC3_Hywind**) between HYDROSTAR, real part (blue, solid-line), imaginary part (blue, dashed-dot line) and NEMOH, real part (red, dashed-line), imaginary part (red, dotted-line).

.. _`fig:QTFM_diag_OC3_HYWIND_heave`:
.. figure:: figures/OC3_HYWIND/QTFM_Heave_diag.svg

   Comparison of the off-diagonal heave difference frequency QTF for the OC3-HYWIND platform (test case **11_QTF_OC3_Hywind**) between HYDROSTAR, real part (blue, solid-line), imaginary part (blue, dashed-dot line) and NEMOH, real part (red, dashed-line), imaginary part (red, dotted-line).

.. _`fig:QTFM_diag_OC3_HYWIND_pitch`:
.. figure:: figures/OC3_HYWIND/QTFM_Pitch_diag.svg

   Comparison of the off-diagonal pitch difference frequency QTF for the OC3-HYWIND platform (test case **11_QTF_OC3_Hywind**) between HYDROSTAR, real part (blue, solid-line), imaginary part (blue, dashed-dot line) and NEMOH, real part (red, dashed-line), imaginary part (red, dotted-line).

Full description of the QTF test-cases results is reported in :raw-latex:`\citep{Kurnia22_JH,Kurnia22}`. Note that the QTF comparisons between NEMOH and HYDROSTAR for the bidirectional case are in good agreement only if the direction is switched, in NEMOH :math:`\beta=(\beta_1,\beta_2)` and in Hydrostar :math:`\beta=(\beta_2,\beta_1)`; further investigation regarding this is needed. The imaginary part of QTFs have also a difference sign between NEMOH and HYDROSTAR that may be due to different conventions of the incident potential.

Acknowledgement
===============

This work was done within the framework of the FLOATECH project. This project has received funding from the European Union’s Horizon 2020 research and innovation program under grant agreement No 101007142.

We thank Moran Charlou from LHEEA, ECN for his help in finalizing NEMOH v3.0 release in the gitlab.
