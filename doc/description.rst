
###########
Description
###########

*******************************
Modelling and Numerical aspects
*******************************

This section provides background information on NEMOH v3.0 with a focus on basic scientific ideas.

NEMOH v3.0 contains two main modules. First, NEMOH1 solves linear diffraction and radiation problems of wave-structure interaction using a 3-D boundary element method in the frequency domain. Second, NEMOH2, an extended module for computing difference- and sum- frequencies Quadratic Transfer Functions (QTFs) for fixed or floating structures.
The following subsections describe the underlying modelling and numerical approaches used in NEMOH.

Notations
=========

Symbols
^^^^^^^

.. table:: List of Symbols
   :name: tab:symbols

   ================= ================== ====================================================================================
   Symbol            Unit               Description
   =========================================================================================================================
   :math:`D`         [m]                Water depth
   :math:`T`         [s]                Wave period
   :math:`\omega`    [rad/s]            Wave radial frequency, :math:`\omega = 2\pi/T`
   :math:`\lambda`   [m]                Wave length
   :math:`k`         [rad/m]            Wave number, :math:`k = 2\pi/\lambda`
   :math:`\beta`     [rad]              Wave propagation direction
   :math:`t`         [s]                Time
   :math:`p`         [Pa]               Pressure
   :math:`\phi`      [m:sup:`2`/s]      Velocity potential
   :math:`\vec{V}`   [m/s]              Velocity, :math:`\vec{V} = \vec{\nabla}\phi`
   :math:`g`         [m.s:sup:`-2`]     Gravitational acceleration
   :math:`\rho`      [kg/m:sup:`3`]     Water density
   :math:`a`         [m]                Wave amplitude
   :math:`\eta`      [m]                Free surface elevation
   ================= ================== ====================================================================================


Frames
^^^^^^

As sketched in :numref:`fig:sketch`, we consider fluid domain in the Cartesian coordinate :math:`\boldsymbol x=(\mathbf{x},z)` with :math:`\mathbf{x}=(x,y)` the horizontal coordinates perpendicular to the :math:`z` axis in the opposite direction of gravity :math:`\boldsymbol g`. Free-surface boundary :math:`S_F` is defined by the free surface elevation at time :math:`t`, denoted as :math:`\eta(\mathbf{x},t)` with respect to the mean water level at :math:`z=0`. The fluid velocity potential is denoted as :math:`\Phi(\boldsymbol x,t)` with :math:`\boldsymbol x` in fluid domain :math:`V_{\Omega}`.

.. figure:: figures/Sketch.png
   :align: center
   :name: fig:sketch

   Sketch definition of the system

The floating body has 6 degrees of freedom (DoF), :math:`\boldsymbol\xi=(\boldsymbol{X},\boldsymbol{\theta})` where the positions, :math:`\boldsymbol{X}=(X,Y,Z)` and the orientations, :math:`\boldsymbol{\theta}=(\theta_1,\theta_2,\theta_3)` are determined at the center of gravity (COG). Displacements of points at the hull are specified by a body of vector :math:`\boldsymbol r` with respect to the COG as :math:`\boldsymbol{\mathcal{X}}=\boldsymbol{X}+R(\boldsymbol{r})`. :math:`R` is a rotation operator where :math:`R(\boldsymbol r)\approx \boldsymbol\theta \times \boldsymbol r`. The velocity of the points at the hull is expressed as :math:`\dot{\boldsymbol{\mathcal{X}}}`.

On the body hull :math:`S_B`, the wetted part is defined as a function :math:`z=\zeta(\mathbf{x},t)`. The normalized normal vector is defined as directed toward the fluid domain, :math:`\boldsymbol n=-\boldsymbol N/|\boldsymbol N|` with :math:`\boldsymbol N=\left(-\nabla_2\zeta,1 \right)` where :math:`\nabla_2` is the two-dimensional gradient in :math:`\mathbf{x}`. Then the six-dimensional generalized normal vector is defined as :math:`\boldsymbol\nu=(\boldsymbol n,\boldsymbol r \times \boldsymbol n)^T`, with :math:`( )^T` the matrix transpose operator.


Expression of physical quantities
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All time-varying physical quantities are represented in the frequency domain as harmonic complex numbers for the first-order part, proportional to the wave amplitude :math:`a`, e.g. quantity :math:`X` is written:

.. math::

   X(\mathbf{x},t) = a X(\mathbf{x}, \omega) \mathrm{e}^{-i \omega t}

:math:`X(\mathbf{x}, \omega)` is conveniently written in the exponential form:

.. math::

   X(\mathbf{x}, \omega) = |X(\mathbf{x}, \omega)| \mathrm{e}^{i \angle X(\mathbf{x}, \omega)}

and its time-domain equivalent :math:`\underline{X}` is its real part:

.. math::

   \underline{X}(\mathbf{x},t) = \Re\left\lbrace X(\mathbf{x},t)\right\rbrace = \Re\left\lbrace a |X(\mathbf{x}, \omega)| \mathrm{e}^{-i(\omega t - \angle X(\mathbf{x}, \omega))} \right\rbrace = a |X(\mathbf{x}, \omega)| \cos(\omega t - \angle X(\mathbf{x}, \omega))

.. note::

   :math:`X` may depend on spatial coordinates :math:`\mathbf{x}` or not, respectfully if it represents a field (free surface elevation, pressure, velocity potential) or a scalar (forces, motions).


Modelling
=========

NEMOH1, the first-order solver, is based on the following modelling principles:

-  Potential flow theory for an inviscid and incompressible fluid in irrotational motion. The set of equations solved looks for unknowns satisfying free-surface conditions, impermeable bottom condition, diffraction and radiation conditions on the body hull and radiation wave condition in the far field.

-  The harmonic fluid potential is defined as

   .. math::
      :label: Eq:PhiHarm

      \begin{aligned}
      \Phi(\mathbf{x},t)=\Re\left\lbrace a \Phi^{(1)}(\mathbf{x}, \omega)\mathrm{e}^{-i\omega t}\right\rbrace.
      \end{aligned}

   The total potential, :math:`\Phi`, is the sum of the incident potential, the diffraction potential and the radiation potential.

-  The incident wave is modelled using Airy's theory. In this framework, the incident wave elevation is written:

   .. math::
      :label: eq:eta_I

      \eta_I(\mathbf{x}, \omega) = \mathrm{e}^{\mathbf{k}\cdot\mathbf{x}}

   where :math:`a` is the wave amplitude and :math:`\beta` is the wave direction.
   :math:`\mathbf{k}=k(\cos(\beta),\sin(\beta),0)` is the wave number vector, related to the radial frequency :math:`\omega` by the dispersion relation:

   .. math::
      :label: eq:dispersion

      \omega^2 = g k \tanh(k D) \xrightarrow[D \longrightarrow \infty]{} g k

   The incident potential is defined as:

   .. math::
      :label: Eq:PhiI

      \Phi_I^{(1)}(\mathbf{x}, \omega)= - i \frac{g}{\omega} f_0(z) \mathrm{e}^{i \mathbf{k} \cdot \mathbf{x}}

   where:

   .. math::
      :label: eq:f_0

      f_0(z) = \frac{\cosh(k(D + z))}{\cosh(k D)} \xrightarrow[D \longrightarrow \infty]{} \mathrm{e}^{k z}

   The resulting incident pressure is:

   .. math::
      :label: eq:pressure

      p_I(\mathbf{x}, \omega) = \rho g f_0(z) \mathrm{e}^{\mathbf{k}\cdot\mathbf{x}}

-  The radiation potential is defined as :math:`\Phi_R(\boldsymbol x,t)=Re\left\lbrace \dot{\boldsymbol\xi}^{(1)}(t) \cdot \boldsymbol\psi(x)\right\rbrace` where :math:`\boldsymbol\psi(\boldsymbol x)` is the normalized vector radiation potential.

-  The three-dimensional linear potential flow problem around arbitrary body condition is reformulated in the Boundary Integral Equation (BIE) and transformed into the two-dimensional problem of the source distribution, :math:`\sigma`, on the body surface, :math:`S_B`, using Green’s second identity and the appropriate Green function, :math:`G(\boldsymbol x,\boldsymbol x')`.

-  The Green function is based on Delhommeau’s formulation and is available for finite and infinite water-depth, see :cite:t:`Delhommeau`.

-  The source distribution depends on the considered boundary condition problem. For each frequency and wave direction, the diffraction source distribution, :math:`\sigma_D(\boldsymbol x)`, depends on the position of the panels while the radiation source distribution, :math:`\sigma_{R_j}(\boldsymbol x)`, depends on the position of the panels and the considered degree of freedom :math:`j`.

-  Then, the BIE for :math:`\boldsymbol x \in S_B`, is expressed as, with flow points :math:`\boldsymbol x` and source points :math:`\boldsymbol x'`,

   .. math::
      :label: Eq:BIE_source_distribution

      \begin{aligned}
      \frac{1}{2}\sigma_{D,R_j}(\boldsymbol x)-\frac{1}{4\pi}\int_{S_B} \partial_n G(\boldsymbol x, \boldsymbol x') \sigma_{D,R_j}(\boldsymbol x') dS'=\mathcal{N}_{D,R_j}(\boldsymbol x).
      \end{aligned}

   where :math:`\mathcal{N}(\boldsymbol x)` is the body normal condition. The diffraction normal condition is defined as :math:`\mathcal{N}_D (\boldsymbol x)=-\partial_{n} \Phi_I^{(1)}(\boldsymbol x)`, the normalized radiation condition, :math:`\mathcal{N}_R (\boldsymbol x)=\partial_{n} \Phi_{R_j}(\boldsymbol x)`, with :math:`\Phi_{R_j}(\boldsymbol x)` is the vector component-:math:`j` of the normalized radiation potential :math:`\boldsymbol\psi(\boldsymbol x)`, explicitly :math:`\boldsymbol\psi=(\Phi_{R_1},\Phi_{R_2},\cdots,\Phi_{R_{Ndof}})`.

-  The diffraction potential, :math:`\Phi^{(1)}_{D}`, the normalized radiation potential vector component-:math:`j`, :math:`\Phi_{R_j}` and the corresponding velocities are then computed as follows, for the flow points in the fluid domain :math:`\boldsymbol x \in S_B \cup V_{\Omega_F}`,

   .. math::
      :label: Eq:BIE_Sol_Pot_Sb

      \begin{aligned}
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
      :label: Eq:RAO

      \begin{aligned}
      \left[-[\boldsymbol M+\boldsymbol M^a(\omega)]\omega^2-i\omega[\boldsymbol B(\omega)+\boldsymbol B_{add}]+[\boldsymbol K_h+\boldsymbol K_M]\right]\mathcal{\boldsymbol\xi}(\omega)=\boldsymbol F_{exc}(\omega)
      \end{aligned}

   where :math:`[\boldsymbol B_{add}]` and :math:`[\boldsymbol K_M]` are user-specified additional damping and stiffness matrices.


NEMOH2, the second-order QTF module, is based on the following principles

-  The second-order loads are composed of the quadratic part and the potential part, the detailed formulation is given in :cite:t:`Kurnia22_JH,Kurnia22`.

-  The quadratic part is based on the near-field method :cite:p:`CHEN88`.

-  The potential part is based on the
   indirect method :cite:p:`CHEN88,MOLIN79`.

Numerical Methods
=================

NEMOH1 uses the following numerical approach:

-  The BIE, Eq. :eq:`Eq:BIE_source_distribution`, is discretised using the constant panel method with quadrilateral mesh. This leads to a linear system with the influence coefficients matrix. The mesh is user-specified with the normal direction towards fluid.

-  Numerical implementation of the Green function is described in :cite:t:`Babarit15`.

-  Free-surface Green function integrands are pre-calculated with the discretized :math:`\omega^2r/g\in [0,100]` with 676 points in a constant scale and :math:`\omega^2(z+z')/g \in [-251,-1.6\, 10^{-6}]` with 130 points in logarithmic scale. A polynomial surface interpolation with the :math:`5^{th}` order Lagrange formula is used for interpolating any values in the specified interval.

-  The specified points for the interpolation of the Green function are finer than in the previous release. However, an option to switch the two different tabulated Green function data is available in the source file ``/Solver/Core/INITIALIZE_GREEN.f90`` with the parameter FLAG_IGREEN=1 or 2, 2 being the default.

-  Influence coefficients, the integration of :math:`\partial_n G(\boldsymbol x, \boldsymbol x')` over a body panel, is computed using Gauss-quadrature integration with a user-input number of Gauss-quadrature points.

-  The source distributions on body panels are then obtained after solving the corresponding linear system.

-  The linear system is solved using a user-choice solver among the available ones, which are Gauss elimination, LU-decomposition (default) and GMRES-iterative solvers.

-  The GMRES solver code :cite:p:`GMRES` from `CERFACS <https://www.cerfacs.fr/algor/Softs/GMRES/index.html>`__ is embedded in NEMOH solver module. For using the GMRES solver, the user has to obtain a license at https://www.cerfacs.fr/algor/Softs/GMRES/license.html.

-  For free-surface piercing bodies problem, the irregular frequencies removal (IRR) method is applied by the user providing lid panels at :math:`z=0`. Then, the extended boundary integral equation will be solved :cite:p:`Babarit15,Malenica98`. As in :cite:t:`Malenica98`, the IRR may be influenced by the input parameter :math:`\epsilon` in ``input_solver.txt`` that shifts the lid panels from :math:`z=0` to :math:`z=-\epsilon d_B` where :math:`d_B` is a maximum horizontal distance of points on the body. :math:`d_B` is computed by the software.

-  RAO in Eq. :eq:`Eq:RAO` is obtained by applying the inverse matrix using LU-decomposition.

-  The software can solve multi-bodies problems, as well as multi-directional waves.


NEMOH2 uses the following numerical approach

-  The QTF module can be run only after the first order-hydrodynamic coefficients are computed in NEMOH1.

-  In the potential part, the computation of the free-surface integral is an option:

   -  For the difference-frequency QTFs, it is in general acceptable not to compute the free-surface integral terms.

   -  For the sum-frequency QTFs, it is necessary to compute the free-surface integrals.

-  Important notice: the computation with the free-surface integral still has an issue if the lid body panels exist (cf. IRR method). For now, the user is suggested not to specify the lid body panels in the mesh file input for NEMOH1 computation if he wants to compute the full QTFs with the free surface integral.

-  For the free-surface integral, a quadrilateral free-surface mesh has to be specified.

-  The computation can be done for bi-directional or uni-directional wave for the specified multiple wave direction.

-  QTF computations have not been tested yet for the multi-bodies problem.


NEMOH related publications to be referred are :cite:t:`Babarit15` for the first order NEMOH and :cite:t:`Philippe15,Kurnia22_JH,Kurnia22,Kurnia23` for the QTF module.

*****
Units
*****

NEMOH expects all quantities to be expressed in S.I. units: :math:`m, kg, s, rad` (meter, kilogram, seconds, radian, respectively). But some of the phase outputs may be expressed in :math:`deg` or :math:`^{\circ}`, in this case it will be indicated in the file header.

The force unit is [:math:`N`], the moment unit is [:math:`Nm`], added Mass [:math:`kg`], damping coefficient [:math:`kg/s`]. As the force output is normalized with the unit wave amplitude :math:`a` :math:`[m]`, then the normalized force unit is [:math:`N/m`] and the normalized moment is [:math:`N`].

Response amplitude operator for translation motion has unit [:math:`m/m`] and for rotation it is [:math:`deg/m`].

The force quadratic transfer function (QTF) has unit [:math:`N/m^2`] and for the moment QTF it is [:math:`N/m`]. The QTF output is normalized by :math:`\rho g` where the fluid density :math:`\rho,\ [kg/m^3],` and the gravitation constant :math:`g,\ [m/s^2]`.

*****************
Software features
*****************

.. _`fig:flowchart`:
.. figure:: figures/FlowChart.png
   :align: center

   Global flowchart of NEMOH software

:numref:`fig:flowchart` shows a global overview of the software. There are three main programs: a mesh preprocessor, NEMOH1 and NEMOH2. The program features and capabilities are described as follows.

Mesh Preprocessor
=================

NEMOH mesh preprocessor, the executable file ``mesh``, is for generating the NEMOH mesh file with a given geometry input file and an input ``Mesh.cal`` file. This ``mesh`` is not a meshing code but allows the user to refine an existing mesh and to calculate properties such as displacement, buoyancy center, and hydrostatic stiffness. It also makes estimates of masses and inertia matrix. The concept with this program is to write by hand a coarse description of the body under consideration in a ``GeomInput`` file and to have ``mesh`` make the refined mesh for NEMOH calculations.

NEMOH1: 1st-order solver
========================

NEMOH1 solves the first-order potential flow problem. There are four modules: ``preProc``, ``hydrosCal``, ``solver`` and ``postProc``, described as follows.

-  ``preProc``: processes the input mesh file and generates the body condition for each calculation case (diffraction and radiation). The outputs are used as input for ``solver``.

-  ``hydrosCal``: computes hydrostatic parameters, i.e. stiffness matrix and inertia matrix. The output file will be used in the ``postProc`` for computing the RAOs. If the input mesh is generated by the NEMOH mesh preprocessor, ``mesh``, the hydrostatic parameters are already computed and then it is not necessary to execute this program.

-  ``solver``: solves the boundary value problems for each problem, diffraction and radiation, defined in the file ``Normalvelocities.dat``, provided by the ``preProc``.

   -  The influence coefficients matrix is constructed with the infinite/finite depth Green function.

   -  If a finite depth is specified, then the finite depth green function is applied only for :math:`\frac{\omega^2}{g}D<20`, otherwise infinite depth case is applied.

   -  The integration of the Green function on a panel for the influence coefficients is obtained by the Gauss-quadrature integration. The number of Gauss quadrature points is a user input.

   -  The minimum distance, :math:`\epsilon`, between the flow and source points for the influence coefficient computation is user-specified.

   -  The source distributions are then obtained by solving the linear system. There are three options for the solver: Gauss elimination, LU-decomposition and GMRES. If the GMRES solver :cite:p:`GMRES` is used and the target tolerance is not achieved after the maximum number of iterations, the problem is automatically solved by LU-decomposition. License for using GMRES has to be obtained in https://www.cerfacs.fr/algor/Softs/GMRES/license.html.

-  ``postProc``: post-processes the ``solver``\ ’s output files. The results are the excitation forces, added mass and damping coefficients. Optionally, the program computes

   -  the radiation damping impulse response function, the infinite frequency added mass and the excitation force impulse response function,

   -  the Kochin coefficient,

   -  the free-surface elevation,

   -  the motion response amplitude operator (RAO). For the RAO computation, additional stiffness matrix :math:`[\boldsymbol K_m]` and additional damping :math:`[\boldsymbol B_{add}]` can be user-specified in the ``Mechanics/`` folder.

NEMOH2: 2nd-order QTF module
============================

NEMOH2 computes the second-order wave loads that are expressed as Quadratic Transfer Function (QTF). It is suggested to verify the first-order results before running the QTF module. There are three modules in this program: ``QTFpreProc``, ``QTFsolver`` and ``QTFpostProc``, described as follows

-  ``QTFpreProc``: computes the perturbed potential, the total potential, the normalized radiation potential and the corresponding velocities on the body panels, the water-line and the free-surface panels.

   -  The computation on free-surface panels requires possibly long computational time. Then, it is suggested not to compute the free-surface integral for the first execution of NEMOH2. This is controlled by the flag HASFS, which is available in the input file ``Nemoh.cal``.

   -  In general, the free-surface integral may be negligible for the difference-frequency QTFs computation.

   -  The potential on the waterline is rather sensitive with the :math:`\epsilon` value. For default, :math:`\epsilon=0.001`, it can be adjusted in ``input_solver.txt``. The :math:`\epsilon` can be set differently for NEMOH1 and NEMOH2. Further investigation into this is needed.

   -  In case the body lid panels exist, the influence coefficients are affected and give a somewhat larger error for higher frequencies on the free-surface potentials and velocities. This also needs to be investigated.

   -  For now, in the case of full-QTFs computation, the user is suggested not to specify the lid body panels in a mesh file input for NEMOH1 computation.

-  ``QTFsolver``: computes the quadratic part and the potential part of the second order loads. The free-surface integrals in the potential part QTF are optionnally computed (flag HASFS in ``Nemoh.cal``).

-  ``QTFpostProc``: adds all the computed QTF parts and produces the total QTF. The option to sum only some parts of the QTF is available in ``Nemoh.cal``.
