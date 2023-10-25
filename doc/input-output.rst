############
Input/Output
############

Following is the list of the user’s input files and the output files for each executable file:

-  **``mesh``**

   -  User’s input files: ``projdir/Mesh.cal``, ``projdir/geomInput``,

   -  Output files: ``projdir/meshfile``, ``projdir/mesh/[*.dat, *.tec]``,

-  **``preProc``**

   -  User’s input files: ``projdir/Nemoh.cal``, ``projdir/meshfile``,

   -  Output files: ``projdir/Normalvelocities.dat``,
      ``projdir/results/[FKForce.dat, FKForce.tec, index.dat]``,
      ``projdir/mesh/[L10.dat, L12.dat]``,

-  **``hydrosCal``**

   -  User’s input files: ``projdir/Nemoh.cal``, ``projdir/mesh.cal``,

   -  Output files: ``projdir/mesh/[*.dat, *.tec]``, ``projdir/Mechanics/``,

-  **``solver``**

   -  User’s input files: ``projdir/Nemoh.cal``, ``projdir/input_solver.txt``

   -  Output files: ``projdir/results/Forces.dat``,
      ``projdir/results/sources`` (if QTF will be computed, indicated in ``Nemoh.cal``)

-  **``postProc``**

   -  User’s input files: ``projdir/Nemoh.cal``,
      ``projdir/Mechanics/[Km.dat,Badd.dat]``

   -  Output files: ``projdir/results/[ExcitationForce.tec,``
      ``DiffractionForce.tec,RadiationCoefficients.tec,``\ :math:`\cdots`\ ``]``,
      ``projdir/Motion/RAO.dat``,

-  **``QTFpreProc``**

   -  User’s input files: ``projdir/Nemoh.cal``, ``projdir/FSmeshfile`` (If the free-surface integral, HASFS flag, is computed),
      ``projdir/Mechanics/[Km.dat,Badd.dat]``

   -  Output files: ``projdir/QTFPreprocOut/*.bin``,

-  **QTFsolver**

   -  User’s input files: ``projdir/Nemoh.cal``

   -  Output files: ``projdir/results/QTF/*.dat``,

-  **QTFpostproc**

   -  User’s input files: ``projdir/Nemoh.cal``

   -  Output files: ``projdir/results/QTF/[OUT_QTFM_N.dat,OUT_QTFP_N.dat]``


As summary, following files are needed for the input, some depends on user-specified choice in ``Nemoh.cal``:

-  ``Nemoh.cal`` contains all NEMOH computation parameters

-  ``Mesh.cal`` contains information of ``geomInput`` file. It is an input for **``mesh``** and **``hydrosCal``**.

-  a ``meshfile``, input for **``preProc``**, or ``geomInput`` file, input for **``mesh``**

-  ``input_solver.txt`` contains **``solver``** parameters

-  ``Km.dat`` and ``Badd.dat``, are the additional stiffness and damping matrices. These optional input are for **``postProc``**/**``QTFpreProc``**

-  ``FSmeshfile`` contains the free-surface mesh if contrib=3 in ``Nemoh.cal``.


NEMOH produces the following main output files, some depending on user-specified choice in ``Nemoh.cal``,

-  hydrostatic files: inertia and stiffness matrices,

-  hydrodynamics coefficients: Froude-Krylov force, excitation force, added-mass, damping coefficient,

-  Kochin function, free-surface elevation,

-  Response Amplitude Operator (RAO),

-  Total difference- and sum-frequencies QTFs.


Detail descriptions of the input/output files are discussed in the next subsections.

******************
User’s input files
******************

.. table:: ``Nemoh.cal`` input file
   :name: tab:NemohCal

   =============== ===== ===== ===== ===== ===== ===== ==================================================================
   File contents                                       Signification
   =================================================== ==================================================================
   \--- Environment ---------------------------------- *Section header*
   --------------------------------------------------- ------------------------------------------------------------------
   1025\.                                              Fluid density :math:`\rho` :math:`[kg/m^3]`
   9.81                                                Gravitional acceleration :math:`g` :math:`[m/s^2]`
   200\.                                               Water depth :math:`d` :math:`[m]`
   0\.             0\.                                 Wave measurement location :math:`(x_{eff},y_{eff})` :math:`[m]`
   \--- Description of floating bodies --------------- *Section header*
   --------------------------------------------------- ------------------------------------------------------------------
   1                                                   Number of bodies
   \--- Body 1 --------------------------------------- *Section header*
   --------------------------------------------------- ------------------------------------------------------------------
   mesh.dat                                            Name of mesh file
   657             610                                 Number of nodes and number of panels in mesh
   6                                                   Number of degrees of freedom
   1               1\.   0\.   0\.   0\.   0\.   0\.   Surge
   1               0\.   1\.   0\.   0\.   0\.   0\.   Sway
   1               0\.   0\.   1\.   0\.   0\.   0\.   Heave
   2               1\.   0\.   0\.   0\.   0\.   -5\.  Roll about a point (here :math:`(0,0,-5.)`)
   2               0\.   1\.   0\.   0\.   0\.   -5\.  Pitch about a point (here :math:`(0,0,-5.)`)
   2               0\.   0\.   1\.   0\.   0\.   -5\.  Yaw about a point (here :math:`(0,0,-5.)`)
   ...                                                 *This line is repeated for each degree of freedom*
   6                                                   Number of resulting generalised forces
   1               1\.   0\.   0\.   0\.   0\.   0\.   Force in x direction
   1               0\.   1\.   0\.   0\.   0\.   0\.   Force in y direction
   1               0\.   0\.   1\.   0\.   0\.   0\.   Force in z direction
   2               1\.   0\.   0\.   0\.   0\.   -5\.  Moment force in x direction about a point (here :math:`(0,0,-5.)`)
   2               0\.   1\.   0\.   0\.   0\.   -5\.  Moment force in y direction about a point (here :math:`(0,0,-5.)`)
   2               0\.   0\.   1\.   0\.   0\.   -5\.  Moment force in z direction about a point (here :math:`(0,0,-5.)`)
   ...                                                 *This line is repeated for each generalised forces*
   0                                                   Number of lines of additional information
   ...                                                 *This section is repeated for each body*
   \--- Load cases ----------------------------------- *Section header*
   --------------------------------------------------- ------------------------------------------------------------------
   1               100   0.062 6.28                    Input frequency unit: 1 for *rad/s*, 2 for *Hz* and 3 for *s*, Number of wave frequencies, Min, and Max
   2               0     30                            Number of wave directions, Min and Max (degrees)
   \--- Post processing ------------------------------ *Section header*
   --------------------------------------------------- ------------------------------------------------------------------
   0               0.1   10\.                          IRF (Impulse Response Function) flag (0/1), time step and duration
   0                                                   Pressure output flag (0/1)
   0               0\.   180\.                         Kochin functions: number of directions of calculation (0 to deactivate), Min and Max (degrees)
   0               50    400\. 400\.                   Free surface elevation output: number of points in x direction (0 to deactivate) and y direction and (x,y) dimensions of domain
   1                                                   RAO (Response Amplitude Operator) flag (0/1)
   1                                                   Output requency unit, 1 for *rad/s*, 2 for *Hz* and 3 for *s*
   \-- QTF ------------------------------------------- *Section header*
   --------------------------------------------------- ------------------------------------------------------------------
   1                                                   QTF (Quadratic Transfer Function) flag (0/1)
   65              0.062 4.082                         Number of radial frequencies, Min, and Max values for the QTF computation
   1                                                   Bidirectional QTF computation flag (0/1)
   2                                                   Contributing terms: 1 DUOK, 2 DUOK+HASBO, 3 Full QTF (DUOK+HASBO+HASFS+ASYMP)
   NA                                                  Name of free surface meshfile (**only for full QTF**), 'NA' if not applicable
   0               0     0                             Free surface QTF parameters: Re, Nre and NBessel (**only for full QTF**)
   0                                                   Include hydrostatic terms of the quadratic first order motion (:math:`-[\boldsymbol K] \tilde{\boldsymbol\xi}^{(2)}`), flag (0/1)
   1                                                   For QTFposProc: output frequency unit, 1 for *rad/s*, 2 for *Hz* and 3 for *s*
   1                                                   For QTFposProc: include DUOK in total QTF, flag (0/1)
   1                                                   For QTFposProc: include HASBO in total QTF, flag (0/1)
   0                                                   For QTFposProc (**only for full QTF**): include HASFS+ASYMP in total QTF, flag (0/1)
   =============== ===== ===== ===== ===== ===== ===== ==================================================================

``Nemoh.cal``: contains all computation parameters with the format as in :numref:`tab:NemohCal`. The following parameters have to be specified:

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

      Note that this term is optional and needed only in QTFsolver. In other software this term is not always included, *e.g.* HYDROSTAR :cite:p:`HYDROSTAR` does not included it.

   -  Frequency output option: 1 for the radial frequency [rad/s], 2 for the frequency [Hz] and 3 for the period [s].

   -  In **``QTFpostproc``**, QTFs total is calculated with summation of all the terms. Option to exclude/include the terms are available with the corresponding switch for DUOK, HASBO and HASFS+ASYMP terms. Switch 1 to include, 0 to exclude. If Contrib=2, then HASFS+ASYMP switch has to be 0.



``meshfile``: contains all the mesh information with a format as shown in :numref:`tab:meshfile`. Lid panels (:math:`z=0`) of the structure may be included in this file to activate the irregular frequencies removal method. This mesh file may be generated by NEMOH **``mesh``** preprocessor or by an external mesh generator.
External mesh generators, *e.g.* the open-source software GMSH :cite:p:`GMSH`, may be used to generate mesh files but they must be adapted to the NEMOH format. A Matlab file for converting GMSH mesh file to the NEMOH format is provided in the dedicated repository. The Matlab file will be described in the next section.

.. table:: ``meshfile`` format
   :name: tab:meshfile

   ======= ============= ============= ============= ============================================================================================================
   File contents                   Signification
   ================================================= ============================================================================================================
   2       1                                         First column must be a 2. Second column is 1 for a symmetric (about :math:`xOz`) body half-mesh, 0 otherwise.
   1       :math:`x_1`   :math:`y_1`   :math:`z_1`   Table of nodes: first column is the node ID, other 3 are the coordinates :math:`(x,y,z)` of each node, listed as rows.
   ...     ...           ...           ...
   0       0\.           0\.           0\.           Last line of the table of nodes.
   1       2             3             4             Table of connectivities: node IDs of each panel listed as rows.
   ...     ...           ...           ...
   0       0             0             0             Last line of the table of connectivities.
   ======= ============= ============= ============= ============================================================================================================

``geomInput``: contain coarse description of mesh, that are number of nodes, number of panels, table of nodes and table of connectivities. The input file has to follow the format as shown in :numref:`tab:geomInput`.

.. table:: ``geomInput`` file format
   :name: tab:geomInput

   ============= ============= ============= ==== ==================================================================
   File contents                Signification
   ============================================== ==================================================================
   100                                            Total number of nodes.
   25                                             Total number of panels.
   :math:`x_1`   :math:`y_1`   :math:`z_1`        Table of nodes: coordinates :math:`(x,y,z)` of each node listed as rows.
   ...           ...           ...           ...
   1             2             3             4    Table of connectivities: node IDs of each panel listed as rows.
   ...           ...           ...           ...
   ============= ============= ============= ==== ==================================================================

``Mesh.cal:`` contains mesh and environmental parameters with a format as in :numref:`tab:meshcal`. This file is used as input for **``mesh``** and **``hydroCal``**. All the parameters are used in **``mesh``**. Only center of gravity, water density, and gravity are used in **``hydroCal``**.

.. table:: ``Mesh.cal`` file format
   :name: tab:meshcal

   =============== === === ==================================================================
   File contents           Signification
   ======================= ==================================================================
   geomInput_name          Name of the geomInput file.
   0                       1 for a symmetric (about :math:`xOz`) body half-mesh, 0 otherwise.
   0\.             0\.     Translation about x and y axis (respectively)
   0\.             0\. -7  Coordinates of gravity centre
   500\.                   Target for the number of panels in refined mesh
   2\.
   0\.
   1\.
   1025                    Water density :math:`(kg/m^3)`
   9.81                    Gravity acceleration :math:`(m/s^2)`
   =============== === === ==================================================================

``input_solver.txt`` contains solver parameters with format as in Table :numref:`tab:input_solver`. The parameters are described as follows.

-  Number of Gauss Quadrature points, :math:`N^2`, is used for the surface integration in the influence coefficients. User specifies an integer value of :math:`N\in [1,4]`, default :math:`N=2`.

-  Minimum z of flow and source points is defined with a factor :math:`\epsilon_{zmin}` multiplied by the maximal horizontal distance between two point of the mesh, default :math:`\epsilon_{zmin}=0.001`.

-  Three linear-system solvers are available; 1 Gauss elimination, 2 LU Decomposition, 3 GMRES iterative solver.

-  If GMRES solver is chosen then the three parameters, the restart parameter, the relative tolerance and the maximum number of iterations, have to be specified. If the tolerance is not achieved after the maximum iteration exceeded then LU decomposition solves the system directly.

.. table:: ``input_solver.txt`` file format
   :name: tab:input_solver

   =========== ===== ===== =====================================================================================
   File contents           Signification
   ======================= =====================================================================================
   2                       Gauss quadrature order N=\[1,4\] for surface integration, resulting in :math:`N^2` nodes.
   0.001                   eps_zmin for determining minimum z of flow and source points of panel.
   1                       Solver option: 0 for GAUSS ELIM., 1 for LU DECOMP., 2 for GMRES.
   10          1e-5  1000  GMRES parameters: restart parameter, relative tolerance and max number of iterations.
   =========== ===== ===== =====================================================================================

``Km.dat`` and ``Badd.dat`` are additional stiffness matrix and damping coefficient matrix. The files contains the matrix components with size :math:`(Nbody\cdot Nradiation)\times (Nbody\cdot Nradiation)`.
``FSmeshfile`` contains all the free-surface mesh information with a format as shown in :numref:`tab:FSmeshfile`. Quadrilateral panels discretized free-surface area in between the body waterline, :math:`R_B`, and the exterior radius :math:`R_e`. Waterline on :math:`R_B` and :math:`R_e` has to discretized by line segments.

.. table:: ``FSmeshfile`` format (Free surface mesh file)
   :name: tab:FSmeshfile

   ======= ============= ============= ============= ============================================================================================================
   File contents                   Signification
   ================================================= ============================================================================================================
   1       5000          4900           400          Free-surface computation parameters: first column is 1 for a symmetric (about :math:`xOz`) body half-mesh, 0 otherwise. Column 2-4 are number of nodes, number of panels and number of segments for the waterline, respectively.
   1       :math:`x_1`   :math:`y_1`   :math:`z_1`   Table of nodes: first column is the node ID, other 3 are the coordinates :math:`(x,y,z)` of each node, listed as rows.
   ...     ...           ...           ...
   0       0\.           0\.           0\.           Last line of the table of nodes.
   1       2             3             4             Table of connectivities: node IDs of each panel listed as rows.
   ...     ...           ...           ...
   4901    4902                                      Table of connectivities for the waterline: node IDs of each segment listed as rows.
   ...     ...           ...           ...
   0       0             0             0             Last line of the table of connectivities.
   ======= ============= ============= ============= ============================================================================================================

************
Output files
************

Hydrostatic output files such as inertia and stiffness matrices are produced by **``mesh``**, if ``geomInput`` is prescribed, or by **``hydroCal``**, if ``meshfile`` is prescribed. The files contain the matrix components with size :math:`(Nbody\cdot Nradiation)\times (Nbody\cdot Nradiation)`.

The following hydrodynamic coefficients are produced in Tecplot format, which can be opened by the Tecplot program or by a simple text-editor program,

-  ``FKForce.tec``, ``DiffractionForce.tec`` and ``ExcitationForce.tec`` are the output files of the Froude-Krylov, the diffraction and the excitation forces respectively. The output file format is given in :numref:`tab:WaveForce`. The file contains the absolute value and the phase [deg] of the force for each ’frequency’ :math:`f`. The force is given for each specified force axis (i.e. surge, heave, pitch) for each body. The ’frequency’ is given based on the chosen type, [rad/s, Hz, s], of the post-processing parameter in ``Nemoh.cal``, except the Froude-Krylov force, which is only in the radial frequency [rad/s].

-  ``RadiationCoefficients.tec`` is the output file for added mass and damping coefficients with format as in Table :numref:`tab:addedmass_damping_coeffs`. The radiation coefficients are given for each :math:`DoF`, each force axis and for each frequency. The frequency is given based on the chosen ’frequency’ type, [rad/s, Hz, s], of the post-processing parameter in ``Nemoh.cal``.

The hydrodynamic coefficients are also produced in the *.dat* files, i.e. *CA.dat* for the damping coefficients, *CM.dat* for the added mass coefficients, *Fe.dat* for the excitation force and *FKForce.dat* for the excitation force. The frequency type of the output files is only radial frequency [rad/s]. These output files are used as input files for the QTF module.

.. table:: Output file format of Froude-Krylov, diffraction and excitation forces
   :name: tab:WaveForce

   ================ ======================= ============================ ================ ================ ============================== ================================
   :math:`f_1`      :math:`|F_1(f_1)|`      :math:`\angle F_1(f_1)`      :math:`\cdots`   :math:`\cdots`   :math:`|F_{Ninteg}(f_1)|`      :math:`\angle F_{Ninteg}(f_1)`
   :math:`f_2`      :math:`|F_1(f_2)|`      :math:`\angle F_1(f_2)`      :math:`\cdots`   :math:`\cdots`   :math:`|F_{Ninteg}(f_2)|`      :math:`\angle F_{Ninteg}(f_2)`
   :math:`\vdots`   :math:`\vdots`          :math:`\vdots`               :math:`\vdots`   :math:`\vdots`   :math:`\vdots`                 :math:`\vdots`
   :math:`f_{Nf}`   :math:`|F_1(f_{Nf})|`   :math:`\angle F_1(f_{Nf})`   :math:`\cdots`   :math:`\cdots`   :math:`|F_{Ninteg}(f_{Nf})|`   :math:`\angle F_{Ninteg}(f_{Nf})`
   ================ ======================= ============================ ================ ================ ============================== ================================


.. table:: Output file format of the radiation coefficients
   :name: tab:addedmass_damping_coeffs

   ==================== ================================= ============================== ================ ================ ====================================== ================================
   :math:`f_1`          :math:`M^a_{11}(f_1)`             :math:`B_{11}(f_1)`            :math:`\cdots`   :math:`\cdots`   :math:`M^a_{1Ninteg}(f_1)`             :math:`B_{1Ninteg}(f_1)`
   :math:`f_2`          :math:`M^a_{11}(f_2)`             :math:`B_{11}(f_2)`            :math:`\cdots`   :math:`\cdots`   :math:`M^a_{1Ninteg}(f_2)`             :math:`B_{1Ninteg}(f_2)`
   :math:`\vdots`       :math:`\vdots`                    :math:`\vdots`                 :math:`\vdots`   :math:`\vdots`   :math:`\vdots`                         :math:`\vdots`
   :math:`f_{N_f}`      :math:`M^a_{11}(f_{N_f})`         :math:`B_{11}(f_{N_f})`        :math:`\cdots`   :math:`\cdots`   :math:`M^a_{1Ninteg}(f_{N_f})`         :math:`B_{1Ninteg}(f_{N_f})`
   :math:`f_1`          :math:`M^a_{21}(f_1)`             :math:`B_{21}(f_1)`            :math:`\cdots`   :math:`\cdots`   :math:`M^a_{2Ninteg}(f_1)`             :math:`B_{2Ninteg}(f_1)`
   :math:`\vdots`       :math:`\vdots`                    :math:`\vdots`                 :math:`\vdots`   :math:`\vdots`   :math:`\vdots`                         :math:`\vdots`
   :math:`f_{N_f}`      :math:`M^a_{21}(f_{N_f})`         :math:`B_{21}(f_{N_f})`        :math:`\cdots`   :math:`\cdots`   :math:`M^a_{2Ninteg}(f_{N_f})`         :math:`B_{2Ninteg}(f_{N_f})`
   :math:`\vdots`       :math:`\vdots`                    :math:`\vdots`                 :math:`\vdots`   :math:`\vdots`   :math:`\vdots`                         :math:`\vdots`
   :math:`f_{N_f}`      :math:`M^a_{N_{DoF}1}(f_{N_f})`   :math:`B_{N_{DoF}1}(f_{N_f})`  :math:`\cdots`   :math:`\cdots`   :math:`M^a_{N_{DoF}Ninteg}(f_{N_f})`   :math:`B_{N_{DoF}Ninteg}(f_{N_f})`
   ==================== ================================= ============================== ================ ================ ====================================== ================================

``RAO.dat`` is the output file of the response amplitude operator with the file format as in Table. :numref:`tab:RAO`. The output file gives the absolute value and the phase of RAO for each degree of freedom and each frequency. The frequency is given based on the chosen ’frequency’ type, [rad/s, Hz, s], of the post-processing parameter in ``Nemoh.cal``. Only radial frequency output file will be produced in the case of the QTF computed.

.. table:: Output file format of ``RAO.dat``
   :name: tab:RAO

   ====================== ========================== ================ ========================== =============================== ================ =========================
   :math:`f_1`            :math:`|\xi_1(f_1)|`       :math:`\cdots`   :math:`|\xi_6(f_1)|`       :math:`\angle \xi_1(f_1)`       :math:`\cdots`   :math:`\angle \xi_6(f_1)`
   :math:`\vdots`         :math:`\vdots`             :math:`\vdots`   :math:`\vdots`             :math:`\vdots`                  :math:`\vdots`   :math:`\vdots`
   :math:`f_{N_f}`        :math:`|\xi_1(f_{N_f})|`   :math:`\cdots`   :math:`|\xi_6(f_{N_f})|`   :math:`\angle \xi_1(f_{N_f})`   :math:`\cdots`   :math:`\angle \xi_6(f_{N_f})`
   ====================== ========================== ================ ========================== =============================== ================ =========================

``IRF.tec`` and ``IRF_excForce.tec`` are the impulse response functions for the radiation damping and the excitation force, respectively. The radiation damping IRF has the file format as in :numref:`tab:IRF` and the excitation force IRF as in :numref:`tab:IRFExcF`.

.. table:: Output file format of ``IRF.tec``
   :name: tab:IRF

   ================== ================================ ============================= ================ =================== ===================================== ================================
   :math:`t_1`        :math:`M^a_{11}(\infty)`         :math:`IRF_{11}(t_1)`         :math:`\cdots`   :math:`\cdots`      :math:`M^a_{1Ninteg}(\infty)`         :math:`IRF_{1Ninteg}(t_1)`
   :math:`t_2`        :math:`M^a_{11}(\infty)`         :math:`IRF_{11}(t_2)`         :math:`\cdots`   :math:`\cdots`      :math:`M^a_{1Ninteg}(\infty)`         :math:`IRF_{1Ninteg}(t_2)`
   :math:`\vdots`     :math:`\vdots`                   :math:`\vdots`                :math:`\vdots`   :math:`\vdots`      :math:`\vdots`                        :math:`\vdots`
   :math:`t_1`        :math:`M^a_{21}(\infty)`         :math:`IRF_{21}(t_1)`         :math:`\cdots`   :math:`\cdots`      :math:`M^a_{2Ninteg}(\infty)`         :math:`IRF_{2Ninteg}(t_1)`
   :math:`\vdots`     :math:`\vdots`                   :math:`\vdots`                :math:`\vdots`   :math:`\vdots`      :math:`\vdots`                        :math:`\vdots`
   :math:`t_N`        :math:`M^a_{N_{DoF}1}(\infty)`   :math:`IRF_{N_{DoF}1}(t_N)`   :math:`\cdots`   :math:`\cdots`      :math:`M^a_{N_{DoF}Ninteg}(\infty)`   :math:`IRF_{N_{DoF}Ninteg}(t_N)`
   ================== ================================ ============================= ================ =================== ===================================== ================================


.. table:: Output file format of ``IRF_excForce.tec``
   :name: tab:IRFExcF

   ================ ====================== ================ ========================
   :math:`t_1`      :math:`IRF_{1}(t_1)`   :math:`\cdots`   :math:`IRF_{Ninteg}(t_1)`
   :math:`t_2`      :math:`IRF_{1}(t_2)`   :math:`\cdots`   :math:`IRF_{Ninteg}(t_2)`
   :math:`\vdots`   :math:`\vdots`         :math:`\vdots`   :math:`\vdots`
   :math:`t_N`      :math:`IRF_{1}(t_N)`   :math:`\cdots`   :math:`IRF_{Ninteg}(t_N)`
   ================ ====================== ================ ========================

``pressure.00XXX.dat``, ``kochin.00XXX.dat`` and ``freesurface.00XXX.dat`` are output files of pressure, Kochin and free surface, respectively, for a specific problem-XXX. The problem number is defined as in order of the diffraction problem (:math:`Nbeta`), the radiation problem (:math:`Ndof`) and for each frequency. So problem-001 is the, first frequency and first wave direction, diffraction problem. Suppose :math:`Nbeta=1`, then problem-002 is the first frequency radiation problem DoF 1. If :math:`Ndof=6` then problem-008 is the second frequency diffraction problem.

-  ``pressure.00XXX.dat`` is a pressure output file for the problem-XXX. In each file, the absolute value of pressure, :math:`|P|`, (Pa) and the phase, :math:`\angle P`, (rad) are given for each panel. The format of the output file is given in :numref:`tab:pressure`.

.. table:: Output file format of ``pressure.00XXX.dat``
   :name: tab:pressure

   ==================== ==================== ==================== ===================================== ====================================
   :math:`x_1`          :math:`y_1`          :math:`z_1`          :math:`|P(\boldsymbol x_1)|`          :math:`\angle P(\boldsymbol x_1)`
   :math:`\vdots`       :math:`\vdots`       :math:`\vdots`       :math:`\vdots`                        :math:`\vdots`
   :math:`x_{Npanel}`   :math:`y_{Npanel}`   :math:`z_{Npanel}`   :math:`|P(\boldsymbol x_{Npanel})|`   :math:`\angle P(\boldsymbol x_{Npanel})`
   ==================== ==================== ==================== ===================================== ====================================

-  ``kochin.00XXX.dat`` is an output file of the Kochin function on a prescribed direction for the problem-XXX. In each file, depending on the diffraction/radiation problem, the computed absolute value of the Kochin, :math:`|\mathcal{H}|`, and the phase, :math:`\angle \mathcal{H}`, (rad) are saved for each direction, :math:`\vartheta`. The format of the output file is given in :numref:`tab:kochin`.

.. table:: Output file format of *kochin.00XXX.dat*
   :name: tab:kochin

   ================================ =============================================== ===============================================
   :math:`\vartheta_1`              :math:`|\mathcal{H}(\vartheta_1)|`              :math:`\angle \mathcal{H}(\vartheta_1)`
   :math:`\vdots`                   :math:`\vdots`                                  :math:`\vdots`
   :math:`\vartheta_{N\vartheta}`   :math:`|\mathcal{H}(\vartheta_{N\vartheta})|`   :math:`\angle \mathcal{H}(\vartheta_{N\vartheta})`
   ================================ =============================================== ===============================================

-  ``freesurface.00XXX.dat`` is an output file of the free-surface elevation on a prescribed free-surface domain for the problem-XXX. In each file, depending on the diffraction/radiation problem, the computed absolute value of the free-surface elevation, :math:`|\eta|`, and the phase, :math:`\angle \eta`, (rad) are saved for each free-surface panel position. The format of the output file is given in :numref:`tab:freesurface`.

.. table:: Output file format of ``freesurface.00XXX.dat``
   :name: tab:freesurface

   ==================== ==================== ================================== ======================================= ===================================== ==================================
   :math:`x_1`          :math:`y_1`          :math:`|\eta(\vec{x}_1)|`          :math:`\angle \eta(\vec{x}_1)`          :math:`Re[ \eta(\vec{x}_1)]`          :math:`Im[ \eta(\vec{x}_1)]`
   :math:`\vdots`       :math:`\vdots`       :math:`\vdots`                     :math:`\vdots`                          :math:`\vdots`                        :math:`\vdots`
   :math:`x_{Npanel}`   :math:`y_{Npanel}`   :math:`|\eta(\vec{x}_{Npanel})|`   :math:`\angle \eta(\vec{x}_{Npanel})`   :math:`Re[ \eta(\vec{x}_{Npanel})]`   :math:`Im[ \eta(\vec{x}_{Npanel})]`
   ==================== ==================== ================================== ======================================= ===================================== ==================================

``OUT_QTFM_N.dat`` and ``OUT_QTFP_N.dat`` are the output files of difference- and sum-frequencies QTF. The QTF results are either the total QTF or parts of the QTF terms that depend on the user choice QTF post-processing parameters in ``Nemoh.cal``. The QTF values are given in the absolute value with the phase in deg and real-imaginary parts. The QTF values are normalized by :math:`\rho g`. The ’frequency’ type, [rad/s, Hz, s], depends on the user choice in the ``Nemoh.cal``. The format of the output file is given in :numref:`tab:QTF`. Only the lower triangular part of the QTF matrix is saved in the file. The full difference-frequency QTF matrix can be constructed with the lower triangular part of the matrix and the upper triangular part which is in conjugate-symmetric with the lower part. The upper triangular part of the sum-frequency QTF is symmetric with the lower triangular part. A Matlab file for reading this output file is provided in ``matlabRoutines/`` and will be described in the next section.

.. table:: Output file format of ``OUT_QTFM_N.dat`` and ``OUT_QTFP_N.dat``
   :name: tab:QTF

   ==================== ====================== =========================== =========================== ==================== ====================== ==================== ======================== ======================
   :math:`f_{1_1}`      :math:`f_{2_1}`        :math:`\beta_{1_1}`         :math:`\beta_{2_1}`         :math:`DoF_1`        :math:`|QTF|/\rho g`   :math:`\angle QTF`   :math:`Re[QTF]/\rho g`   :math:`Im[QTF]/\rho g`
   :math:`f_{1_2}`      :math:`f_{2_1}`        :math:`\beta_{1_1}`         :math:`\beta_{2_1}`         :math:`DoF_1`        :math:`|QTF|/\rho g`   :math:`\angle QTF`   :math:`Re[QTF]/\rho g`   :math:`Im[QTF]/\rho g`
   :math:`\vdots`       :math:`\vdots`         :math:`\vdots`              :math:`\vdots`              :math:`\vdots`       :math:`\vdots`         :math:`\vdots`       :math:`\vdots`           :math:`\vdots`
   :math:`f_{1_{Nf}}`   :math:`f_{2_{Nf}}`     :math:`\beta_{1_{Nbeta}}`   :math:`\beta_{2_{Nbeta}}`   :math:`DoF_{NDof}`   :math:`|QTF|/\rho g`   :math:`\angle QTF`   :math:`Re[QTF]/\rho g`   :math:`Im[QTF]/\rho g`
   ==================== ====================== =========================== =========================== ==================== ====================== ==================== ======================== ======================
