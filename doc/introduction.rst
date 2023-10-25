############
Introduction
############

This document is the Manual of NEMOH Software in version v3.0 Released on 2nd December 2022. It serves as a guide for using, installing and running the software.

NEMOH in its original version is an open-source potential flow boundary element solver for computing first-order hydrodynamic coefficients in the frequency domain. NEMOH v3.0 includes the treatment of irregular frequencies, has an extended module to post-process the first-order hydrodynamic results, and compute complete Quadratic Transfer Functions (QTFs).
The new developments/extensions in this NEMOH v3.0 are summarised as follows.

-  GNU General Public v3 license is used instead of Apache.

-  Newly written source codes for better readability and coding practice.

-  Green function uses finer interpolation points.

-  The influence coefficient is constructed by applying Gauss-quadrature integration over a panel.

-  Irregular frequency removal method is available.

-  Two new options of the linear system solvers, LU decomposition and GMRES iterative solver, are available for enhanced computational efficiency.

-  Full difference- and sum-frequencies QTF module.

The present User Manual is organized as follows. Section `2 <#Sec:Descrip_NEMOH>`__ describes briefly the mathematical background and the capabilities of the software; it is advised to read this Section before continuing to the rest of the manual. Section `3 <#Sec:Getstarted>`__ describes the installation procedure, how to run the codes, and includes a description of the input and output files. Section `4 <#Sec:MatlabFiles>`__ describes supporting Matlab files such as a NEMOH wrapper, a mesh converter and the routines to post-process the test cases. Finally, Section `5 <#Sec:Testcase>`__ describes briefly test cases that will show the capabilities of the code and its use. For new NEMOH users, they can be used as the first elements for setting up a new configuration.
