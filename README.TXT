==============
=== CREDIT ===
==============

Basic 2D monocular visual SLAM assuming a bottom-looking camera and flat floor.

Loop detection is achieved by means of RANSAC and SIFT features. IEKF is used to optimize the motions graph.

Author : Antoni Burguera
e-mail : antoni.burguera@uib.es

This software was developed and used in a paper currently under review. Please, let us know if you plan to use this software so we can cite us with the proper reference.

=========================
=== BUILD THE PACKAGE ===
=========================

The software is pure matlab code. Just put it wherever you want and use it.

However, the code uses the VLFeat library. If you don't have VLFeat, download it from:

http://www.vlfeat.org/

Depending on how you install VLFeat, it may be necessary to run (from matlab) vl_setup.m from the vlfeat/toolbox directory.

====================
=== RUN THE CODE ===
====================

The main funcion is:

[Xslam,Pslam,theLoops]=slam(odoData,drawRobot,drawLoops,onlinePlot,finalPlot)

Except for odoData, the remaining input parameters control what to plot during execution. Check slam.m for more information.

As for odoData, it is an array of structures containing odometric estimates and SIFT features. The exact format is the one provided by the function compute_odometry from the VISUAL_ODOMETRY_2D package. Check this link for more information:

https://github.com/aburguera/VISUAL_ODOMETRY_2D

Having a properly formatted odoData, running the software is as simple as typing the following from the matlab prompt:

>> [Xslam,Pslam,theLoops]=slam(odoData,1,1,1,1);

The output is:

- Xslam: An N*3 sized column vector, being N the number of frames in odoData. Each item Xslam(i*3-2:i*3,1) contains the relative motion (x displacement, y displacement and change in heading) from frame i-1 to frame i. The first item Xslam(1:3,1) is set to zero by definition.
- Pslam: The N*3xN*3 covariance matrix of Xslam.
- theLoops: List of detected loops. It is an 2xM matrix, M being the number of loops. The meaning is that odoData(theLoops(1,k)) and odoData(theLoops(2,k)) close a loop for all 1<=k<=M.

The function draw_slam can be used to plot this data.

=======================
=== TROUBLESHOOTING ===
=======================

This software package has been tested using Matlab 9.0.0.341360 (R2016a) running on Ubuntu 16.04LTS. Other configurations are untested.
