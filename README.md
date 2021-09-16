# Laminate
Effective behavior and stress distribution in a laminate with arbitrary number of phases

This MATLAB code computes the effective behavior (effective stiffness tensor and the effective plastic strain tensor) and the stress fields in an N-phase laminate with planar interface. 
The laminate is submitted to a macroscopic homogeneous and remotely applied stress as well as to piecewise uniform plastic strains. 
The phases of the laminate can have different grain volume fractions and can correspond to different materials or crystals
The phases are assumed perfectly bonded  along the planar interface whose normal is along e2. 
The code below is an application for a beta-titanium alloy with elongated grains loading in elasticity only,  modeled as laminate made of 100,000 different orientations and equal volume fraction.
At the end, the maximal von Mises stress is plotted with respect to theta, a rotation around e3 of the uniaxial stress of magnitude 300 MPa.
A reference for the model can be found in:

The contracted Voigt notation (11:1, 22: 2, 33: 3, 23: 4, 31:5, 12: 6) is used for stresses and strains in vector notation and stiffnesses and compliances in 6x6 matrix notation.
