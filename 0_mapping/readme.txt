This folder contains the code used for the mapping procedure. In this context, the term “mapping” means referencing the coordinate systems of the two spectral channels in the fluorescence microscopy setup to each other.

The concept is explained in Chapter 4 of ‘Single-molecule techniques: a laboratory manual’, edited by Paul R. Selvin & Taekjip Ha (Cold Spring Harbor Laboratory Press, Cold Spring Harbor, NY, USA, 2008)

In TMFM, fluorescent beads with dyes that emit in both spectral channels were immobilized on the sample surfaces by incubation in aqueous solution containing 1 mM Magnesium Chloride. Images in both channels were acquired and taken as input for the Matlab script ‘TWO_COLOR_MAPPING.m’, which makes use of the in-built Matlab function fitgeotrans. (Large portions of the script ‘TWO_COLOR_MAPPING.m’ were written by Jonas Funke and Letizia Meregalli, at the time PhD students in the group of Prof. Dietz, TUM).

The output variables needed for transforming x/y coordinates from one channel to the other in subsequent analysis steps are called ‘tform_2TO1’, ’tform_1TO2’ or simply ’tform’.

- - - - -
Matthias Schickinger
October 2020