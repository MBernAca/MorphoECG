# MorphoECG
This repository includes Matlab scripts used to perform morphological analyses on fjords in East Central Greenland (65°N-75°N).
It is associated with the companion paper 'The impact of lithology on fjord morphology' by M. BERNARD (Rennes University), P. STEER (Rennes University), K. Gallagher (Rennes University) and D.L. Egholm(Aarhus University).

The script 'Main.m' allows for the extraction of fjord widths and depths, with drainage area computed as a proxy for ice discharge. The code 'plot_results.m' computes the power-law fits for width or depth against drainage area, for each lithological domain.

The 'GrayC.mat' colormap used in this script has been created by Fabio Crameri. Crameri, F. (2018), Scientific colour-maps, Zenodo, doi:10.5281/zenodo.1243862

Authors: Maxime Bernard (Rennes University, France), Philippe Steer (Rennes University, France).

IMPORTANT: to perform the analyses you will need to download first the BedMachine V3 DEM (Morlighem et al., 2017), and the geological map data (Henriksen et al., 2000)

Bedmachine V3: https://sites.uci.edu/morlighem/dataproducts/bedmachine-greenland/
Geological Map of Greenland 1/2500000: https://frisbee.geus.dk/webshop/?customer=nanoq&lang=en
