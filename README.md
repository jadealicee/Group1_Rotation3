# Group1_Rotation3
This repository consists of the code used to complete Rotation 3 of the Bioinformatics Group Project module. Used by both Franklin Brown and Jade Gregory.

## Overview
Studies on both ‘pure’ lyrata and arenosa have already been undergone. Instead, the purpose of this project was to look into various lyrata/arenosa hybrids. To conclude, we wanted to look into the divirsity of the hybrids - and whether they are more diverse than 'pure' arenosa and lyrata.

## Objectives
The first aim of this project was to sense-check the dataset. This was done by utilising three separate methods of analysis: fastSTRUCTURE, PCA, and SplitsTree. The second aim was to determine an ideal cohort by looking at the allele frequency spectra of the data provided, and to decide whether hybrid-types have a higher rate of diversity than their ‘pure’ counterparts.

## Expected Outcome
The expected outcome of this project was that hybrids would be more diverse than regular lyrata or arenosa. This is due to the fact that hybridisation is known to stabilise meiosis in these species, which would improve chances of survival.

## Contents
Dependencies required to run the code in this repository are noted in the file '[Dependencies](Dependencies.json)'. Code credits are provided in the file labelled '[Bibliography](Bibliography)'.

__The scripts in this Repository are to be run in the following order:__
1) [fastSTRUCTURE_prep.py](fastSTRUCTURE_prep.py) which utilises both 'structure.py' and 'Cochlearia_create_structure_file.py'
2) [PCA.R](PCA.R) which allows for both eigen analysis as well as a PCA scatter graph to be produced
3) [tree.R](tree.R) which will produce a neighbour-joining tree
