# Kohonen
Fortran2003 library for conventional and Two Levels Self-Organizing Maps.

## Motivation

Self-Organizing Maps or Kohonen Maps are powerful computational tools to cluster multivariate data using a topology preservation approach, that is, the clustering obtained by using this methodology is designed to preserve neighboring relationships between samples (closer samples in the input space remain closer in the output space).

Two Level Self-Organizing Maps are important in hydrological regionalization but currently there is not a proper implementation that suited my needs. So I decided to implement this clustering approach in Fortran using the features included in the new standard Fortran2003. This ensures computational efficiency and software extensibility, in addition to having to learn more about Fortran. 



## Structure 

