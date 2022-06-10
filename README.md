# Surface Water and Ocean Topography Data Assimilation With Correlated Error Reduction
The Surface Water and Ocean Topography (SWOT) satellite altimeter will provide unprecedented high-resolution two-dimensional wide-swath sea surface height (SSH) data. The SWOT altimetric data is expected to be affected by specially correlated errors. This study proposes a procedure that reduces the correlated error by embedding the Correlated Error Reduction (CER) into data assimilation scheme, and solve for the correlated SWOT errors as part of the assimilation. We compare the the performance of linear regression and data-driven machine learning approach in reducing correlated errors. With a Rossby wave test model, we design a series of experiments that compare linear regression (LR) and machine learning (ML) data assimilation techniques. We incorporate white noise and correlated along-track noise into the Rossby wave model, and reconstructing SSH field with LR model and Artificial Neural Networks (ANNs), respectively.


## Discretized wave projection problem
This project solves the discretized wave projection problem, given the vertical profiles of Temperature, Salinity, Pressure and depth inteval length.

VERT_FSFB2.m: 
MATLAB script that solves the discretized wave projection problem. Gabriel A. Vecchi - May 12, 1998

VERT_FSFB3.m: 
MATLAB script that solves the discretized wave projection problem. Youran Li (yol039@ucsd.edu) - 2022

VERT_FSFB3.py: Yu Gao (yug032@ucsd.edu) - April 18 2022

CER00_VERT_FSFB3-N2-constant.ipynb: Vertical mode decomposition of constant buoyancy frequency squared, with barotropic mode included.

CER00_VERT_FSFB3-N2-constant-removeBT.ipynb: Vertical decomposition mode of constant buoyancy frequencysquared, with barotropic mode removed.

CER01_VERT_FSFB3.ipynb: Vertical decomposition mode of a sample buoyancy frequency squared in the California coastal region, with barotropic mode removed.

CER02_VERT_FSFB3-function-usage.ipynb :  Demonstration of how to use the VERT_FSFB3.py function.
