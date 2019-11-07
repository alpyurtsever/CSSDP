-------------------------------------------------------------------------------------------
					CSSDP
		An Optimal-Storage Approach to Semidefinite Programming 
			 using Approximate Complementarity
-------------------------------------------------------------------------------------------

This toolbox includes the code to reproduce numerical results in [DYCTU2019].

Follow <https://github.com/alpyurtsever/CSSDP> for the updates. 

Please contact "alp.yurtsever@epfl.ch" or "ld446 at cornell dot edu" for your questions/comments.


IMPORTANT NOTE: "Test1_MaxCut.m" script downloads "G1" Sparse Matrix Dataset from 
"https://www.cise.ufl.edu/research/sparse/mat/Gset/G1.mat". 
Please read the README file for this dataset, for the description of data and the terms of 
use. You can find the README file at 
"https://www.cise.ufl.edu/research/sparse/mat/Gset/README.txt".
By running "Test1_MaxCut", you agree all the terms and conditions for downloading and using 
this dataset. 

-------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------

DEPENDENCIES

This toolbox requires CVX and MOSEK. We tested the implementation on MATLAB R2018a. 
It might not work on other (especially earlier) versions of MATLAB.

Visit "http://cvxr.com/cvx/" to download CVX.

Visit "http://www.mosek.com/" to download MOSEK.

-------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------

INSTRUCTIONS

This package includes 3 main files for the experiments. Descriptions below:

"Test1_MaxCut" performs the experiments for reproducing Figure 1 in [DYCTU2019].

"Test2_MaxCut" performs the experiments for MaxCut plots of Figure 2 in [DYCTU2019].

"Test2_MatrixCompletion" performs the experiments for Matrix Completion plots of Figure 2 
in [DYCTU2019].

All these three files save the data to reproduce the plots under "./results" folder. Once 
the experiments are completed, you can generate the plots by running "PlotFig1_MaxCut", 
"PlotFig2_MaxCut", and "PlotFig2_MatrixCompletion" scripts. 

-------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------

INDEX

"./data/MatrixCompletionData.mat"
- contains the random data that we generate and use for the Matrix Completion experiments.

"./data/GenerateMatrixCompletionData.m"
- script to generate "MatrixCompletionData.mat". 

"./methods/AdaGradSDP.m"
- function to run AcceleGrad [LYC2018], for solving the model SDP template.

"./methods/AdaGradSDP.m"
- function to run AdaGrad [DHS2011] as described in Algorithm 1 in [L2017], for solving 
the model SDP template.

"./methods/AdaNGDSDP.m"
- function to run AdaNGD [L2017], for solving the model SDP template.

"./utils/spmult.m"
- function to compute sparse matrix - vector multiplication.

"./utils/spmult.c"
- mex function to replace "spmult.m". You need to compile this function for your system.

"./Test1_MaxCut.m"
- script to perform the experiments for Fig 1 [DYCTU2019]. 

"./Test2_MaxCut.m"
- script to perform the MaxCut experiments for Fig 2 [DYCTU2019]. 

"./Test2_MatrixCompletion.m"
- script to perform the Matrix Completion experiments for Fig 2 [DYCTU2019]. 

"./PlotFig1_MaxCut.m"
- script to generate Fig 1 [DYCTU2019]. Run "Test1_MaxCut.m" first. 

"./PlotFig2_MaxCut.m"
- script to generate 1st part of Fig 2 [DYCTU2019]. Run "Test2_MaxCut.m" first. 

"./PlotFig2_MatrixCompletion.m"
- script to generate 2nd part of Fig 2 [DYCTU2019]. Run "Test2_MatrixCompletion.m" first. 

"./README.txt" 
- this README file. 

- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

References

* [DYTCU2019] L. Ding, A. Yurtsever, V. Cevher, J.A. Tropp, M. Udell,
"An Optimal-Storage Approach to Semidefinite Programming using Approximate Complementarity"
arXiv:1902.03373, 2019.

[DHS2011] J. Duchi, E. Hazan, Y. Singer,
"Adaptive Subgradient Methods for Online Learning and Stochastic Optimization" 
Journal of Machine Learning Research 12, 2011.

[L2017] K. Levy, 
"Online to Offline Conversions, Universality and Adaptive Minibatch Sizes" 
Advances in Neural Information Processing Systems 30, (NeurIPS 2017).

[LYC2018] K. Levy, A. Yurtsever, V. Cevher,
"Online Adaptive Methods, Universality and Acceleration" 
Advances in Neural Information Processing Systems 31, (NeurIPS 2018).

-------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------

COPYRIGHT

Copyright (C) 2019

This package (CSSDP) implements the numerical experiments for the Optimal-Storage Approach 
to Semidefinite Programming using Approximate Complementarity proposed in [DYTCU2019]. 

CSSDP is free software: you can redistribute it and/or modify it under the terms of the 
GNU General Public License as published by the Free Software Foundation, version 3 of the 
License.

CSSDP is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. 
If not, see <http://www.gnu.org/licenses/>.

-------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------

CITATION

If you find this toolbox useful in your research, please cite our work:

[DYTCU2019] L. Ding, A. Yurtsever, V. Cevher, J. Tropp, M. Udell,
"An Optimal-Storage Approach to Semidefinite Programming using Approximate Complementarity"
arXiv:1902.03373, 2019.

-------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------

Last edit: Alp Yurtsever - November 6, 2019