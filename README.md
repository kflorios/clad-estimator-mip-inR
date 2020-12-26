# clad-estimator-mip-inR

Exact computation of Censored Least Absolute Deviations estimator with Mixed Integer Programming via R code

Use handy R code in main.R in order to exactly compute the CLAD estimator with MIP.
The main code sets up the matrices A,b,c,Aeq,beq,lb,ub of the MIP model,
and then relies on a MIP solver to solve it. We can use CPLEX and the Rcplex package.
All used functions are supplied in the file functions.R

The dataset is read in readXyw.R function via the files X.txt and ys.txt which can be adopted as desired.
Currently, left censoring at zero is supposed, as is in most applications of CLAD.
In order to have a more flexible modeling approach, readers are suggested to consult the GAMS version
of the same model in https://www.gams.com/modlib/libhtml/clad.htm.

In case you encounter an issue while trying to install Rcplex, do not hesitate to contact me for support.

This code is a translation of the MATLAB code supplied by me at the repository https://github.com/kflorios/clad-estimator-mip

For completeness, I also supply the Manual of the MATLAB version of the software in this repository also.

Feedback for the R code at: cflorios@central.ntua.gr, cflorios@aueb.gr

Suggested Reference:  

Bilias, Y., Florios, K., & Skouras, S. (2019). Exact computation of Censored
Least Absolute Deviations estimator. Journal of Econometrics, 212(2), 584-606.

https://www.sciencedirect.com/science/article/abs/pii/S030440761930140X
