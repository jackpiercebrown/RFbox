RFbox is a Matlab toolbox for the simulation and estimation of random fields.

For theoretical details see PhD thesis 'Effective properties of random heterogeneous beams'

-----
CMD1D
-----
CMD1D.m
CMD1Dcond.m
corMat1D.m
example.m

-----
CMD2D
-----
corMat2D.m
CMD2D.m
example.m
points2D.m

-----
FFT1D
-----
FFT1D
example.m

-----
FFT2D
-----
FFT2D.m
example.m

------
Images
------
All of example images



*README.md

=====
RFbox
=====

-----
CMD
-----
Simulate realisations of a 1D or 2D random field using the covariance matrix decomposition method.
Features: Normal or lognormal marginal distribution, four correlation functions, conditional or unconditional.

<Gaussian correlation function conditional simulation picture>

-----
FFT
-----
Simulate realisations of a 1D or 2D Gaussian random field using the fast Fourier transform method.
Features: Four correlation functions. Much fast than CMD for large fields.

<Markov 2D random field picture>

--------
Estimate
--------
Estimate the parameters of a random field

*File layout

-----
CMD
-----
corMat.m
CMD1D.m
CMD2D.m
example.m
points2D.m

-----
FFT
-----
FFT1D.m
FFT2D.m
example.m

-----
Est
-----
corFun.m
corFunFFT.m
fitCorFun.m
krige.m
lag1cor.m
lag1theta.m
lfunGauss.m
MLE.m

