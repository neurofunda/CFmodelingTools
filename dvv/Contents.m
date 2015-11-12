% The Delay Vector Variance (DVV) method uses local predictibility  in phase space to examine
% the determinism and nonlinearity within signals. 
% It yields DVV plots in which the target variance is plotted as a standardised distance.
% DVV scatter plots can also be generated, where the horizontal 
% axis corresponds to the original time series while the vertical axis to that of surrogate series;
% nonlinearity within the signal can be detected if the 
% scatter plot deviates from the bisector line. This toolbox implements the DVV method in terms of
% nonlinearity testing for some basic signals. 
%
%
%   A Delay Vector Variance (DVV) toolbox for MATLAB
%   (c) Copyright Danilo P. Mandic 2008
%   http://www.commsp.ee.ic.ac.uk/~mandic/dvv.htm
%
%
% The following Matlab files are used to perform DVV and demonstrate its operation:
%
% - analysis.m      Matlab code to perform the DVV analysis for the case studies stored in *.mat files
% - dvv.m           Implementation of the DVV method for real valued and complex signals
% - surrogate.m     Matlab code which generates surrogate data for real and complex signals
%
%
% ar.mat          Stored linear AR2 Signal
% wind.mat        Stored 2D wind Signal
% henon.mat       Stored nonlinear henon map signal
