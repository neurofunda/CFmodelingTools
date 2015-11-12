% Matlab code to perform the DVV analysis for the case studies stored in *.mat files
%
% An example program to perform DVV analysis on some predefined real-valued and complex signals
% Choose from linear AR(2) signal, nonlinear henon signal, and the real world wind data.
%
%
%   A Delay Vector Variance (DVV) toolbox for MATLAB
%   (c) Copyright Danilo P. Mandic 2008
%   http://www.commsp.ee.ic.ac.uk/~mandic/dvv.htm
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   http://www.gnu.org/copyleft/gpl.html or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc; clear all; close all 

% Initializations 
% DVV specific parameters 
m=3;            % embedding parameter
N_DV=200;       % number of reference DVs to consider
nd=3;           % Span over which to perform DVV   
Ntv=50;         % number of points on horizontal axes

% Number of surrogates to consider
Ns=25;

% Input: Choose any one of the following input signal for DVV analysis 
%load ar.mat;                    % Stored linear AR2 signal
%load henon.mat;               % Stored henon signal
load wind.mat;                % Stored 2D/complex wind signal

% DVV analysis of the original time series
disp('Starting DVV analysis on original series ...')
dvv_orig = dvv(X, m, N_DV, nd, Ntv);
% dvv_orig = dvv(X);

% Generates surrogate data
disp('Generating surrogate time series ...')
Xs = surrogate(X, Ns);

% DVV analysis of surrogate data
disp('Starting DVV analysis on Surrogate series ...')
for acc = 1:size(Xs,2)
    dvv_surr(:,:,acc) = dvv(Xs(:,acc),m, N_DV, nd, Ntv);
%     dvv_surr(:,:,acc) = dvv(Xs(:,acc));
end

% DVV Plots 
figure(1);
plot(dvv_orig(:,1),dvv_orig(:,2),'b-'); grid on; shg; hold on;
title('DVV Plot'); xlabel('std. distance'); ylabel('variance');

% Mean variance values of surrogates
avg_dvv_surr = mean(dvv_surr,3);
plot(avg_dvv_surr(:,1), avg_dvv_surr(:,2),'r-'); grid on; legend('Original','Surrogates'); shg

% DVV scatter Plots 
bisector = 0:0.01:1;
figure(2); 
plot(bisector, bisector, 'k-.'); axis([ 0 1 0 1]); 
title('DVV Scatter Plot'); xlabel('Original'); ylabel('Surrogates'); grid on; shg; hold on
errorbar(dvv_orig(:,2),avg_dvv_surr(:,2),std(dvv_surr(:,2,:),0,3));