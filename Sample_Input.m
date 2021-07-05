% sample input file which contains user defined parameters for simulation
%
% Tianxiang Wu 2021/06/30
% wtx@zju.edu.cn

%-----------------------------------------------------------------------------------
np = 1e4; % Number of photons to be simulated.
nl = 2; % Number of layers (not include the medium above and below)
n = [1,...      % medium above
    1.5,1.5,...% Refractive index for each layer.
    1];  %The first and last elements are for 
         % medium above and below, respectively. So, there should have (nl+2) elements.
mua = [0.2 0.2]; % Absorption coefficients for each layer. [cm-1]
mus = [2 2]; % Scattering coefficients for each layer. [cm-1]
g = [0.91 0.93]; % Anisotropy for each layer.
d = [0.1 0.2]; % Depth for each layer. [cm]

nz = 200; % Number of grids considered in z coordinate, array ranges from 1 to nz.
nr = 200; % Number of grids considered in r coordinate, array ranges from 1 to nr.
na = 50; % Number of grids considered in alpha coordinate, array ranges from 1 to na.

rlim = 1; %Maximum tissue radius [cm]
%------------------------------------------------------------------------------------
