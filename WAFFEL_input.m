%%%%  WafFEL simulation code, User entered parameters %%%%
clearvars
close all

param.accuracy = 1e-4;

param.progress = 0;
param.tapering = 0;

%% physical constants
global c mu0 e0 me eps0 IA
c = 2.99792458.*10^8;                                             % speed of light
e0 = 1.60217657e-19;                                              % electron charge
me = 9.10938291e-31;                                              % electron mass
eps0 = 8.85418782e-12;                                            % eps_0
mu0 = 1.256637e-6;                                                % mu_0
IA = 17045;                                                       % Alfven current

%% radiation parameters
param.frozenfield = 0;                                            % Frozen field approximation 1 = On/ 0 = Off
param.order2 = 0;
param.nfreq = 51;                                                 % number of spectral points, keep at 1 for genesis template 
param.lambda0 = 3.57e-4;                                          % seed wavelength
param.E0 = 0*12e6;                                                   % Peak field (V/m)  
param.waist = .001;                                               % radiation waist size, focused at entrance to undulator/waveguide  
param.deltaBW = .95;                                              % set spectral bandwidth i.e. from w0-deltaBW*w0, to w0+deltaBW*w0
param.sw = 2e12;                                                  % Input spectrum width for gaussian
param.radphase = pi;
param.timingTHz = -.1e-12; % 1e-12;                                       % THz time of arrival offset

%% beam parameters
param.gaussian_beam = 1;                                          % set to 1 for gaussian current profile
param.spacecharge = 1;                                            % Longitudinal space charge field On/Off
param.prebunching = 0;                                            % set to 1 to start from a pre-bunched beam. 
param.Np = 2000;                                                  % # of macroparticles (500-1000 well) 
param.Ee = 5.8e6;                                                 % Total e-beam energy (eV)
param.bunchlength = .2E-12;                                      % RMS bunchlength (s)
param.bunch = 0.0;                                                % Initial bunching factor
param.bunchphase = 0.0;                                           % Initial bunching phase
param.sigmax = 200.0e-6;                                          % beam radius
param.charge = 8e-12;                                                   % beam charge 
param.deltagammarel = .05;                                              % Percent energy spread
param.Edist = 2;                                                  % 0: gaussian, 1: chirped gaussian, 2: measured spectrum
param.tfocus = .3;
param.chirp = 0;
%% undulator and waveguide parameters
param.lambdau = 3e-2;       % undulator period
param.Kset = 1.27;
%[linspace(2,2,14),linspace(1.7,1,param.Nsteps+1-14)]; %1.27; %   

param.und_periods = 10;                                           % number of undulator periods to simulate
param.b = .0024; %.00215;                                                   % PPWG plate spacing
param.Rppwg = .002;                                               % PPWG plate radius of curvature
%%

WAFFEL_THz_field;
WAFFEL_beam_dist;
WAFFEL_core;
WAFFEL_output;

