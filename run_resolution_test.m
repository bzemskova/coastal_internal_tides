% File to accompany manuscript by Zemskova, V.E., Musgrave, R.C. and
% Lerczak, J. A., 
%titled "Internal tides at the coast: energy flux of baroclinic tides propagating into
%the deep ocean in the presence of supercritical shelf topography"

%File to run a resolution test (for vertical and horizontal resolutions)
%       for an example simulation with specified grid spacing, topography,
%       and stratification parameters 

%   Also can specify Coriolis parameter, forcing frequency, along-shore
%   wavenumber, forcing type (Baines vs. isolated), surface boundary
%   condition (rigid lid vs. free linear surface)
%        
%   Outputs and plots vertically-integrated cross-shore energy flux for
%       different vertical and horizontal resolutions to assess convergence.


%  Barbara Zemskova
%  College of Earth, Ocean, and Atmospheric Sciences
%  Oregon State University
%  barbara.zemskova@oregonstate.edu
%  November, 2022

%% Set-up parameters

addpath('./matlab_functions')
load('colorblind_colormap.mat')

%grid points in x
Nx_vec = [300,600,900,1200, 1500,2000, 2500,...
          3000, 3600, 4200, 4800, 5400  ];
%grid points in z
Nz = [30, 40, 50, 70, 90, 110, 150, 200, 250, 300]; 

R = 0; %rigid lid
Zpyc = -600; %pycnocline depth
mupyc = 200; %pycnocline width
force_type = 0; %Baines body force
l = 7e-9; %along-shore wavenumber; small value close to 0
L = 800e3; %cross-shore domain extent (m)
h0 = 3100; %max depth (m)
W = 32e3; %slope width (m)
xs = 80e3; %shelf width (m)
hc = 100; %depth at the coast (m)
hs = 150; %depth at shelf break
f = 9.3e-5; %Coriolis parameter (1/s)
sigma = 1.41e-4; %forcing frequency (1/s), M2 tide
rho0 = 1000; %background density (kg/m^3)
g = 9.81; %gravity (m/s^2)

%% Sweep over Nx (horizontal resolution), keep Nz constant

figure;
subplot(1,2,1)
hold on
colororder(colorblind)
for inx = 1:length(Nx_vec)
    
    [XFlux,~,~,~,~,~,~,...
    ~,~,~,...
    xx,~,~,~,~,~] = ...
    func_simulation(Nx_vec(inx), max(Nz_vec), R, Zpyc, mupyc, force_type, l,...
    L,h0,W,xs,hc,hs,f,sigma);

    plot(xx(1,1:2:end)/1e3,XFlux/1e3,'LineWidth',2)
    
end
legend('300','600','900','1200', '1500','2000', '2500',...
          '3000', '3600', '4200', '4800', '5400');
xlabel('Cross-shore distance (km)')
ylabel('Cross-shore energy flux (kW/m)')
title('Vary Nx')

%% Sweep over Nz (vertical resolution), keep Nx constant

subplot(1,2,2)
hold on
colororder(colorblind)
for inz = 1:length(Nz_vec)
    
    [XFlux,~,~,~,~,~,~,...
    ~,~,~,...
    xx,~,~,~,~,~] = ...
    func_simulation(max(Nx_vec), Nz_vec(inz), R, Zpyc, mupyc, force_type, l,...
    L,h0,W,xs,hc,hs,f,sigma);

    plot(xx(1,1:2:end),XFlux,'LineWidth',2)
    
end
legend('30', '40', '50', '70', '90', '110', '150', '200', '250', '300');
xlabel('Cross-shore distance (km)')
ylabel('Cross-shore energy flux (kW/m)')
title('Vary Nz')