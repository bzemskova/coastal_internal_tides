% File to accompany manuscript by Zemskova, V.E., Musgrave, R.C. and
% Lerczak, J. A., 
%titled "Internal tides at the coast: energy flux of baroclinic tides propagating into
%the deep ocean in the presence of supercritical shelf topography"

%File to run an example simulation with specified grid spacing, topography,
%       and stratification parameters 

%   Also can specify Coriolis parameter, forcing frequency, along-shore
%   wavenumber, forcing type (Baines vs. isolated), surface boundary
%   condition (rigid lid vs. free linear surface)
%        
%   Outputs velocity and pressure fields plus energy flux fields
%   (cross-shore, along-shore; 2D and vertically integrated values)

%   Plots relevant fields (topography, stratification, topographic
%   criticality, 2D forcing, cross-shore energy flux) 

%  Barbara Zemskova
%  College of Earth, Ocean, and Atmospheric Sciences
%  Oregon State University
%  barbara.zemskova@oregonstate.edu
%  November, 2022


%% Set-up parameters

addpath('./matlab_functions')

Nx = 4800; %grid points in x
Nz = 300; %grid points in z
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

filename_mat = 'example_data_600_200.mat';

%% Run simulation
[XFlux,AFlux,fu,fv,uModeF,pModeF,vModeF,...
    up_modal_offshore,up_modal_all,up_modal_integrated,...
    xx,z,N2,h,hx,F] = ...
    func_simulation(Nx, Nz, R, Zpyc, mupyc, force_type, l,...
    L,h0,W,xs,hc,hs,f,sigma);



%% Plot Baines forcing function and cross-shore energy flux
figure;
subplot(3,1,1)
pcolor(xx(:,1:2:end)/1e3,z(:,1:2:end),F(:,1:2:end))
shading interp
caxis([-max(F(:)) max(F(:))])
colorbar
colormap(gca,bluewhitered(256))
set(gca,'color',[0.8 0.8 0.8])
xlabel('Cross-shore distance (km)')
ylabel('Forcing')

subplot(3,1,2)
pcolor(xx(:,1:2:end)/1e3,z(:,1:2:end),fu/1e3)
shading interp
caxis([min(fu(:)/1e3) 0])
colorbar
colormap(gca,bluewhitered(256))
set(gca,'color',[0.8 0.8 0.8])
xlabel('Cross-shore distance (km)')
ylabel({'Cross-shore', 'energy flux (kW/m^2)'})

subplot(3,1,3)
plot(xx(1,1:2:end)/1e3,XFlux/1e3,'LineWidth',2)
xlabel('Cross-shore distance (km)')
ylabel({'Integrated cross-shore', 'energy flux (kW/m)'})

%% Save variables to .mat 

save(filename_mat)
