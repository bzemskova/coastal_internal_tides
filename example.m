% File to accompany manuscript by Zemskova, V.E., Musgrave, R.C. and
% Lerczak, J. A., 
%titled "Internal tides at the coast: energy flux of baroclinic tides propagating into
%the deep ocean in the presence of supercritical shelf topography"


%  Barbara Zemskova (barbara.zemskova@oregonstate.edu)
%  College of Earth, Ocean, and Atmospheric Sciences
%  Oregon State University
%  May, 2023

%  EXAMPLE TO RUN /matlab_functions/func_simulations.m
%       for constant stratification
%           & non-constant stratification (pycnocline depth=800m,
%           width=200m)


%% COMMON INPUTS: 

addpath('./matlab_functions/')

Nx = 1800; %grid points in x
Nz = 300; %grid points in z
R = 0;  %rigid lid
force_type = 0; %Baines body force
l = 0; %along-shore wavenumber; 
L = 300e3; %cross-shore domain extent (m)
h0 = 3100; %max depth (m)
xW = 32e3; %slope width (m)
xs = 80e3; %shelf width (m)
hc = 100;  %depth at the coast (m)
hs = 150;  %depth at shelf break
f = 0.93e-4;  %Coriolis parameter (1/s)
sigma = 1.41e-4; %forcing frequency (1/s), M2 tide
rho0 = 1000; %background density (kg/m^3)
g = 9.81; %gravity (m/s^2)
N2back = (2*pi/(0.5*60*60))^2 ; %background linear stratification N^2 (1/s^2)
eta0 = 0.1; %sea surface elevation at the coast (m)


%% EXAMPLE 1: CONSTANT STRATIFICATION
filename1 = 'constantN2_example.mat';

% SET BOTH TO ZERO FOR CONSTANT STRATIFICATION
Zpyc = 0;  % Pycnocline width (m) MUST BE NEGATIVE
mupyc = 0;  % Pycnocline depth (m) MUST BE POSITIVE


[uModeF,vModeF,pModeF,XFlux,AFlux, fu, fv,...
                Cg, up_modal_offshore, up_modal_integrated, E_u2, Eu2_integrated,...
                xx, z, N2, h, hx, F] = ...
                        func_simulation(Nx, Nz, R, Zpyc, mupyc, force_type, l,...
                        L,h0,xW,xs,hc,hs,f,sigma,N2back,eta0);

figure; 
subplot(4,1,1)
pcolor(xx/1000,z/1e3,sqrt(N2)) ;
shading interp
colorbar
set(gca,'color',[0.8 0.8 0.8])
ylabel('Depth (km)')
title('Buoyancy frequency $N\, (\rm{s}^{-1}$)', ...
    'FontSize',10,'Interpreter','latex')

subplot(4,1,2)
pcolor(xx(:,1:2:end)/1e3,z(:,1:2:end)/1e3,F(:,1:2:end))
shading interp
clim([-max(F(:)) max(F(:))])
colorbar
colormap(gca,bluewhitered(256))
set(gca,'color',[0.8 0.8 0.8])
ylabel('Depth (km)')
title('Forcing $\tilde{F}_B\, (\rm{s}^{-2}$)', ...
    'FontSize',10,'Interpreter','latex')

subplot(4,1,3)
pcolor(xx(:,1:2:end)/1e3,z(:,1:2:end)/1e3,fu)
shading interp
clim([min(fu(:)) 0])
colorbar
colormap(gca,bluewhitered(256))
set(gca,'color',[0.8 0.8 0.8])
ylabel('Depth (km)')
title({'Cross-shore energy flux (W/m^2)'})

subplot(4,1,4)
plot(xx(1,1:2:end)/1e3,XFlux,'LineWidth',2)
xlabel('Cross-shore distance (km)')
ylabel({'Integrated cross-shore', 'energy flux (W/m)'})

save(filename1)


%% EXAMPLE 2: NON-CONSTANT STRATIFICATION
% pycnocline width = 800m,  depth = 200m
filename2 = 'Zpyc800_mupyc200_example.mat';

% SET BOTH TO ZERO FOR CONSTANT STRATIFICATION
Zpyc = -800;  % Pycnocline width (m) MUST BE NEGATIVE
mupyc = 200;  % Pycnocline depth (m) MUST BE POSITIVE


[uModeF,vModeF,pModeF,XFlux,AFlux, fu, fv,...
                Cg, up_modal_offshore, up_modal_integrated, E_u2, Eu2_integrated,...
                xx, z, N2, h, hx, F] = ...
                        func_simulation(Nx, Nz, R, Zpyc, mupyc, force_type, l,...
                        L,h0,xW,xs,hc,hs,f,sigma,N2back,eta0);

figure; 
subplot(4,1,1)
pcolor(xx/1000,z/1e3,sqrt(N2)) ;
shading interp
colorbar
set(gca,'color',[0.8 0.8 0.8])
ylabel('Depth (km)')
title('Buoyancy frequency $N\, (\rm{s}^{-1}$)', ...
    'FontSize',10,'Interpreter','latex')

subplot(4,1,2)
pcolor(xx(:,1:2:end)/1e3,z(:,1:2:end)/1e3,F(:,1:2:end))
shading interp
clim([-max(F(:)) max(F(:))])
colorbar
colormap(gca,bluewhitered(256))
set(gca,'color',[0.8 0.8 0.8])
ylabel('Depth (km)')
title('Forcing $\tilde{F}_B\, (\rm{s}^{-2}$)', ...
    'FontSize',10,'Interpreter','latex')

subplot(4,1,3)
pcolor(xx(:,1:2:end)/1e3,z(:,1:2:end)/1e3,fu)
shading interp
clim([min(fu(:)) 0])
colorbar
colormap(gca,bluewhitered(256))
set(gca,'color',[0.8 0.8 0.8])
ylabel('Depth (km)')
title({'Cross-shore energy flux (W/m^2)'})

subplot(4,1,4)
plot(xx(1,1:2:end)/1e3,XFlux,'LineWidth',2)
xlabel('Cross-shore distance (km)')
ylabel({'Integrated cross-shore', 'energy flux (W/m)'})
