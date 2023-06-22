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
%  May, 2023

%% Set-up parameters

addpath('./matlab_functions')

%grid points in x
Nx_vec = 600:200:3400;
%grid points in z
Nz = [10, 20, 30, 40, 50, 75, 100, 150, 200, 250, 300]; 

Nx = 3400; %grid points in x
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

f_vec = linspace(0.5,1.35,50)*1e-4;

Zpyc = -600; %pycnocline depth (m) %must be negative
mupyc = 200; %pycnocline width (m) %must be positive

figure;

%% Sweep over Nx (horizontal resolution), keep Nz constant
subplot(2,2,1)
hold on
for inx = 1:length(Nx_vec)
    
    [~,~,~,XFlux,~, ~, ~,...
                ~, ~, ~, ~, ~,...
                xx, ~, ~, ~, ~, ~] = ...
                        func_simulation(Nx_vec(inx), Nz, R, Zpyc, mupyc, force_type, l,...
                        L,h0,xW,xs,hc,hs,f,sigma,N2back,eta0);


    plot(xx(1,1:2:end)/1e3,XFlux/1e3,'LineWidth',2)
    
end
legend('600','800','1000','1200', '1400','1800', '2000',...
          '2200', '2400');
xlabel('Cross-shore distance (km)')
ylabel('Cross-shore energy flux (kW/m)')
title('Vary Nx')

%% Sweep over Nz (vertical resolution), keep Nx constant

subplot(2,2,2)
hold on
for inz = 1:length(Nz_vec)
    
    [~,~,~,XFlux,~, ~, ~,...
                ~, ~, ~, ~, ~,...
                xx, ~, ~, ~, ~, ~] = ...
                        func_simulation(Nx, Nz_vec(inz), R, Zpyc, mupyc, force_type, l,...
                        L,h0,xW,xs,hc,hs,f,sigma,N2back,eta0);


    plot(xx(1,1:2:end)/1e3,XFlux/1e3,'LineWidth',2)
    
end
legend('10', '20', '30', '40', '50', ...
    '75', '100','150', '200', '250', '300');
xlabel('Cross-shore distance (km)')
ylabel('Cross-shore energy flux (kW/m)')
title('Vary Nz')

%% Sweep over Nx (horizontal resolution), keep Nz constant
% this time varying Coriolis parameter f
subplot(2,2,3)
hold on
for inx = 1:length(Nx_vec)
    for ifx = 1:length(f_vec)
    [~,~,~,XFlux,~, ~, ~,...
                ~, ~, ~, ~, ~,...
                xx, ~, ~, ~, ~, ~] = ...
                        func_simulation(Nx_vec(inx), Nz, ...
                        R, Zpyc, mupyc, force_type, l,...
                        L,h0,xW,xs,hc,hs,f_vec(ifx), ...
                        sigma,N2back,eta0);

    XF(ifx) = XFlux(1);
    end
    plot(sqrt(sigma^2-f_vec.^2),XF/1e3,'LineWidth',2)
end
legend('600','800','1000','1200', '1400','1800', '2000',...
          '2200', '2400','2600','2800','3000','3200','3400');
xlabel('$\sqrt{\omega^-f^2}\,(\rm{s}^{-1})$', ...
    'Interpreter','latex')
ylabel('Cross-shore energy flux (kW/m)')
title('Vary Nx')

%% Sweep over Nz (vertical resolution), keep Nx constant

subplot(2,2,4)
hold on
for inz = 1:length(Nz_vec)
    for ifx = 1:length(f_vec)
    [~,~,~,XFlux,~, ~, ~,...
                ~, ~, ~, ~, ~,...
                xx, ~, ~, ~, ~, ~] = ...
                        func_simulation(Nx, Nz_vec(inz), ...
                        R, Zpyc, mupyc, force_type, l,...
                        L,h0,xW,xs,hc,hs,f_vec(ifx), ...
                        sigma,N2back,eta0);

    XF(ifx) = XFlux(1);
    plot(sqrt(sigma^2-f_vec.^2),XF/1e3,'LineWidth',2)
    end
end
legend('10', '20', '30', '40', '50', ...
    '75', '100','150', '200', '250', '300');
xlabel('$\sqrt{\omega^-f^2}\,(\rm{s}^{-1})$', ...
    'Interpreter','latex')
ylabel('Cross-shore energy flux (kW/m)')
title('Vary Nz')