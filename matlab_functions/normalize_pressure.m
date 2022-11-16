function [uModeF,vModeF,pModeF,XFlux,AFlux,fu,fv] = ...
    normalize_pressure(PP,l,f,sigma,xx,z,h,np,nu,eta_coast)

%function to dimensionalize and normalize cross-shore velocity and pressure
% output and also compute energy fluxes

% dimensions: [Nx, Nz]

%INPUTS: 
%     PP[Nx, Nz]: solution for pressure (on odd-numbered grid points in x)
%         and cross-shore velocity (on even-numbered grid points in x)
%     l[1,1]: along-shore wavenumber (1/m)
%     f[1,1]: Coriolis parameter
%     sigma[1,1]: forcing frequency
%     xx[Nx,Nz]: x-coordinates (cross-shore direction)
%     z[Nx,Nz]: z-coordinates (depth, sigma coordinates)
%     h[Nx,1]: depth at each cross-shore location
%     np[Nx/2,1]: odd-numbered grid points in x
%     nu[Nx/2,1]: even-numbered grid points in x
%     eta_coast[1,1]: sea-surface elevation at the coast (=0.1 m)


% OUTPUT:
%     uModeF[Nx/2,Nz]: dimensionalized cross-shore velocity (m/s)
%                    (on even-numbered grid points in x)
%     vModeF[Nx/2,Nz]: dimensionalized along-shore velocity (m/s)
%                    (on even-numbered grid points in x)
%     pModeF[Nx/2,Nz]: dimensionalized pressure 
%                    (on even-numbered grid points in x)
%     fu[Nx/2,Nz]: cross-shore flux (W/m^2)
%     fv[Nx/2,Nz]: along-shore flux (W/m^2)
%     XFlux[Nx/2,1]: vertically-integrated cross-shore flux (W/m)
%     AFlux[Nx/2,1]: vertically-integrated along-shore flux (W/m)


% background density and gravity
rho0 = 1000;
g = 9.8;

% normalize raw pressure output such that sea-surface elevation
%   at the coast = eta_coast
%   (take eta_coast = 0.1 m)
pModeC = PP(:,np);
scl = mean(abs(pModeC(:,end))); 
pModeC = eta_coast*pModeC/scl;

%computes cross-shore and off-shore fluxes
[uModeF,vModeF,pModeF] = calcModeFunc(pModeC,l,f,sigma,xx,z,h,np,nu,rho0);
[XFlux,AFlux,fu,fv] = calcFluxFunc(uModeF,vModeF,pModeF,h,nu);

% convert flux values to W/m^2 and vertically-integrated flux values to W/m
XFlux = g*rho0*XFlux;
AFlux = g*rho0*AFlux;
fu = g*rho0*fu;
fv = g*rho0*fv;


