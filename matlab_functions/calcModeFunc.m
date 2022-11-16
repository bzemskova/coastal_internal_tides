function [uModeF,vModeF,pModeF] = ...
    calcModeFunc(pModeC,l,f,sigma,xx,z,h,np,nu,rho0)
% function [uModeF,vModeF] = calcModeFunc(pModeC,l,f,sigma,xx,z,h,np,nu,rho0)
% uses pressure wave field to derive along and 
% across shore velocity components in m/s.
% Ruth Musgrave Oct. 2018
% Barbara Zemskova Nov. 2022

% dimensions: [Nx, Nz]

% INPUT:
%     pModeC[Nx, Nz]: solution for pressure (on odd-numbered grid points in x)
%         and cross-shore velocity (on even-numbered grid points in x)
%     l[1,1]: along-shore wavenumber (1/m)
%     sigma[1,1]: forcing frequency
%     xx[Nx,Nz]: x-coordinates (cross-shore direction)
%     z[Nx,Nz]: z-coordinates (depth, sigma coordinates)
%     h[Nz,1]: depth at each cross-shore location
%     np[Nx/2,1]: odd-numbered grid points in x
%     nu[Nx/2,1]: even-numbered grid points in x
%     rho0[1,1]: background density

% OUTPUT:
%     uModeF[Nx/2,Nz]: dimensionalized cross-shore velocity (m/s)
%                    (on even-numbered grid points in x)
%     vModeF[Nx/2,Nz]: dimensionalized along-shore velocity (m/s)
%                    (on even-numbered grid points in x)
%     pModeF[Nx/2,Nz]: dimensionalized pressure 
%                    (on even-numbered grid points in x)



I = sqrt(-1);
% initialize some variables
pModeF = NaN*pModeC; 
dpdxF = NaN*pModeC;  
dpdzF = NaN*pModeC;
dhdx = NaN*h;
if mean(pModeC(:,end),1)<0
    pModeC=-pModeC;
end

x = xx(1,:);
xp = x(np); dxp = mean(diff(xp));

% get pressure on faces
pModeF(:,1:end-1) = 0.5*(pModeC(:,1:end-1)+pModeC(:,2:end)); 
% interpolate for end point
pModeF(:,end) = pModeF(:,end-1)+(l*f/sigma)*pModeC(:,end)*dxp;

% on u grid, uncorrected for curvature
dpdxF(:,1:end-1) = (pModeC(:,2:end) - pModeC(:,1:end-1))/dxp; 
% impose the boundary condition
dpdxF(:,end) = (l*f/sigma)*pModeF(:,end); 

% metric terms for horizontal derivative
dpdzF(2:end-1,:) = (pModeF(3:end,:)-pModeF(1:end-2,:))./(z(3:end,nu)-z(1:end-2,nu));
dpdzF(1,:) = (pModeF(2,:)-pModeF(1,:))./(z(2,nu)-z(1,nu));
dpdzF(end,:) = (pModeF(end,:)-pModeF(end-1,:))./(z(end,nu)-z(end-1,nu));        

dhdx(2:end-1) = (h(3:end)-h(1:end-2))./(x(3:end)-x(1:end-2));
dhdx(1) = (h(2)-h(1))/(x(2)-x(1));
dhdx(end) = (h(end)-h(end-1))/(x(end)-x(end-1));

dzdx = bsxfun(@times,z,dhdx./h);

dpdxF2 = dpdxF - dpdzF.*dzdx(:,nu); % correct for curvature

% calculate cross and along-shore modes
% Cross shore on u grid
uModeF = (I*sigma/(f^2-sigma^2))*(dpdxF2 - (f/sigma)*l*pModeF); 
% Along shore on u grid
vModeF = (sigma/(f^2-sigma^2))*((f/sigma)*dpdxF2 - l*pModeF);   

