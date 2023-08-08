function [XFlux,AFlux,fu,fv] = calcFluxFunc(uModeF,vModeF,pModeF,h,nu)
% function [XFlux,AFlux] =
%       calcFluxFunc(uModeF,vModeF,pModeF,h,nu)
% uses pressure wave field to derive along and across shore velocity
% components, then computes energy flux in each direction. Flux is in W/m.
% Ruth Musgrave Oct. 2018 Barbara Zemskova Nov. 2022

% dimensions: [Nx, Nz]

%INPUTS:
%     uModeF[Nx/2,Nz]: cross-shore velocity (from calcModeFunc)
%     vModeF[Nx/2,Nz]: along-shore velocity (from calcModeFunc)
%     pModeF[Nx/2,Nz]: pressure (from calcModeFunc)
%     h[Nx/2,1]: depth 
%     nu[Nx/2,1]: even-numbered points

%OUTPUTS:
%     fu[Nx/2,Nz]: cross-shore flux (m^2/s)
%     fv[Nx/2,Nz]: along-shore flux (m^2/s)
%     XFlux[Nx/2,1]: vertically-integrated cross-shore flux (m^3/s)
%     AFlux[Nx/2,1]: vertically-integrated along-shore flux (m^3/s)
    


% Cross shore flux [kg /s^3]/[m / s^2]/[kg / m^3] = [m^2 / s]
fu = 0.5*(real(uModeF).*real(pModeF) +...
    imag(uModeF).*imag(pModeF)); 
% Along shore flux [kg /s^3] = W/m^2
fv = 0.5*(real(vModeF).*real(pModeF) +...
    imag(vModeF).*imag(pModeF)); 


% depth integral, [m^3 / s]
XFlux = mean(fu,1).*h(nu);  %cross-shore
AFlux = mean(fv,1).*h(nu);  %along-shore


