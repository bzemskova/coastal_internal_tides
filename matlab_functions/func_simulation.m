%  Barbara Zemskova
%  College of Earth, Ocean, and Atmospheric Sciences
%  Oregon State University
%  barbara.zemskova@oregonstate.edu
%  August, 2022
%
%  This code solves trapped and scattered internal waves in a rotating,
%  stratified basin with continental margin that is infinitely long and
%  does not vary in the along-coast direction.  The method follows that of
%  Dale et al. (2001), solving a system of two equations for pressure
%  (normalized by rho_0) and cross-shore velocity (u).  The assumed
%  solution propagates along the coast as a plane wave -- exp(-i*sigma*t + i*l*y).
%  A positive along-coast wavenumber, l, corresponds to wave propagation in
%  the positive y-direction (the positive x-direction is toward the coast).
%
%  Unlike Dale01, I do not nondimensionalize the equations.  Also, unlike 
%  Dale01, I allow non-evanescent interal Poincare wave modes to radiate
%  away from the coast (at super-inertial frequencies).
%
%  The code works for both sub- and super-inertial frequencies and can be
%  forced by an obliquely incident Poincare wave mode (set the 'force'
%  variable to one, incident wave mode is set by variable 'p').  This will
%  only work at super-inertial frequencies and for a wave that is
%  consistent with the free Poincare wave dispersion relationship (i.e.,
%  the mode can't be evanescent).
%
%  The code will also solve for the ocean response to weak forcing on the
%  shelf (set the 'force' variable to zero).  This will work for both sub-
%  and super-inertial frequencies.  The structure of this forcing is set by
%  the variable 'F' (see the get_initial.m script).

% Specifically, the set-up here follows the test code from J. Klymak,
% which in itself follows Dale and Sherwin (1996).
% Stratification is constant with coastal shelf bathymetry.

function [XFlux,AFlux,fu,fv,uModeF,pModeF,vModeF,...
    up_modal_offshore,up_modal_all,up_modal_integrated,...
    xx,z,N2,h,hx,F] = ...
    func_simulation(Nx, Nz, R, Zpyc, mupyc, force_type, l,...
    L,h0,W,xs,hc,hs,f,sigma)


%INPUTS:
%     Nx[1,1]: number of grid points in x (cross-shore) direction
%     Nz[1,1]: number of grid points in z (depth) direction
%     R[1,1]: rigid lid or linear free surface
%         (R=0, rigid lid approx
%          R=1, linear free surface approx.)
%     Zpyc[1,1]: pycnocline depth (m) --> negative number!
%     mupyc[1,1]: pycnocline width (m)
%     force_type[1,1]: forcing type for the problem
%         (force_type=0, Baines body force 
%          force_type=1, isolated topographic force)
%     l[1,1]: along-shore wavenumber (1/m)
%     L[1,1]: cross-shore domain extent (m)
%     h0[1,1]: maximum ocean depth (m)
%     W[1,1]: width of topographic slope (m)
%     xs[1,1]: width of topographic shelf (m)
%     hc[1,1]: depth at the coastal boundary (m)
%     hs[1,1]: depth at the shelf break (m)
%     f[1,1]: Coriolis parameter (1/s)
%     sigma[1,1]: forcing frequency (1/s)

% OUTPUTS:
%     XFlux[Nx/2,1]: vertically-integrated cross-shore flux (W/m)
%     AFlux[Nx/2,1]: vertically-integrated along-shore flux (W/m)
%     fu[Nx/2,Nz]: cross-shore flux (W/m^2)
%     fv[Nx/2,Nz]: along-shore flux (W/m^2)
%     up_modal_offshore[Nz,Nz]: cross-shore energy flux (W/m^2) 
%         in each vertical mode 
%         at the off-shore boundary as a function of depth
%     up_modal_all[Nz,Nz]: vertically-integrated cross-shore energy flux
%               (W/m) in each mode and cross-mode
%     up_modal_integrated[Nz,1]: vertically-integrated cross-shore energy
%               flux (W/m) in each vertical mode at the off-shore boundary
%     xx[Nx,Nz]: x-coordinates (cross-shore)
%     z[Nx,Nz]: z-coordinates (depth
%     N2[Nx,Nz]: buoyancy frequency N^2 (1/s^2)
%     h[Nx,1]: depth at each cross-shore location (m)
%     hx[Nx,1]: horizontal gradient of depth, dh/dx
%     F[Nx,Nz]: forcing in pressure (odd x-points) and u (even x-points)
%     uModeF[Nx/2,Nz]: dimensionalized cross-shore velocity (m/s)
%     pModeF[Nx/2,Nz]: dimensionalized pressure (m^2/s^2)
%     vModeF[Nx/2,Nz]: dimensionalized along-shore velocity (m/s)

%%  Essentials
I = sqrt(-1) ;
eta_coast = 0.1;
                 
g = 9.81 ;
rho0 = 1000;

gamma = f/sigma;

ds = 1/(Nz-1) ;
s = (-1:ds:0)' ;

%%  Set up the grid

dx = L/Nx ;
x = dx*(0:1:Nx-1) ;
x = x-max(x) ;            %  odd grid points are for p equation, even grid points
                          %  are for u equation.  The grid ends at the
                          %  coastal wall at x = 0 ;
                          
np = 1:2:Nx ;
nu = 2:2:Nx ;
xx = ones(length(s),1)*x ;
                                                    
%  water depth and z grid.  The bottom depth and its derivative with respect to x
%  need to be specified.

h = h0+0*x;
for i=1:length(x)
   if x(i)>-xs-W
       h(i) = hs+(h0-hs)*(0.5*(1-cos((pi/W)*(x(i)+xs))))^0.75;
   end
end
[~,ind] = min(h);
h(ind:end) = linspace(h(ind),hc,length(h)-ind+1);
hx = gradient(h,dx);

z = 0*xx;
zx = 0*xx;
for ii = 1:Nx
    z(:,ii) = s*h(ii) ;
    zx(:,ii) = s*hx(ii);
end


%%  Stratification, buoyancy frequency and other related variables.  

%  stratification parameters
r1 = 1027 ;           %  density of upper layer
r2 = 1030 ;           %  density of lower layer
r0 = (r1+r2)/2 ;      %  prescribed rho0
dr = r2-r1 ;
%  Background linear density gradient
N2back = (2*pi/(0.5*60*60))^2 ;
rz0 = -N2back*r0/g ;

arg = (z - Zpyc)/mupyc ;
rho = r2 - 0.5*dr*(1 + tanh(arg)) ;
% 
if  sqrt(N2back) < 1e-6
    N2 = 0.5*(dr/r0)*(g/mupyc)*(sech(arg).^2) ;
    N2z_N2 = -(2/mupyc)*tanh(arg) ;
else
    rho = rho + rz0*(z + h0/2) ;
    N2 = 0.5*(dr/r0)*(g/mupyc)*(sech(arg).^2) + N2back ;
    N2z_N2 = -(dr/r0)*(g/(mupyc^2))*(sech(arg).^2).*tanh(arg)./N2 ;
end

%for non-hydrostatic term
N2z_N2is = N2z_N2.*N2./(N2-sigma^2);

%%
figure
clf
subplot(3,1,1)
plot(x/1000,-h) ;
axis([min(x)/1000 0 -h0 0]) ;
title('water depth (m) -- -h(x)')

subplot(3,1,2) 
plot(x/1000,-hx) ;
hold on
axis([min(x)/1000 0 0 max(-hx)]) ;
title('-dh/dx')

subplot(3,1,3) ;
plot(x/1000,-hx.*sqrt((N2(1,:) - f^2)/((sigma)^2 - f^2)),'r') ;
title('bottom slope / critical slope (>1 is supercritical)')
xlabel('cross-shore distance (km)')



figure
clf
subplot(2,1,1)
pcolor(x/1000,z,rho) ;
shading flat
colorbar
ylabel('density')

subplot(2,1,2)
pcolor(x/1000,z,sqrt(N2)) ;
shading flat
colorbar
ylabel('Depth (m)')
xlabel('Cross-shore distance (km)')
title('Buoyancy frequency N (s^{-1})')

%%  Here is the solution engine.


eyeNz = eye(Nz) ;
S = diag(s) ;  % s is the vertical coordinate
DS = ddz(s,1) ;  % vertical derivatives
DS2 = ddz2(s,1) ;

%  set up the forcing type
%
force = 0 ;        %  = 1 if the shelf is forced by an incident wave of mode p
                   %      the routine will check whether an incident
                   %      Poincare wave can exist with the specified
                   %      frequency, alongshore wavenumber and mode number.
                   %      If it can not, a warning is issued and the
                   %      program is stopped
                   %  = 0 for random forcing at the coast (only the u
                   %      equation is forced).
p = 1 ;            %  Mode number for incident wave.  Only used when force = 1.


% Calculate Baines body force
%first use barotropic SSH that matches Dale and Sherwin model
%and shelf width to calculate Baines' body force
eta0=0.1;
Fin = -eta0*(-xs)*zx.*N2./z;
Fin(end,:) = Fin(end-1,:);

%also need derivative of that F term (see manuscript)
dFdz = Fin*0;
for ix=1:Nx
    Z = h(ix)*s;
    DP = ddz(Z,0)/h(ix);
    dFdz(:,ix) = DP*Fin(:,ix);
end

% Now solve the coupled u-p equations at each (x,z) grid point
getMode;
get_initial;
scatter_invisc_solve;

%%  Post-processing the solution
% Dimensionalize cross-shore velocity and pressure, find along-shore velocity
[uModeF,vModeF,pModeF,XFlux,AFlux,fu,fv] = ...
normalize_pressure(PP,l,f,sigma,xx,z,h,np,nu,eta_coast);


%% Modal decomposition at the off-shore boundary
%  solve for the scattering vector (amplitude of reflected or evanescent
%  modes).
P1 = E1*E2inv ;
AAscat = E2inv*PP(:,2) ;
E1_amp = E1*diag(AAscat);

%normalize pressure such that SSH at the coast = eta_coast
pModeC = PP(:,np);
scl = max(abs(pModeC(:,end)));
E1_amp1 = eta_coast*E1_amp/scl;
%find cross-shore velocity (u) at the off-shore boundary
% from u-p relation
U1 = E1_amp1*((-I/(sigma*(1 - (gamma^2))))*diag((-I*k )));

%find cross-shore energy flux (up) for each mode at the off-shore boundary
%(W/m^2)
up_modal_offshore = g*rho0*(0.5*(real(E1_amp1).*real(U1) + ...
    imag(E1_amp1).*imag(U1)));

%find total cross-shore energy flux in u_i*p_j (includes modal energy i=j
% and cross-term modal energy i \neq j) (W/m^2)
up_modal_all = zeros(Nz-1,Nz-1);
for i=1:Nz-1
for j=1:Nz-1
up_modal_all(i,j) = h(1)*g*rho0*(mean(0.5*(real(E1_amp1(:,i)).*real(U1(:,j)) + ...
        imag(E1_amp1(:,i)).*imag(U1(:,j)))));
end
end

%vertically-integrated cross-shore energy flux at the off-shore boundary
%(W/m)
up_modal_integrated = h(1)*(mean(up_modal_offshore,1));
end
