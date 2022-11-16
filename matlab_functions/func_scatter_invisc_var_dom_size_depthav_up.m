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

function [XFlux,AFlux,AAscat,E1,P1,E1_amp1,U1, UP1,fu,fv,up,up_mode,xx,z,k,N2,uModeF,vModeF,pModeF,h] = ...
    func_scatter_invisc_var_dom_size_depthav_up(Nx, Nz, R,  ...
    Zpyc, mupyc, force_type, l,L,f)

%%  Essentials
I = sqrt(-1) ;
eta_coast = 0.1;

%  Horizontal and vertical grid size
%Nx = 200 ;
%Nz = 200 ;

%  Rigid lid or linear free surface
%R = 0 ;           %  R = 0, rigid-lid approx.
                  %  R = 1, linear free surface approx.

%  ocean and shelf dimensions
h0 = 3100 ;       
%L = 300e3;
W = 32e3;
xs = 80e3;
hc = 100;
hs = 150;
hmatch = 40;
xmatch = 13e3;
                 
%f = 9.3e-5 ; %2.91e-4;%2.18e-4;%9.3e-5 ;
g = 9.81 ;

sigma = 1.407e-4; %1.25*f;
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

% figure
% clf
% subplot(2,1,1)
% plot(x/1000,-h) ;
% axis([min(x)/1000 0 -h0 0]) ;
% title('water depth (m) -- -h(x)')
% 
% subplot(2,1,2) 
% plot(x/1000,-hx) ;
% hold on
% axis([min(x)/1000 0 0 max(-hx)]) ;
% title('-dh/dx')

% subplot(4,1,3) ;
% N2bot = 0.5*(dr/r0)*(g/mupyc)*(sech((-h - Zpyc)/mupyc).^2) + N2back ;
% plot(x/1000,sqrt(2*pi./N2bot)/60,'r') ;
% title('buoyancy-frequency squared (N^2) along the bottom')
% 
% subplot(4,1,4) ;
% plot(x/1000,-hx.*sqrt((N2bot - f^2)/((1.1*f)^2 - f^2)),'r') ;
% title('bottom slope / critical slope for sigma = 1.1f (>1 is supercritical)')
%xlabel('cross-shore distance (km)')

%%  Stratification, buoyancy frequency and other related variables.  
% here buoyancy frequence N2 = constant
% subsequently, dN2_dz/N2 = 0

%  stratification parameters
% depth and thickness of pycnocline are function parameters
%Zpyc = -600 ;         %  depth of pycnocline (m)
%mupyc = 300 ;         %  thickness of pycnocline (m)
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
%N2 = 0*z +(2.72e-3)^2;
%N2z_N2 = 0*N2;

%for non-hydrostatic term
N2z_N2is = N2z_N2.*N2./(N2-sigma^2);

% figure
% clf
% 
% % subplot(2,1,1)
% % pcolor(x/1000,z,rho) ;
% % shading flat
% % colorbar
% % ylabel('density')
% 
% %subplot(2,1,2)
% pcolor(x/1000,z,sqrt(N2)) ;
% shading flat
% colorbar
% ylabel('Depth (m)')
% xlabel('Cross-shore distance (km)')
% title('Buoyancy frequency N (s^{-1})')

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


%  In this example, we set the frequency to a super-inertial value (1.25*f)
%  and search for super-intertial coastal trapped waves at that frequency.
%sigma = 0.75*f ;

% Calculate Baines body force
%first use barotropic SSH that matches Dale and Sherwin model
%and shelf width to calculate Baines' body force
eta0=0.1;
%Fin = -eta0*(-L)*zx.*N2./z;
Fin = -eta0*(-xs)*zx.*N2./z;
Fin(end,:) = Fin(end-1,:);

%also need derivative of that F term (see write-up)
dFdz = Fin*0;
for ix=1:Nx
    Z = h(ix)*s;
    DP = ddz(Z,0)/h(ix);
    dFdz(:,ix) = DP*Fin(:,ix);
end


%sigma = 1.25*f; gamma = f/sigma ;
%l = 6.2062e-6;
%force_type = 3;
getMode;
get_initial;
scatter_invisc_solve;

[uModeF,vModeF,pModeF,XFlux,AFlux,fu,fv] = ...
normalize_pressure(PP,l,f,sigma,xx,z,h,np,nu,eta_coast);

P1 = E1*E2inv ;
AAscat = E2inv*PP(:,2) ;
E1_amp = E1*diag(AAscat);

pModeC = PP(:,np);%*rho0;
scl = max(abs(pModeC(:,end)));

E1_amp1 = eta_coast*E1_amp/scl;
U1 = E1_amp1*((-I/(sigma*(1 - (gamma^2))))*diag((-I*k )));

UP1 = 0.5*(real(E1_amp1).*real(U1) + ...
    imag(E1_amp1).*imag(U1));

up = zeros(Nz-1,Nz-1);
for i=1:Nz-1
for j=1:Nz-1
up(i,j) = mean(0.5*(real(E1_amp1(:,i)).*real(U1(:,j)) + ...
imag(E1_amp1(:,i)).*imag(U1(:,j))));
end
end

up_mode = mean(UP1,1);
end
