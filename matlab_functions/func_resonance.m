%  Barbara Zemskova
%  College of Earth, Ocean, and Atmospheric Sciences
%  Oregon State University
%  barbara.zemskova@oregonstate.edu
%  May, 2023
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
%     N2back[1,1]: background linear stratification N^2 (1/s^2)
%     eta0[1,1]: sea surface elevation at the coast (m)

% OUTPUTS:
%     P0u_sd[1,1]: domain (x,z)-integrated velocity response
%     P0p_sd[1,1]: domain (x,z)-integrated pressure response


function [P0u_sd,P0p_sd] = ...
                func_resonance(Nx, Nz, R, Zpyc, mupyc, force_type, l,...
                L,h0,xW,xs,hc,hs,f,sigma,N2back,eta0)

%%  Essentials
I = sqrt(-1) ;
                 
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
   if x(i)>-xs-xW
       h(i) = hs+(h0-hs)*(0.5*(1-cos((pi/xW)*(x(i)+xs))))^0.75;
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
rz0 = -N2back*r0/g ;

arg = (z - Zpyc)/mupyc ;
rho = r2 - 0.5*dr*(1 + tanh(arg)) ;

if  sqrt(N2back) < 1e-6
    N2 = 0.5*(dr/r0)*(g/mupyc)*(sech(arg).^2) ;
    N2z_N2 = -(2/mupyc)*tanh(arg) ;
else
    rho = rho + rz0*(z + h0/2) ;
    N2 = 0.5*(dr/r0)*(g/mupyc)*(sech(arg).^2) + N2back ;
    N2z_N2 = -(dr/r0)*(g/(mupyc^2))*(sech(arg).^2).*tanh(arg)./N2 ;
end

% Constant linear stratification
if Zpyc==0 && mupyc==0
    N2 = N2back*ones(size(xx));
    N2z_N2 = 0*N2;
end

% For non-hydrostatic term
N2z_N2is = N2z_N2.*N2./(N2-sigma^2);

% Topographic criticality parameter
alpha = -hx.*sqrt((N2(1,:) - f^2)/((sigma)^2 - f^2));


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
Fin = 0*xx;
for i=1:length(x)
    Fin(:,i) = -eta0*(-xs)*hx(i).*N2(:,i).*z(:,i)/(h(i)^2);
end

%also need derivative of that F term (see manuscript)
dFdz = Fin*0;
for ix=1:Nx
    Z = h(ix)*s;
    DP = ddz(Z,0);%/h(ix);
    dFdz(:,ix) = DP*Fin(:,ix);
end

% Now solve the coupled u-p equations at each (x,z) grid point
getMode;
get_initial;
scatter_invisc_solve;

%  cross-shore current variance across the domain in response to
    %  forcing
Ar = 0 ;
vr = 0 ;
for iu = nu ;
    Ar = Ar + 2*dx*h(iu) ;
    vr = vr + h(iu)*ds*2*dx*(PP(:,iu)'*PP(:,iu)) ;
end
P0u_sd = vr/Ar ;
%  pressure variance across the domain in response to forcing
Ar = 0 ;
vr = 0 ;
for ip = np ;
    Ar = Ar + 2*dx*h(ip) ;
    vr = vr + h(ip)*ds*2*dx*(PP(:,ip)'*PP(:,ip)) ;
end
P0p_sd = vr/Ar ;

end