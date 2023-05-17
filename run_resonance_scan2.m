% File to accompany manuscript by Zemskova, V.E., Musgrave, R.C. and
% Lerczak, J. A., 
%titled "Internal tides at the coast: energy flux of baroclinic tides propagating into
%the deep ocean in the presence of supercritical shelf topography"

%File to run a resonance scan over a range of frequency and along-shore
%       wavenumbers with specified grid spacing, topography,
%       and stratification parameters 

%   Also can specify Coriolis parameter, forcing type (Baines vs. isolated), 
%       surface boundary
%       condition (rigid lid vs. free linear surface)
%        
%   Outputs domain-integrated values of pressure and cross-shore velocity
%       (u) for each value of frequency and along-shore wavenumber.
%       Modes are resonant if these values are amplified.

%   Plots resonance scan (domain-integrated pressure values) sweeping over 
%       frequency and along-shore wavenumbers and dispersion curves for
%       other relevant modes.

%  Barbara Zemskova
%  College of Earth, Ocean, and Atmospheric Sciences
%  Oregon State University
%  barbara.zemskova@oregonstate.edu
%  May, 2023

%% Set-up parameters

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

%% Range of frequencies and along-shore wavenumbers to consider
Nsigma = 1000;
Nl = 400;

SI = linspace(0.001*f, 2*f,Nsigma); %range of frequencies
LI = linspace(0,1e-4,Nl); %range of along-shore wavenumbers

%% Run resonance scan

%domain-integrated cross-shore velocity response
P0u_sweep = zeros(Nsigma,Nl);

%domain-integrated pressure response
P0p_sweep = zeros(Nsigma,Nl);



for is = 1:Nsigma
    disp(is)
    for il = 1:Nl
            [P0p,P0u] = ...
                 func_resonance(Nx, Nz, R, Zpyc, mupyc, force_type, l,...
                L,h0,xW,xs,hc,hs,f,sigma,N2back,eta0);
            %  cross-shore current variance across the domain in response to
            %  forcing
            P0u_sweep(il) = P0u ;
            %  pressure variance across the domain in response to forcing
            P0p_sweep(il) = P0p ;
    end
    save('resonance_scan_600_200_R_0.mat')
end

%% Plot the resulting resonance scan
% values at the resonant along-shore wavenumbers will have amplified
% response


figure;
pcolor(LI,SI/f,log10(abs(P0p_sweep))); shading flat
cmap = colormap(gray) ; 
colormap(flipud(cmap)) ;
colorbar

hold on
ds = 2*f/1000 ;
ss = f+ds/2:ds:max(SI) ;
ll = sqrt((ss.^2 - f^2)/(g*h0)) ;
plot(ll,ss/f,'b','linewidth',1)
plot(-ll,ss/f,'b','linewidth',1)
ss = 0:ds:max(SI) ;
ll = ss/sqrt(g*h0) ;
plot(ll,ss/f,'r--','linewidth',1) ;
ll = ss/sqrt(g*hs) ;
plot(ll,ss/f,'g--','linewidth',1) ;
plot(LI,1 + 0*LL,'k--') ;

ss = f+ds/2:ds:max(SI) ;
for ii = 1:4
    scl = 4/(2*ii - 1) ;
    kk = 2*pi/(scl*Lsh) ;
    ll = sqrt((ss.^2 - f^2)/(g*Hs) - kk^2) ;
    mm = find(real(ll)>0) ;
    plot(ll(mm),ss(mm)/f,'m--') ;
    plot(-ll(mm),ss(mm)/f,'m--') ;
end


ylabel('\sigma/f','fontsize',16) ;
xlabel('Alongshore wavenumber (m^{-1})','FontSize',16) ;