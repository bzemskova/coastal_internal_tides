function [c,P,W,ML,MR]=nmodesP_f_01(zin,N2in,N2z_N2in,sigma,f,R)
%
% USAGE: [c,P,W]=nmodesP_f_01(zin,N2in,N2z_N2in,sigma,f,R)
%
% Computes baroclinic modes for arbitrary N2 profile using the vertical
% pressure structure function equation
% INPUTS:
% z = vertical coordinate vector (evenly spaced, and starting at the bottom
% and ending at the surface)
% N2 = squared buoyancy frequency
% N2_N2z = (1/N2)(dN2/dz) evaluated at the grid points
% sigma = frequency of wave
% f = Coriolis frequency
% R = 0 for rigid lid
%   = 1 for linear, free-surface condition
% 
% OUTPUTS:
% c = phase speed (real or imaginary)
% P = pressure eigenfunctions
% W = vertical velocity eigenfunction
%
%  Jim Lerczak
%  College of Earth, Ocean, and Atmospheric Sciences
%  Oregon State University
%  jlerczak@coas.oregonstate.edu
%  November 2016
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
g = 9.81 ;

if abs(sigma-f)<1e-12 
    disp('Can not solve at the inertial frequency') ;
    c = [] ;
    W = [] ;
    P = [] ;
    return
end

z = zin ;
N2 = N2in ;
N2z_N2 = N2z_N2in ;

dz = nanmean(diff(z)) ;

M = length(z) ;

% Assemble matrices
D2=ddz2(z,0);                     %2nd derivative matrix for rigid lid condition
                                  %Correct the matrix for the surface and
                                  %bottom BCs
%  zero gradient at the bottom
D2(1,1:5) = [-2 2 0 0 0]/(dz^2) ;
%  zero gradient or linear surface condition at surface
D2(end,end-4:end) = [0 0 0 2 -2*(1 + R*dz*N2(end)/g)]/(dz^2) ;

D1=ddz(z,0);                     %1st order derivative matrix for rigid lid condition
                                  %Correct the matrix for the surface and
                                  %bottom BCs
%  zero gradient at the bottom
D1(1,1:5) = [0 0 0 0 0] ;
%  zero gradient or linear surface condition at surface
D1(end,end-4:end) = [0 0 0 0 -R*N2(end)/g] ;

ML = D2 - diag(N2z_N2)*D1 ;
% MR = -diag(N2)/(1 - (f/sigma)^2) ;
MR = -diag(N2-sigma^2)/(1 - (f/sigma)^2) ; %replaced with (N2-sigma^2)

% Solve generalized eigenvalue problem
[v,e]=eig(ML,MR);
 
% Convert eigvals to wave speeds for superinertial  frequencies
if sigma>f
    cin=diag(e);
    cs=sqrt(1./cin);
 
% Sort speeds and eigvecs, but be sure to keep the full eigenvalue as
% determined from the eigenvalue problem
    [val,ind]=sort(real(cs),1,'descend');
    c = cs(ind) ;
    c = c(1:M) ;
    
else  %  subinertial
 
    % Look for negative eigenvalues
    c2=1./diag(e);
    [vl,ind] = sort(-c2,'descend') ;
    % Sort eigenvalues and eigvecs
    c2 = c2(ind) ;
    c = sqrt(c2(1:M)) ;
    
end

p=v(:,ind);
P = zeros(length(zin),M) ;
W = zeros(length(zin),M) ;
for j=1:M
    ap=abs(p(:,j));
    con=max(p(abs(p(:,j))==max(ap),j));
    P(:,j)=p(:,j)/con;
    W(:,j)=(1./N2).*(D1*P(:,j)) ; 
    aw= abs(W(:,j)) ;
    con=max(W(abs(W(:,j))==max(aw),j));
    W(:,j) = W(:,j)/con ;
end

return
