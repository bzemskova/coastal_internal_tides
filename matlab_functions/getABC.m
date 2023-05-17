%  getABC.m
%
%  Jim Lerczak
%  College of Earth, Ocean, and Atmospheric Sciences
%  Oregon State University
%  jlerczak@coas.oregonstate.edu
%  November 2016
%
%  Get the A, B and C matrices for the viscid scattering problem
%
%   A(i)P(i-1) + B(i)P(i) + C(i)P(i+1) + F(i) = 0
%

dz = ds*h(ix) ;
DX = 2*dx ;
Z = h(ix)*s ;

%  for odd ix, solve for the p-equation.  ix+1 and ix-1 are the u values
%  surrounding p on the staggered grid
if mod(ix,2) == 1                           %  odd ix: p equation
    
    %N2_is = N2(:,ix)/(I*sigma) ;
    N2_is = (N2(:,ix)-sigma^2)/(I*sigma) ; %replaced with (N2-sigma^2)/(i sigma)
%     N2zN2 = N2z_N2(:,ix) ;
    N2zN2 = N2z_N2is(:,ix) ;  %replaced with (dN2/dz)/(N2-sigma^2)
    
    du = ddz(Z,1) ;
    
    A =     -diag(N2_is)/DX ...
          -  (1/2)*diag(hx(ix)*N2_is.*s)*du ...
          +  (1/2)*diag(N2_is)*gamma*l ;
    %  add terms to A for the bottom BC of pressure
    A(1,1) = A(1,1) + (1/2)*N2zN2(1)*N2_is(1)*hx(ix) ...
                    + (1/dz)*N2_is(1)*hx(ix) ;
    
    C =      diag(N2_is)/DX ...
          -  (1/2)*diag(hx(ix)*N2_is.*s)*du ...
          +  (1/2)*diag(N2_is)*gamma*l ;
    %  add terms to A for the bottom BC of pressure
    C(1,1) = C(1,1) + (1/2)*N2zN2(1)*N2_is(1)*hx(ix) ...
                    + (1/dz)*N2_is(1)*hx(ix) ;
        
    %  first derivative for pressure with appropriate boundary conditions
    DP = ddz(Z,0);%/h(ix) ;
    DP(Nz,Nz-4:Nz) = [0 0 0 0 -R*N2(Nz,ix)/g];
    DP(1,1:5) = [0 0 0 0 0];
    
    %  second derivative for pressure with appropriate boundary conditions
    DP2 = ddz2(Z,0) ;
    DP2(Nz,Nz-4:Nz) = [0 0 0 2 -2*(1 + R*dz*N2(Nz,ix)/g)]/(dz^2) ;
    DP2(1,1:5) = [-2 2 0 0 0]/(dz^2) ;
    
    B =   diag(((N2(:,ix)-sigma^2)/(sigma^2))*(l^2)) ... %replaced with (N2-sigma^2)/sigma^2
        - diag(N2zN2)*DP                   ...
        + DP2 ;
    %diag((N2(:,ix)/(sigma^2))*(l^2)) ...

else                           %  even ix (u equation)
    
    %  first derivative for pressure with appropriate boundary conditions
    DP = ddz(Z,0) ;
    DP(Nz,Nz-4:Nz) = [0 0 0 0 -R*N2(Nz,ix)/g] ;
    DP(1,1:5) = [0 0 0 0 0] ;
    
    A =  -(1/DX)*eyeNz ...
        - (1/2)*diag(hx(ix)*s)*DP ...
        - (1/2)*gamma*l*eyeNz ;
    
    C =   (1/DX)*eyeNz ...
        - (1/2)*diag(hx(ix)*s)*DP ...
        - (1/2)*gamma*l*eyeNz ;
    
    B = ((1 - gamma^2)*sigma/I)*eyeNz ;
    
%     B(1,1) = B(1,1) - (hx(ix)^2)*N2(1,ix)/(I*sigma) ;
    B(1,1) = B(1,1) - (hx(ix)^2)*(N2(1,ix)-sigma^2)/(I*sigma) ;
            %replaced with (N2-sigma^2)
    
end

return
   