%  scatter_invisc_solve.m
%
%  Jim Lerczak
%  College of Earth, Ocean, and Atmospheric Sciences
%  Oregon State University
%  jlerczak@coas.oregonstate.edu
%  November 2016



%  Now iterate forward to get all alpha's and beta's
for ix = 2:Nx ;
    getABC ; %  Get the A, B and C matrices for the viscid scattering problem:
             %  A(i)P(i-1) + B(i)P(i) + C(i)P(i+1) + F(i) = 0
    MM = A*al{ix-1} + B ;
    al{ix} = MM\(-C) ;
    bt{ix} = MM\(F(:,ix) - A*bt{ix-1}) ;
end

PP = zeros(Nz,Nx) ;
ix = Nx-1 ;
% getABC ;        
% MM = A*al{ix-1} + B ;
% PP(:,ix) = MM\(F(:,ix) - A*bt{ix-1}) ;
PP(:,ix) = bt{ix} ;

for ix = Nx-2:-1:1 ;
    PP(:,ix) = al{ix}*PP(:,ix+1) + bt{ix} ;
end

%  solve for the scattering vector (amplitude of reflected or evanescent
%  modes).
if force == 0
    AAscat = E2inv*PP(:,2) ;
else
    AAscat = E2inv*PP(:,2) - E2inv*(a*Pp) ;
end
