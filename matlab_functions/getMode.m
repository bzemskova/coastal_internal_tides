%  getMode
%
%  Jim Lerczak
%  College of Earth, Ocean, and Atmospheric Sciences
%  Oregon State University
%  jlerczak@coas.oregonstate.edu
%  November 2016

%  get the vertical structure at the deep ocean
% [c,P,W]=nmodesP_f_01(z(:,1),N2(:,1),N2z_N2(:,1),sigma,f,R) ;
[c,P,W,ML,MR]=nmodesP_f_01(z(:,1),N2(:,1),N2z_N2is(:,1),sigma,f,R) ;
        %replaced with (dN2/dz)/(N2-sigma^2)

%  When R = 0, (rigid lid) the eigenvalue analysis sometimes picks out a barotropic
%  mode with 'infinite' phase speed.  I think this happens when the
%  verticial derivative of the buoyancy frequency is non-zero at the
%  surface or bottom.  I will remove this mode.

if (R==0)
    if (abs(c(1)) > 100)
        P = P(:,2:end) ;
        W = W(:,2:end) ;
        c = c(2:end) ;
    end
    
    if abs(c(end))>100
         P = P(:,1:end-1) ;
         W = W(:,1:end-1) ;
         c = c(1:end-1) ;
      end
  
end


        
%  number of modes allowed at the offshore end
M = length(c); 
P = P(:,1:M) ;
W = W(:,1:M) ;
c = c(1:M) ;

