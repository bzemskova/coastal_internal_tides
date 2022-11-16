%  Get the forcing matrix F and alpha(1) and beta(1) for the calculation
%
%  Jim Lerczak
%  College of Earth, Ocean, and Atmospheric Sciences
%  Oregon State University
%  jlerczak@coas.oregonstate.edu
%  November 2016

%  Get the across-shore length scales for each vertical mode
K2 = (sigma^2)./(c.^2) ;
k = sqrt(K2 - l^2) ;

E1 = P ; % P are the pressure vertical modes
E2 = P*((-I/(sigma*(1 - (gamma^2))))*diag((-I*k - gamma*l).*exp(-I*k*dx))) ;

E2inv = ((E2')*E2)\(E2') ;

F = 0*N2 ;


if force == 0
    if force_type==0
        %Baines forcing term only applies to p-equation
        F(:,np) = -((dFdz(:,np)-...
        Fin(:,np).*N2z_N2(:,np).* N2(:,np)./(N2(:,np)-sigma^2))) ;
  
    elseif force_type==1
        %isolated topographic forcing, applied only to p-equation
            FB(:,np) = -((dFdz(:,np)-...
        Fin(:,np).*N2z_N2(:,np).* N2(:,np)./(N2(:,np)-sigma^2))) ;
        for i=1:length(np)   
            if abs(x(np(i))+100e3)<1e3
                F(:,np(i))=max(abs(FB(:))); 
            end
        end
    end
 
    
    al{1} = E1*E2inv ;
    bt{1} = zeros(Nz,1) ;
    
else
    %  
    if R == 0
        mode = p+1 ;
    else
        mode = p+1 ;

    end
    if ~isreal(k(mode))
        disp('force = 1') ;
        disp('incident mode is evanescent') ;
    else
        Pp = 0*P(:,1);
        a = 0;

        for ii = 1:mode
            Pp = Pp + P(:,ii) ;
            a = a + (-I/(sigma*(1 - (gamma^2))))*(I*k(ii) - gamma*l)*exp(I*k(ii)*dx) ;
        end
        %        Pp = P(:,mode);
        %        a = a+ (-I/(sigma*(1 - (gamma^2))))*(I*k(mode) - gamma*l)*exp(I*k(mode)*dx) ;
        %        Pp = sum(P(:,1:mode),2); %FORCE With many modes, but a needs to be modified...

        
        al{1} = E1*E2inv ;
        bt{1} = (eyeNz - (E1*E2inv)*a)*Pp ;
        
    end
        
end
