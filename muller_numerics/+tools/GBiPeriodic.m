%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GBiPeriodic
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Remarks:
%   The quasi-biperiodic Green's function for the Helmholtz
%   equation.
%
% References:
%   Mathematical and Computational Methods in Photonics - Tutorial Notes
%   Mathematical and Computational Methods in Photonics - Lecture Notes
%
% Authors:
%   Habib Ammari, Brian Fitzpatrick, Matias Ruiz, Sanghyeon Yu.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ G GSpec GSpat ] = GBiPeriodic(k, x1, x2, d, h, alpha)

E = pi/d/h;

%%%%%%%%%%%%%%%%%%%%%%%%%
% GSpec
%%%%%%%%%%%%%%%%%%%%%%%%%
pqMax = 5;

GSpec = zeros(size(x1));
parfor p=-pqMax:pqMax
    for q=-pqMax:pqMax
        kxp = -alpha(1) + 2*pi*p/d;
        kyq = -alpha(2) + 2*pi*q/h;
        
        kpqSqr = kxp.^2 + kyq.^2;
        gammapqSqr = kpqSqr - k^2;
        
        GSpec = GSpec + exp(-gammapqSqr/4/E).*exp(-1i*(kxp*x1 + kyq*x2))./gammapqSqr;
    end
end

GSpec = -1/d/h * GSpec;

%%%%%%%%%%%%%%%%%%%%%%%%%
% GSpat
%%%%%%%%%%%%%%%%%%%%%%%%%
mnMax = 5;
Q = 10;

GSpat = 0;
parfor m=-mnMax:mnMax
    for n=-mnMax:mnMax
        rhomn = [ m*d; n*h ];
        rho_rhomnSqr = (x1-rhomn(1)).^2 + (x2-rhomn(2)).^2;
        
        term1 = 0;
        Eqp1 = tools.E1(rho_rhomnSqr*E);
        for q=0:Q
            term1 = term1 + (k/2/sqrt(E)).^(2*q)./factorial(q).*Eqp1;
            Eqp1 = (exp(-rho_rhomnSqr*E) - rho_rhomnSqr*E.*Eqp1)./(q+1);
        end
        
        GSpat = GSpat - 1/4/pi.*exp(1i*(alpha(1)*rhomn(1) + alpha(2)*rhomn(2))).*term1;
    end
end

G = GSpat + GSpec;

end