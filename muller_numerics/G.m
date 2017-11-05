%% G_alpha
%
% Overview:
%   Evaluates the quasi-biperiodic Greens function for the bubble crystal
%
% Input:
%   x,y:        Points of evaluation (vector, point)
%   omega:      Frequency
%   v,v_b:      Wave speed in medium and bubble, respectively
%   delta:      Density fraction
%   R_b         Radius of bubble
%   NN:         Order of Fourier series
%   N1, N2:     Order of truncation for lattice sum and FFT, respectively   
%   N3:         Number of discretization points in quadrature
%
% Output:
%   out:       The Single layer potential

function out = G(x1,x2,y,omega,v,v_b,delta,R_b,NN,N1,N2,N3)
M = length(x1);
out = zeros(M,1);
alpha1D = linspace(0,2*pi*(N3-1)/N3,N3);
d2alpha = (2*pi/N3)^2;
for j = 1:N3 % Trapezoid quadrature method
    for l = 1:N3
        alpha = [alpha1D(j),alpha1D(l)];
        out = out + tools.G_alpha(x1,x2,y,omega,v,v_b,alpha,delta,R_b,NN,N1,N2,N3)*d2alpha;
    end
end

end