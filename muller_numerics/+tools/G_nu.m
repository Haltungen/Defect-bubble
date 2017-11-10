%% G_nu
%
% Overview:
%   Evaluates the radial derivative of the Greens function for the bubble crystal
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

function out = G_nu(x1,x2,y,omega,v,v_b,delta,R_b,NN,N1,N2,N3)
if nargin < 12 % Default parameters
    NN = 3;
    N1 = 3;
    N2 = 6;
    N3 = 10;
end
M = length(x1);
out = zeros(M,1);
alpha1D = linspace(0,2*pi*(N3-1)/N3,N3);
d2alpha = (2*pi/N3)^2;
for j = 1:N3 % Trapezoid quadrature method
    for l = 1:N3
        alpha = [alpha1D(j),alpha1D(l)];
        out = out + tools.G_alpha_nu(x1,x2,y,omega,v,v_b,alpha,delta,R_b,NN,N1,N2,N3)*d2alpha;
    end
end

end