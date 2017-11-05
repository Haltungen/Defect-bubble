%% G_alpha
%
% Overview:
%   Evaluates the quasi-biperiodic Greens function for the bubble crystal
%
% Input:
%   x,y:        Points of evaluation (vector, point)
%   k,k_b:      Wave numbers in medium and bubble, respectively
%   alpha:      Quasi-periodicity
%   delta:      Density fraction
%   R_b         Radius of bubble
%   NN:         Order of Fourier series
%   N1, N2:     Order of truncation for lattice sum and FFT, respectively          
%
% Output:
%   out:       The Single layer potential

function out = G_alpha(x,y,k,k_b,alpha,delta,R_b,NN,N1,N2)

L = 1;
A = MakeA_Bubble(k,k_b,alpha,delta,L,R_b,NN,N1);
F = makeF(R_b, y, k, k_b, alpha, delta, NN, N2);

M = lenght(x);
theta = linspace(0,2*pi*(M-1)/M);
n = -NN:NN;

Phi = A\F;
phi_b = zeros(M,1);
phi = zeros(M,1);
for j = 1:M
    phi_b(j) = sum(Phi(1:end/2).*exp(1i.*n.*theta(j)));
    phi(j) = sum(Phi(end/2:end).*exp(1i.*n.*theta(j)));
end

out = zeros(M,1);
for k = 1:M
    if abs(x(k)) < R_b 
        out(k) = tools.GBiPeriodic(k_b, x(k), y, L, L, alpha) + SBiPeriodic(k_b, alpha, D, phi_b, x);
    else
        out(k) = tools.GBiPeriodic(k, x(k), y, L, L, alpha) + SBiPeriodic(k, alpha, D, phi, x);
    end
end



