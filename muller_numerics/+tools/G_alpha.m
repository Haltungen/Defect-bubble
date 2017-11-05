%% G_alpha
%
% Overview:
%   Evaluates the quasi-biperiodic Greens function for the bubble crystal
%
% Input:
%   x,y:        Points of evaluation (vector, point)
%   omega:      Frequency
%   v,v_b:      Wave speed in medium and bubble, respectively
%   alpha:      Quasi-periodicity
%   delta:      Density fraction
%   R_b         Radius of bubble
%   NN:         Order of Fourier series
%   N1, N2:     Order of truncation for lattice sum and FFT, respectively 
%   N3:         Number of discretazation points for integration of S_D^\alpha 
%
% Output:
%   out:       The Green's function (Mx1 vector)

function out = G_alpha(x1,x2,y,omega,v,v_b,alpha,delta,R_b,NN,N1,N2,N3)
k = omega*v;
k_b = omega*v_b;

L = 1;
A = MakeA_Bubble(omega,v,v_b,alpha,delta,L,R_b,NN,N1);
F = makeF(R_b, y, omega, v, v_b, alpha, delta, NN, N2);

theta = linspace(0,2*pi*(N3-1)/N3,N3);
n = (-NN:NN).';

Phi = A\F;
phi_b = zeros(N3,1);
phi = zeros(N3,1);
for j = 1:N3
    phi_b(j) = sum(Phi(1:end/2).*exp(1i.*n.*theta(j)));
    phi(j) = sum(Phi(end/2+1:end).*exp(1i.*n.*theta(j)));
end

out = zeros(M,1);
for k = 1:M
    if sqrt(x1(k)*x1(k)+x2(k)*x2(k)) < R_b 
        out(k) = tools.GBiPeriodic(k_b, x1(k)-y(1), x2(k)-y(2), L, L, alpha) + ops.SBiPeriodic(k_b, alpha, R_b, phi_b, x1(k), x2(k));
    else
        out(k) = tools.GBiPeriodic(k, x1(k)-y(1), x2(k)-y(2), L, L, alpha) + ops.SBiPeriodic(k, alpha, R_b, phi, x1(k), x2(k));
    end
end