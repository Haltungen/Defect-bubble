%% SBiPeriodic
%
% Overview:
%   Evaluates the quasi-biperiodic single layer potential S_D[F](x).
%
% Input:
%   k:          Wave number
%   alpha:      Quasi-periodicity
%   r:          Radius of disk
%   F:          Layer density
%   x:          Point of evaluation  
%
% Output:
%   out:       The Single layer potential

function out = SBiPeriodic(k, alpha, r, F, x1, x2)
    M = length(x1);
    N = length(F);
    theta = linspace(0,2*pi*(N-1)/N,N)';
    points1 = r*cos(theta);
    points2 = r*sin(theta);
    dsigma = r*2*pi/N;
    G = zeros(N,M);
    for l = 1:M
        G(:,l) = tools.GBiPeriodic(k, points1-x1(l), points2-x2(l), 1, 1, alpha);
    end
    out = (reshape(F, 1, []) .* dsigma) * G; % Dot product
end