%% KstarBiPeriodic
%
% Overview:
%   Evaluates the quasi-biperiodic Neumann-Poincare operator K_D^*[F](x) of
%   a disk D.
%
% Input:
%   k:          Wave number
%   alpha:      Quasi-periodicity
%   r:          Radius of disk
%   F:          Layer density
%   x:          Point of evaluation  
%
% Output:
%   out:       The Neumann-Poincare operator

function out = KstarBiPeriodic(k, alpha, r, F, x1, x2)
    M = length(x1);
    N = length(F);
    theta = linspace(0,2*pi*(N-1)/N,N)';
    points1 = r*cos(theta);
    points2 = r*sin(theta);
    dsigma = r*2*pi/N;
    G = zeros(N,M);
    for l = 1:M
        G(:,l) = dot(tools.GradGBiPeriodic(k, points1-x1(l), points2-x2(l), 1, 1, alpha),[cos(theta),sin(theta)]);
    end
    out = (reshape(F, 1, []) .* dsigma) * G; % Dot product
end