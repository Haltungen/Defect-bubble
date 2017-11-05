clear all
close all

rho = 5000;
rho_b = 1;
kappa = 5000;
kappa_b = 1;
omega = 0.2;

delta = rho_b/rho;
v = sqrt(rho/kappa);
v_b = sqrt(rho_b/kappa_b);

k0 = omega * v;
kb = omega * v_b;

NN = 3;     % Order of fourier series
N1 = 5;     % Order of truncation for lattice sum
N2 = 3;     % Order of trucation for FFT   
N3 = 10;    % Number of discretization points in quadrature for integrating G_alpha
R_b = 0.5;

xMax = 5;
nX = 2^4;
hX = 2*xMax/(nX-1);
[x1, x2] = meshgrid(-xMax:hX:xMax);
r_2 = 0.2;
theta_2 = 0;
y = r_2*[cos(theta_2),sin(theta_2)];
Green = zeros(nX);
for l = 1:nX
    Green(:,l) = G(x1(:,l),x2(:,l),y,omega,v,v_b,delta,R_b,NN,N1,N2,N3);
end

figure('Position', [0 0 1600 800]);
surf(x1, x2, real(G), 'edgecolor', 'none');
xlabel('x_1'); ylabel('x_2'); zlabel('G')
title('Re(G)')
view([0 0 90]); rotate3d on; colormap jet;


