clear all
close all

rho = 5000;
rho_b = 1;
kappa = 5000;
kappa_b = 1;
omega = 0.28;

delta = rho_b/rho;
v = sqrt(rho/kappa);
v_b = sqrt(rho_b/kappa_b);

k0 = omega * v;
kb = omega * v_b;

NN = 2;     % Order of fourier series
N1 = 1;     % Order of truncation for lattice sum
N2 = 1;     % Order of trucation inside Fourier coefficients  
N3 = 20;    % Number of discretization points for integration of S_D^\alpha 
R_b = 0.05;

xMax = 0.06;
nX = 10^2;
hX = xMax/(nX-1);
x1 = 0.04:hX:xMax;
r_2 = 0;
theta_2 = 0;
y = r_2*[cos(theta_2),sin(theta_2)];
Green = tools.G_alpha(x1,zeros(1,nX),y,omega,v,v_b,[pi,pi],delta,R_b,NN,N1,N2,N3);
figure
plot(x1,real(Green),'.');
figure
plot(x1,imag(Green),'.');
