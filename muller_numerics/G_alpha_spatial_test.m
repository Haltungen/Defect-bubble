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

nPoints = 1000;
R_b = 0.05;

xMax = 0.051;
xMin = 0.049;
nX = 10;
hX = (xMax-xMin)/(nX-1);
x1 = xMin:hX:xMax;
r_2 = 0;
theta_2 = 0;
y = r_2*[cos(theta_2),sin(theta_2)];
B = shape.Ellipse(R_b, R_b, nPoints);
Green = ops.G_alpha_spatial(x1,zeros(1,nX),y,omega,v,v_b,[pi,pi],delta,B);
figure
plot(x1,real(Green),'.');
figure
plot(x1,imag(Green),'.');