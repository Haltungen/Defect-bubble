function out = makeRHS(k0, k_b, delta, y, B_d)

green = @(x1,x2,k) -1i/4*besselh(0,k*sqrt((x1-y(1)).^2+(x2-y(2)).^2));
green_nu = @(x1,x2,k) -1i/8*k*(besselh(-1,k*sqrt((x1-y(1)).^2+(x2-y(2)).^2))-besselh(1,k*sqrt((x1-y(1)).^2+(x2-y(2)).^2)));
xd1 = B_d.points(1,:);
xd2 = B_d.points(2,:);

R1 = -green(xd1,xd2,k_b) + green(xd1,xd2,k0);
R2 = zeros(size(xd1));
R3 = delta*green_nu(xd1,xd2,k0)-green_nu(xd1,xd2,k_b);
R4 = zeros(size(xd1));

out = [R1.';R2.';R3.';R4.'];