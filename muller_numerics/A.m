function out = A(k, k_b, delta, B, B_d)
% A represents the matrix A(omega, delta) in equation (2.3) of Ref 1.
%
% Ref 1: Minnaert resonance for acoustic waves in bubbly media

nPoints = length(B.points);

A11 = ops.SingleLayer_H(k_b, B_d, 'P0', 1);
A12 = ops.SingleLayer_H(k, B_d, 'P0', 1);
A13 = ops.SingleLayer_H(k, B_d, 'P0', 1, B, 'P0', 1);    
A22 = ops.SingleLayer_H(k, B, 'P0', 1, B_d, 'P0', 1);
A23 = ops.SingleLayer_H(k, B, 'P0', 1);
A24 = ops.SingleLayer_H_C(k, k_b, B, 'P0', 1);
A31 = ops.Kstar_H(k_b, B_d, 'P0', 1);
A32 = ops.Kstar_H(k, B, 'P0', 1);
A33 = ops.Kstar_H(k, B_d, 'P0', 1, B, 'P0', 1);
A42 = ops.Kstar_H(k, B, 'P0', 1, B_d, 'P0', 1);
A43 = ops.Kstar_H(k, B, 'P0', 1);
A44 = ops.Kstar_H_C(k, k_b, B, 'P0', 1);


out = [ A11.Kmat, -1*A12.Kmat, -1*A13.Kmat, zeros(nPoints); ...
        zeros(nPoints), A22.Kmat, A23.Kmat, -1*A24.Kmat; ...
        -1/2.*eye(nPoints) + A31.Kmat, -delta*(1/2.*eye(nPoints) + A32.Kmat), -delta*A33.Kmat, zeros(nPoints); ...
        zeros(nPoints), A42.Kmat, -1/2.*eye(nPoints) + A43.Kmat, -1*(1/2.*eye(nPoints)+ A44.Kmat)];

end
