function [out, Gamma_k, Gamma_kb] = makeF(r, y, omega, v, v_b, alpha, delta, NN, N)
k = omega*v;
k_b = omega*v_b;

Gamma_k = tools.fftQBiPerG(r, y, alpha, k, NN, N);
Gamma_kb = tools.fftQBiPerG(r, y, alpha, k_b, NN, N);
F1 = Gamma_k - Gamma_kb;
F2 = delta*tools.fftGradQBiPerG(r, y, alpha, k, NN, N) - tools.fftGradQBiPerG(r, y, alpha, k_b, NN, N);
out = [F1, F2].';

end