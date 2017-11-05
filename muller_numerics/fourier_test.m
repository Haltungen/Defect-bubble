clear all
close all

theta = 0:0.25:2*pi;
M = length(theta);
G_fourier = zeros(M,1);
G_nu = zeros(M,1);
G_nu_fft = zeros(M,1);
G_fourier_fft = zeros(M,1);
NN = 5;
m = 50;
n = -NN:NN;
N = 10*NN*m;
R_b = 0.5;
r_p = 0;
theta_p = 0;
y = r_p*[cos(theta_p), sin(theta_p)];
alpha = [pi/2, pi/8];
k = 1;
Ghat_nu = tools.fourierGradQBiPerG(R_b, y, alpha, k, NN, N, m);
Ghat = tools.fourierQBiPerG(R_b, y, alpha, k, NN, N, m);
Ghat_fft = tools.fftQBiPerG(R_b, y, alpha, k, NN, N);
Ghat_nu_fft = tools.fftGradQBiPerG(R_b, y, alpha, k, NN, N);
plot(n,real(Ghat_fft),n,imag(Ghat_fft))
hold on
plot(n,real(Ghat),'--*',n,imag(Ghat),'--*')
legend('real(G_{fft})', 'imag(G_{fft})','real(G_{mpl})', 'imag(G_{mpl})');
hold off
for l = 1:M
    G_fourier(l) = sum(Ghat.*exp(1i.*n.*theta(l)));
    G_nu(l) = sum(Ghat_nu.*exp(1i.*n.*theta(l)));
    G_fourier_fft(l) = sum(Ghat_fft.*exp(1i.*n.*theta(l)));
    G_nu_fft(l) = sum(Ghat_nu_fft.*exp(1i.*n.*theta(l)));
end

x1 = R_b*cos(theta);
x2 = R_b*sin(theta);
d = 1; % Periodicity

G_ewald = ops.GBiPeriodic(k, x1-y(1), x2-y(2), d, d, alpha);
G_ewald_nu = dot(ops.GradGBiPeriodic(k, x1-y(1), x2-y(2), d, d, alpha),[cos(theta);sin(theta)]);
figure
plot(theta,real(G_fourier),'--*',theta,imag(G_fourier),'--*');
hold on
plot(theta,real(G_ewald),theta,imag(G_ewald));
legend('real(G_{fourier})','imag(G_{fourier})','real(G_{ewald})','imag(G_{ewald})');
hold off
figure
plot(theta,real(G_nu),'--*',theta,imag(G_nu),'--*');
hold on
plot(theta,real(G_ewald_nu),theta,imag(G_ewald_nu));
legend('real(G_{\nu,fourier})','imag(G_{\nu,fourier})','real(G_{\nu,ewald})','imag(G_{\nu,ewald})');