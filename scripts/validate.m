pkg load splines
pkg load signal

clear;
clf;

[M N L F_s SPL_ref delta beta] = num2cell(load("parameters")){:};
H_s = load("sampled")(:, 2);
R = load("realized")(:, 2:3);
H_r = complex(R(:, 1), R(:, 2));

f_iso = nthargout (2, @iso226, 90);

pp_ref = contour(SPL_ref);
pp = contour(SPL_ref + delta);

x = linspace(0, F_s / 2, L / 2 + 1);
x(1) = 0.1;

H_0 = exp(0.05 * log(10) * (ppval (pp, log(x)) - ppval (pp_ref, log(x))));
printf ("Max. sampled response deviation: %e\n", max(abs(H_0 - transpose(H_s))));

w_0 = fftshift(ifft([H_0 H_0(end - 1 : -1 : 2)]))((L - M) / 2 + 1 : (L + M) / 2);
printf ("Max. imaginary part of kernel: %e\n", max(abs(imag(w_0))));

w_0 = real(w_0) .* kaiser(M, beta)';
H_r0 = fft(w_0, N)(1 : N / 2 + 1);
printf ("Max. realized response deviation: %e\n", max(abs(H_r0 - transpose(H_r))));

x = logspace(0, 5, 1000);

semilogx(linspace(1, F_s / 2, L / 2 + 1), 20 * log10(H_s), "-.",
         f_iso, iso226(SPL_ref + delta) - iso226(SPL_ref), ".",
         linspace(1, F_s / 2, N / 2 + 1), 20 * log10(abs(H_r)))
