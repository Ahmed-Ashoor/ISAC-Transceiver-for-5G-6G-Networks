% Deck steps 3, 4, 6 for the conjugate + multiplier block.
% Needs sdf_fft.m / sdf_fft_data_types.m (and ideally sdf_fft_mex) on the path.
clear; close all;
P     = radar_params();
nSeeds = 200;

Tfft = sdf_fft_data_types('FxPt');
T    = radar_data_types('FxPt');
Td   = radar_data_types('double');
Ts   = radar_data_types('single');

sqnr_fx  = zeros(nSeeds,1);
err_dbl  = zeros(nSeeds,1);
err_sgl  = zeros(nSeeds,1);

for seed = 1:nSeeds
    rng(seed);
    d     = randi([20, P.N-P.Lc-1]);
    alpha = 0.1 + 0.4*rand;
    [R, S] = radar_make_scene(P, d, alpha, 0.1, Tfft, seed);

    if seed == 1
        buildInstrumentedMex radar_mult_core -args {R, S, T};
    end

    Y = radar_mult_core_mex(R, S, T);                 % fixed point (instrumented)

    % Reference: ideal complex math on the SAME quantised FFT outputs
    Yref = double(R) .* conj(double(S));

    % types-table check (deck step 4): double and single must agree with Yref
    Yd = radar_mult_core(double(R), double(S), Td);
    Ys = radar_mult_core(single(R), single(S), Ts);
    err_dbl(seed) = max(abs(Yd - Yref));
    err_sgl(seed) = max(abs(double(Ys) - Yref)) / max(abs(Yref));

    noise_power   = mean(abs(double(Y) - Yref).^2);
    sqnr_fx(seed) = 10*log10(mean(abs(Yref).^2) / noise_power);
end

fprintf('double vs ideal : max error      = %.3e (expect ~1e-12)\n', max(err_dbl));
fprintf('single vs ideal : max rel. error = %.3e (expect ~1e-6)\n',  max(err_sgl));
fprintf('FxPt   SQNR     : mean = %.2f dB, min = %.2f dB\n', mean(sqnr_fx), min(sqnr_fx));

figure; plot(1:nSeeds, sqnr_fx, 'LineWidth', 2); grid on;
xlabel('Seed', 'FontSize', 14); ylabel('Multiplier SQNR (dB)', 'FontSize', 14);

% Integer widths come from this report (deck step 6); fractional widths from SQNR.
showInstrumentationResults radar_mult_core_mex -proposeFL -defaultDT numerictype(1,32)
