% End-to-end radar RX test:
%   chirp+echo -> FFT (teammate) -> conj -> multiplier -> IFFT (stub) ->
%   power -> running max -> index (bit-reversal applied to the index only) -> distance
% Compares the two power methods ('sq' vs 'ambm') and the two IFFT output
% orders. Replace ifft_stub with the real IFFT model when it is ready.
clear; close all;
P      = radar_params();
nSeeds = 100;
sigma  = 0.5;            % receiver noise per I/Q. Raise it to stress the detector.
 
Tfft   = sdf_fft_data_types('FxPt');
Tm     = radar_data_types('FxPt', 'ambm');   % multiplier / IFFT-stub types
pmodes = {'sq', 'ambm'};
orders = {'natural', 'bitrev'};
 
% ---- build the instrumented MEX files once (deck step 3) ----
[R0, S0] = radar_make_scene(P, 500, 0.3, sigma, Tfft, 1);
buildInstrumentedMex radar_mult_core -args {R0, S0, Tm};
 
y0 = ifft_stub(radar_mult_core_mex(R0, S0, Tm), Tm, 'natural');
for ip = 1:numel(pmodes)
    Tp = radar_data_types('FxPt', pmodes{ip});
    buildInstrumentedMex('radar_peak_core', '-args', ...
        {y0, Tp, coder.Constant(pmodes{ip}), true}, '-o', ...
        sprintf('radar_peak_%s_mex', pmodes{ip}));
end
 
% ---- Monte Carlo ----
for ip = 1:numel(pmodes)
    pmode = pmodes{ip};
    Tp    = radar_data_types('FxPt', pmode);
    peakf = sprintf('radar_peak_%s_mex', pmode);
 
    for io = 1:numel(orders)
        order  = orders{io};
        bitrev = strcmp(order, 'bitrev');
        hit = 0;  derr = zeros(nSeeds,1);
 
        for seed = 1:nSeeds
            rng(seed);
            d     = randi([20, P.N-P.Lc-1]);
            alpha = 0.1 + 0.4*rand;
            [R, S, truth] = radar_make_scene(P, d, alpha, sigma, Tfft, seed);
 
            Y = radar_mult_core_mex(R, S, Tm);
            y = ifft_stub(Y, Tm, order);
 
            % pmode was built in as coder.Constant, so the MEX does not take it
            idx  = feval(peakf, y, Tp, bitrev);
            dist = radar_distance(idx, Tp, P.K);
 
            hit        = hit + (double(idx) == truth.d);
            derr(seed) = abs(double(dist) - truth.dist);
        end
        fprintf('%-5s | %-7s | hit rate = %5.1f %% | max distance error = %.3f m\n', ...
            pmode, order, 100*hit/nSeeds, max(derr));
    end
end
 
% ---- instrumentation reports: integer widths of the registers (deck step 6) ----
showInstrumentationResults radar_mult_core_mex -proposeFL -defaultDT numerictype(1,32)
showInstrumentationResults radar_peak_sq_mex   -proposeFL -defaultDT numerictype(0,32)
showInstrumentationResults radar_peak_ambm_mex -proposeFL -defaultDT numerictype(0,32)
 
% ---- one example plot (natural order for display) ----
[R, S, truth] = radar_make_scene(P, 500, 0.3, 0.1, Tfft, 1);
y = ifft_stub(radar_mult_core_mex(R, S, Tm), Tm, 'natural');
figure; plot((0:P.N-1)*P.K, abs(double(y))); grid on;
xlabel('Range (m)'); ylabel('|y|'); title(sprintf('True target at %.1f m', truth.dist));
 