% Test of the radar RX path (conjugate -> multiplier -> IFFT -> magnitude -> peaks -> distance)
% FFT and IFFT are the teammate's blocks and are NOT used here:
%   - the FFT output (echo spectrum R) is synthesised from the model in gp_basis.pdf
%       Y[k] = C[k] * sum_i a_i exp(-j2*pi*k*df*tau_i) + noise
%   - the IFFT is a placeholder inside radar_rx.m (ideal ifft, same order/types)
close all; clear;

N      = 4096;
nSeeds = 200;

%% ---------------- parameters (gp_basis.pdf values, otherwise PLACEHOLDERS) ----------------
c       = 299792458;
P.fs    = 122.88e6;          % = N*df, 100 MHz channel / 30 kHz row of the NR table (pdf)
P.df    = 30e3;              % subcarrier spacing, mu = 1 (pdf)
P.Ncp   = 288;               % normal CP at 30 kHz = 2.34 us = 288 samples -> R_max ~ 351 m (pdf)
P.K     = c / (2*P.fs);      % metres per IFFT output index (1.22 m)
P.M     = 1024;              % number of sensing bins           (PLACEHOLDER)
P.s     = 2;                 % sensing bins every 2nd bin       (PLACEHOLDER, interleaved)
P.AR    = 64;                % |R| of a unit reflector, FFT output scale          (PLACEHOLDER)
P.sigR  = 4;                 % noise per I and per Q on every sensing bin         (PLACEHOLDER)
P.shift = 8;                 % Y = R*conj(C)*2^-shift, brings Y into the IFFT input range (PLACEHOLDER)
P.pmode = 'ambm';            % 'ambm' | 'max' | 'sq'
P.imin  = 12;                % first index searched (blanks Tx->Rx leakage near 0)
P.imax  = P.Ncp - 1;         % last index searched (echoes beyond the CP break the model)
P.NT    = 8;                 % CFAR training cells per side
P.NG    = 2;                 % CFAR guard cells per side
P.kcfar = 3.25;              % CFAR scale (use 7.5 for 'sq' = power)
P.MAXPK = 16;                % size of the peak list

T  = radar_rx_data_types('FxPt',   P.pmode);
Td = radar_rx_data_types('double', P.pmode);
Ts = radar_rx_data_types('single', P.pmode);

%% ---------------- sensing-bin layout and chirp reference ROM (PLACEHOLDER chirp) ----------------
kf_s   = (-P.M : P.s : P.M-1).';                 % signed bin numbers of the sensing bins
nat    = mod(kf_s, N);                           % natural FFT bin index
kf_all = [0:N/2-1, -N/2:-1].';                   % signed bin number of every bin
sens   = false(N,1);  sens(nat+1) = true;

Cnat   = zeros(N,1);                             % unit magnitude on sensing bins, 0 elsewhere
Cnat(nat+1) = exp(1j*pi*((0:numel(nat)-1).').^2/numel(nat));   % quadratic phase = chirp

perm   = bitrevorder(0:N-1);  perm = perm(:) + 1;  % the FFT outputs bin bitrev(m) at position m
C      = cast(Cnat(perm), 'like', T.C);          % ROM content in visit (stream) order

%% ---------------- Monte Carlo ----------------
nTgtTot = 0;  nHit = 0;  nFalse = 0;  nMismatch = 0;
sqnrY   = zeros(nSeeds,1);

for seed = 1:nSeeds
    rng(seed);

    % 1..3 targets, at least 12 samples apart, delay inside the CP window
    nT = randi(3);
    d  = zeros(nT,1);
    while true
        d = sort(randi([P.imin+8, P.imax-8], nT, 1));
        if nT == 1 || min(diff(d)) >= 12, break; end
    end
    a  = 0.1 + 0.9*rand(nT,1);                   % reflectivity

    % ---- placeholder for the FFT output (random noise + the model above) ----
    Z  = zeros(N,1);
    for t = 1:nT
        Z = Z + a(t)*exp(-1j*2*pi*kf_all*d(t)/N);
    end
    Rn = P.AR*Cnat.*Z + P.sigR*(randn(N,1) + 1j*randn(N,1)).*sens ...
       + (~sens).*P.AR.*(sign(randn(N,1)) + 1j*sign(randn(N,1)));   % comm bins = random data
    R  = cast(Rn(perm), 'like', T.R);

    if seed == 1
        buildInstrumentedMex('radar_rx', '-args', {R, C, T, coder.Constant(P)});
    end

    % fixed point (instrumented MEX; the constant input P is not part of the MEX signature)
    [Y, yi, npk, pk_idx, pk_mag, pk_dist, ovf] = radar_rx_mex(R, C, T); %#ok<ASGLU>

    % double / single versions of the same algorithm (deck step 4 check)
    [~, ~, npk_d, idx_d] = radar_rx(double(R), double(C), Td, P);
    [~, ~, npk_s, idx_s] = radar_rx(single(R), single(C), Ts, P);
    nMismatch = nMismatch + (npk_d ~= npk_s) + ~isequal(idx_d(1:npk_d), idx_s(1:npk_s));

    % SQNR of the multiplier alone (same quantised inputs, ideal complex math)
    Yref      = double(R) .* conj(double(C)) * 2^-P.shift;
    sqnrY(seed) = 10*log10(mean(abs(Yref).^2) / mean(abs(double(Y) - Yref).^2));

    % score
    found    = double(pk_idx(1:npk));
    nTgtTot  = nTgtTot + nT;
    nHit     = nHit + sum(ismember(d, found));
    nFalse   = nFalse + sum(~ismember(found, d));

    % Write inputs and outputs to files (golden vectors for the RTL testbenches)
%     if seed == 1
%         fid = fopen('radar_R_in.txt','w');  fid2 = fopen('radar_Y_out.txt','w');
%         for k = 1:N
%             fprintf(fid,  '%s%s\n', hex(real(R(k))),  hex(imag(R(k))));
%             fprintf(fid2, '%s%s\n', hex(real(Y(k))),  hex(imag(Y(k))));
%         end
%         fclose(fid); fclose(fid2);
%         % detector input/outputs: yi, pk_idx, pk_mag, pk_dist
%     end
end

% write_radar_rom_hex(C);

fprintf('multiplier SQNR           = %.2f dB (mean)\n', mean(sqnrY));
fprintf('targets detected          = %d / %d  (%.1f %%)\n', nHit, nTgtTot, 100*nHit/nTgtTot);
fprintf('extra (false) detections  = %d  (%.3f per frame)\n', nFalse, nFalse/nSeeds);
fprintf('double vs single mismatch = %d frames\n', nMismatch);

figure; plot(1:nSeeds, sqnrY, 'LineWidth', 2); grid on;
xlabel('Seed', 'FontSize', 14); ylabel('Multiplier SQNR (dB)', 'FontSize', 14);

% last frame: range profile and the detections
figure; plot((0:P.imax+20)*P.K, abs(double(yi(1:P.imax+21)))); hold on; grid on;
stem(double(pk_idx(1:npk))*P.K, double(pk_mag(1:npk)), 'r', 'LineWidth', 2);
xlabel('Range (m)', 'FontSize', 14); ylabel('|y|', 'FontSize', 14);
legend('IFFT output magnitude', 'detected peaks');

showInstrumentationResults radar_rx_mex -proposeFL -defaultDT numerictype(1,32)
