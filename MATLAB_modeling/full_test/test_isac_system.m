%% isac_system_test.m
% Main test script for the OFDM ISAC chain.
%
% Builds instrumented MEX files with fixed-point types, runs multiple seeds
% to analyze bit growth with the Fixed-Point Designer, then tests double,
% single, and FxPt data types end-to-end.
%
% Usage:
%   Run directly (F5). The MEX files are built on the first run and reused.
%   After the run, inspect the Fixed-Point Designer windows to check that the
%   proposed word lengths match the actual bit growth.
%
% Parameter notes (check before running):
%   - P.imin = 12 must be >= P.NG + P.NT (here 9), to keep the CFAR window
%     inside the buffer and blank the TX leakage at d = 0.
%   - P.kcfar = 5 (amplitude) / 25 (power) are standard values.
%   - G = 64 fills the FFT input range (Tfft.x has 4 integer bits).
%   - P.shift = 8 keeps the IFFT output below 0.125 (T.yi range).
%
% See also: test_sdf_fft.m, test_isac_system.m

clear; clc; close all;

%% ========================================================================
% 1. System parameters
%% ========================================================================
N    = 4096;                 % FFT / IFFT length
Ncp  = 288;                  % cyclic prefix (30 kHz, normal CP)
Nocc = 3276;                 % occupied subcarriers (273 RB)
M    = Nocc/2;               % comm bins = sensing bins = 1638
fs   = 122.88e6;             % sample rate
c0   = 299792458;            % speed of light
G    = 64;                   % RX gain in front of the ADC

P.pmode = 'ambm';            % 'ambm' | 'max' | 'sq'
P.shift = 8;                 % 2^-shift after the multiplier
P.MAXPK = 8;
P.NG    = 1;
P.NT    = 8;
P.imin  = 12;
P.imax  = Ncp;
if strcmp(P.pmode, 'sq')
    P.kcfar = 25;            % power domain
else
    P.kcfar = 5;             % amplitude domain
end
P.K = c0 / (2*fs);           % metres per delay sample

SNR_dB = 20;
tgt_d  = [0     40            120            230          ];
tgt_a  = [1.0   0.5*exp(0.7i) 0.25*exp(-2.1i) 0.1*exp(1.9i)];

nSeeds = 5000;                 % seeds for Fixed-Point Designer analysis
doFPD  = true;               % run multi-seed FPD analysis?
doPlot = true;               % produce figures?

%% ========================================================================
% 2. Bin mapping (interleaved: even = sensing, odd = comm)
%% ========================================================================
u    = (0:Nocc-1).';
bin  = mod(u - Nocc/2, N) + 1;   % occupied band centred on DC (1-based)
sbin = bin(1:2:end);
cbin = bin(2:2:end);
perm = bit_reverse_perm(N);

%% ========================================================================
% 3. Types for all data types
%% ========================================================================
Tqam_dbl = qam256_mapper_data_types('double');
Trom_dbl = chirp_rom_data_types('double');
Ttx_dbl  = sdf_dit_ifft_data_types('double');
Tfft_dbl = sdf_fft_data_types('double');
Trx_dbl  = radar_rx_data_types('double', P.pmode);

Tqam_sgl = qam256_mapper_data_types('single');
Trom_sgl = chirp_rom_data_types('single');
Ttx_sgl  = sdf_dit_ifft_data_types('single');
Tfft_sgl = sdf_fft_data_types('single');
Trx_sgl  = radar_rx_data_types('single', P.pmode);

Tqam_fx  = qam256_mapper_data_types('FxPt');
Trom_fx  = chirp_rom_data_types('FxPt');
Ttx_fx   = sdf_dit_ifft_data_types('FxPt');
Tfft_fx  = sdf_fft_data_types('FxPt');
Trx_fx   = radar_rx_data_types('FxPt', P.pmode);

%% ========================================================================
% 4. Build instrumented MEX files (FxPt types)
%% ========================================================================
fprintf('\n=== Building instrumented MEX files ===\n');

data0 = uint8(zeros(M, 1));
x0    = cast(complex(zeros(N,1)), 'like', Ttx_fx.x);
adc0  = cast(complex(zeros(N,1)), 'like', Tfft_fx.x);
R0    = cast(complex(zeros(N,1)), 'like', Trx_fx.R);
C0    = cast(complex(zeros(N,1)), 'like', Trx_fx.C);

buildInstrumentedMex qam256_mapper -args {data0, Tqam_fx};
buildInstrumentedMex chirp_rom    -args {Trom_fx};
buildInstrumentedMex sdf_dit_ifft -args {x0, Ttx_fx};
buildInstrumentedMex sdf_fft      -args {adc0, Tfft_fx};

% radar_rx uses P at compile time (switch on P.pmode).  Pass it directly;
% all three branches of the switch fit in the chosen T.mag type.
buildInstrumentedMex radar_rx     -args {R0, C0, Trx_fx, P};

fprintf('MEX files built.\n');

%% ========================================================================
% 5. Multi-seed FxPt analysis (for Fixed-Point Designer)
%% ========================================================================
if doFPD
    fprintf('\n=== Running %d seeds with FxPt MEX ===\n', nSeeds);

    clearInstrumentationResults qam256_mapper_mex;
    clearInstrumentationResults chirp_rom_mex;
    clearInstrumentationResults sdf_dit_ifft_mex;
    clearInstrumentationResults sdf_fft_mex;
    clearInstrumentationResults radar_rx_mex;

    sqnr_tx_s  = zeros(nSeeds, 1);
    sqnr_fft_s = zeros(nSeeds, 1);
    sqnr_rng_s = zeros(nSeeds, 1);
    nerr_s     = zeros(nSeeds, 1);
    nmiss_s    = zeros(nSeeds, 1);
    nfalse_s   = zeros(nSeeds, 1);

    % Reference constellation for the ideal demapper
    [Iref, Qref] = qam256_mapper(uint8(0:255), qam256_mapper_data_types('double'));
    ref = complex(Iref(:), Qref(:));

    for s = 1:nSeeds
        rng(s);

        % --- TX ---
        data   = uint8(randi([0 255], M, 1));
        [I, Q] = qam256_mapper_mex(data, Tqam_fx);
        C      = chirp_rom_mex(Trom_fx);

        X       = cast(complex(zeros(N,1)), 'like', Ttx_fx.x);
        X(cbin) = cast(complex(I, Q), 'like', Ttx_fx.x);
        X(sbin) = cast(C, 'like', Ttx_fx.x);

        x  = sdf_dit_ifft_mex(X(perm), Ttx_fx);
        tx = [x(N-Ncp+1:N); x];

        % --- comm loopback ---
        r_comm = channel(tx, 0, 1, Inf);
        adc    = cast(G * r_comm(Ncp+1 : Ncp+N), 'like', Tfft_fx.x);
        R      = sdf_fft_mex(adc, Tfft_fx);

        Rnat       = complex(zeros(N,1));
        Rnat(perm) = double(R);
        rxsym      = Rnat(cbin) / G;

        [~, k] = min(abs(rxsym - ref.'), [], 2);
        nerr_s(s) = sum(uint8(k-1) ~= data);

        sqnr_tx_s(s)  = sqnr_dB(ifft(double(X)), double(x));
        sqnr_fft_s(s) = sqnr_dB(fft(double(adc)), double(R));

        % --- radar ---
        r_rad = channel(tx, tgt_d, tgt_a, SNR_dB);
        adc   = cast(G * r_rad(Ncp+1 : Ncp+N), 'like', Tfft_fx.x);
        R     = sdf_fft_mex(adc, Tfft_fx);

        Cgrid       = cast(complex(zeros(N,1)), 'like', Trx_fx.C);
        Cgrid(sbin) = cast(C, 'like', Trx_fx.C);

        [Y, yi, npk, pk_idx, ~, ~, ~] = radar_rx_mex(R, Cgrid(perm), Trx_fx, P);

        Ynat       = complex(zeros(N,1));
        Ynat(perm) = double(Y);
        sqnr_rng_s(s) = sqnr_dB(ifft(Ynat), double(yi));

        found = double(pk_idx(1:double(npk)));
        found = found(:);
        truth = tgt_d(tgt_d >= P.imin & tgt_d <= P.imax);
        nm = 0;
        for q = 1:numel(truth)
            if ~any(abs(found - truth(q)) <= 1), nm = nm + 1; end
        end
        nf = 0;
        for q = 1:numel(found)
            if ~any(abs(truth - found(q)) <= 1), nf = nf + 1; end
        end
        nmiss_s(s)  = nm;
        nfalse_s(s) = nf;

        if mod(s, 10) == 0 || s == nSeeds
            fprintf('  seed %3d/%d: SQNR TX %.1f | FFT %.1f | rng %.1f dB | err %d\n', ...
                s, nSeeds, sqnr_tx_s(s), sqnr_fft_s(s), sqnr_rng_s(s), nerr_s(s));
        end
    end

    fprintf('\n=== FxPt summary over %d seeds ===\n', nSeeds);
    fprintf('  SQNR TX  : mean %.1f dB, min %.1f dB\n', mean(sqnr_tx_s), min(sqnr_tx_s));
    fprintf('  SQNR FFT : mean %.1f dB, min %.1f dB\n', mean(sqnr_fft_s), min(sqnr_fft_s));
    fprintf('  SQNR rng : mean %.1f dB, min %.1f dB\n', mean(sqnr_rng_s), min(sqnr_rng_s));
    fprintf('  comm err : %d total, %d seeds with errors\n', sum(nerr_s), sum(nerr_s>0));
    fprintf('  radar    : %d missed, %d false\n', sum(nmiss_s), sum(nfalse_s));

    % --- Fixed-Point Designer results ---
    fprintf('\n=== Fixed-Point Designer results (inspect UI) ===\n');
    showInstrumentationResults sdf_fft_mex -proposeFL -defaultDT numerictype(1,32);
    showInstrumentationResults sdf_dit_ifft_mex -proposeFL -defaultDT numerictype(1,32);
    showInstrumentationResults radar_rx_mex -proposeFL -defaultDT numerictype(1,32);
end

%% ========================================================================
% 6. End-to-end system test for each data type
%% ========================================================================
fprintf('\n=== End-to-end system test ===\n');

dt_list = {'double', 'single', 'FxPt'};
Tqam_list = {Tqam_dbl, Tqam_sgl, Tqam_fx};
Trom_list = {Trom_dbl, Trom_sgl, Trom_fx};
Ttx_list  = {Ttx_dbl, Ttx_sgl, Ttx_fx};
Tfft_list = {Tfft_dbl, Tfft_sgl, Tfft_fx};
Trx_list  = {Trx_dbl, Trx_sgl, Trx_fx};

% Function handles (MEX for FxPt, original for double/single)
f_qam_list   = {@qam256_mapper,     @qam256_mapper,     @qam256_mapper_mex};
f_chirp_list = {@chirp_rom,         @chirp_rom,         @chirp_rom_mex};
f_ifft_list  = {@sdf_dit_ifft,      @sdf_dit_ifft,      @sdf_dit_ifft_mex};
f_fft_list   = {@sdf_fft,           @sdf_fft,           @sdf_fft_mex};
f_radar_list = {@radar_rx,          @radar_rx,          @radar_rx_mex};

res_all = cell(3, 1);

for i = 1:3
    dt = dt_list{i};

    res = run_system_pass(f_qam_list{i}, f_chirp_list{i}, f_ifft_list{i}, ...
        f_fft_list{i}, f_radar_list{i}, ...
        Tqam_list{i}, Trom_list{i}, Ttx_list{i}, Tfft_list{i}, Trx_list{i}, ...
        P, N, Ncp, Nocc, M, G, sbin, cbin, perm, ...
        tgt_d, tgt_a, SNR_dB, 1, doPlot);

    res_all{i} = res;

    fprintf('\n=== %s ===\n', dt);
    fprintf('  block SQNR [dB] : TX IFFT %.1f | RX FFT %.1f | range IFFT %.1f\n', ...
        res.sqnr(1), res.sqnr(2), res.sqnr(3));
    fprintf('  comm            : %d / %d symbol errors, EVM %.1f dB\n', ...
        res.nerr, M, res.evm_dB);
    fprintf('  radar           : %d peaks, %d missed, %d false, overflow %d\n', ...
        numel(res.idx), res.nmiss, res.nfalse, double(res.ovf));
    if res.pass
        fprintf('  RESULT          : PASS\n');
    else
        fprintf('  RESULT          : FAIL\n');
    end
end

%% ========================================================================
% 7. Summary
%% ========================================================================
fprintf('\n=== Summary ===\n');
for i = 1:3
    if res_all{i}.pass
        s = 'PASS';
    else
        s = 'FAIL';
    end
    fprintf('  %-8s : %s\n', dt_list{i}, s);
end

fprintf('\nScript complete. Inspect the Fixed-Point Designer windows for bit growth.\n');

%% ========================================================================
% Local functions
%% ========================================================================

function res = run_system_pass(f_qam, f_chirp, f_ifft, f_fft, f_radar, ...
    Tqam, Trom, Ttx, Tfft, Trx, P, N, Ncp, Nocc, M, G, ...
    sbin, cbin, perm, tgt_d, tgt_a, SNR_dB, seed, doPlot)
% One pass of the full ISAC system using the given function handles.

    rng(seed);

    % --- TX ---
    data   = uint8(randi([0 255], M, 1));
    [I, Q] = f_qam(data, Tqam);
    C      = f_chirp(Trom);

    X       = cast(complex(zeros(N,1)), 'like', Ttx.x);
    X(cbin) = cast(complex(I, Q), 'like', Ttx.x);
    X(sbin) = cast(C, 'like', Ttx.x);

    x  = f_ifft(X(perm), Ttx);
    tx = [x(N-Ncp+1:N); x];

    % --- comm loopback ---
    r_comm = channel(tx, 0, 1, Inf);
    adc    = cast(G * r_comm(Ncp+1 : Ncp+N), 'like', Tfft.x);
    R      = f_fft(adc, Tfft);

    Rnat       = complex(zeros(N,1));
    Rnat(perm) = double(R);
    rxsym      = Rnat(cbin) / G;

    [Iref, Qref] = qam256_mapper(uint8(0:255), qam256_mapper_data_types('double'));
    ref = complex(Iref(:), Qref(:));
    [~, k] = min(abs(rxsym - ref.'), [], 2);
    nerr   = sum(uint8(k-1) ~= data);
    sent   = ref(double(data) + 1);
    evm_dB = 10*log10(sum(abs(rxsym - sent).^2) / sum(abs(sent).^2));

    sqnr_tx  = sqnr_dB(ifft(double(X)), double(x));
    sqnr_fft = sqnr_dB(fft(double(adc)), double(R));

    % --- radar ---
    r_rad = channel(tx, tgt_d, tgt_a, SNR_dB);
    adc   = cast(G * r_rad(Ncp+1 : Ncp+N), 'like', Tfft.x);
    R     = f_fft(adc, Tfft);

    Cgrid       = cast(complex(zeros(N,1)), 'like', Trx.C);
    Cgrid(sbin) = cast(C, 'like', Trx.C);

    [Y, yi, npk, pk_idx, pk_mag, pk_dist, ovf] = f_radar(R, Cgrid(perm), Trx, P);

    Ynat       = complex(zeros(N,1));
    Ynat(perm) = double(Y);
    sqnr_rng   = sqnr_dB(ifft(Ynat), double(yi));

    found = double(pk_idx(1:double(npk)));
    found = found(:);
    truth = tgt_d(tgt_d >= P.imin & tgt_d <= P.imax);
    nmiss = 0;
    for q = 1:numel(truth)
        if ~any(abs(found - truth(q)) <= 1), nmiss = nmiss + 1; end
    end
    nfalse = 0;
    for q = 1:numel(found)
        if ~any(abs(truth - found(q)) <= 1), nfalse = nfalse + 1; end
    end

    pass = (nerr == 0) && (nmiss == 0) && (nfalse == 0) && ~ovf;

    % --- plots ---
    if doPlot
        figure;
        plot(real(rxsym), imag(rxsym), '.', real(ref), imag(ref), 'r+');
        axis equal; grid on;

        figure;
        n  = (0:P.imax).';
        yd = abs(double(yi));
        plot(n*P.K, 20*log10(yd(n+1) + eps)); hold on;
        plot(found*P.K, 20*log10(yd(found+1) + eps), 'ro');
        grid on; xlabel('range [m]'); ylabel('|yi| [dB]');
    end

    res = struct('pass', pass, 'nerr', nerr, 'evm_dB', evm_dB, ...
        'sqnr', [sqnr_tx, sqnr_fft, sqnr_rng], ...
        'idx', found, 'nmiss', nmiss, 'nfalse', nfalse, 'ovf', ovf, ...
        'dist', double(pk_dist(1:double(npk))));
end

function r = channel(tx, d, a, snr_dB)
% Sum of delayed, scaled copies of the TX stream + complex white noise.
    s = double(tx);
    r = complex(zeros(size(s)));
    for q = 1:numel(d)
        r = r + a(q) * [zeros(d(q),1); s(1:end-d(q))];
    end
    sigma = sqrt(mean(abs(s).^2)) * 10^(-snr_dB/20);
    r = r + sigma/sqrt(2) * complex(randn(size(s)), randn(size(s)));
end

function v = sqnr_dB(ref, out)
    v = 10*log10(mean(abs(ref).^2) /mean(abs(out - ref).^2));
end

function perm = bit_reverse_perm(N)
% perm(m+1) = bitrev(m) + 1   (1-based, for indexing)
    nb   = round(log2(N));
    perm = zeros(N,1);
    for m = 0:N-1
        r = 0;  v = m;
        for k = 1:nb
            r = 2*r + mod(v, 2);
            v = floor(v / 2);
        end
        perm(m+1) = r + 1;
    end
end