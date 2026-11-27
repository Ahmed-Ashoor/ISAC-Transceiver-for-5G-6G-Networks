close all; clc; clear;
wlList = 8; %:18;

% double reference
C_ref = chirp_rom(chirp_rom_data_types('double'));

maxErr  = zeros(numel(wlList),1);
sqnr_dB = zeros(numel(wlList),1);

for i = 1:numel(wlList)
    T = chirp_rom_data_types('FxPt', wlList(i));
    C = double(chirp_rom(T));         % quantized ROM values

N = 4096;  M = 1638;  fs = 122.88e6;
p   = (0:3275).';
idx = mod(p - 1638, N) + 1;

X = complex(zeros(N,1));
X(idx(1:2:end)) = C;              % sensing bins only

x = ifft(X) * N/sqrt(M);          % unit RMS
t = (0:N-1).'/fs * 1e6;           % us
f = diff(unwrap(angle(x))) * fs/(2*pi) / 1e6;

figure;
subplot(3,1,1); plot(t, real(x), t, imag(x)); legend('I','Q'); grid on
subplot(3,1,2); plot(t, abs(x)); ylabel('|x|'); grid on
subplot(3,1,3); plot(t(1:end-1), f); ylabel('Freq (MHz)'); xlabel('Time (\mus)'); grid on

    % ---- Delay test and range profile ----
    c0 = 3e8;  R_true = 1000;               % target range (m), try 0..350
    tau = 2*R_true/c0;

    k = [0:N/2-1, -N/2:-1].';              % signed FFT bin index
    r = ifft(fft(x) .* exp(-1i*2*pi*k*(fs/N)*tau));   % r(t) = x(t - tau)

    Y = fft(r) * sqrt(M)/N;                % undo TX scaling
    Z = Y(idx(1:2:end)) .* conj(C_ref);    % dechirp sensing bins

    Nr = 2048;  df_s = 2*fs/N;             % 60 kHz
    h  = ifft(Z, Nr) * Nr/M;               % ideal peak = 1
    R  = (0:Nr-1).' * c0/(2*Nr*df_s);      % range axis (m)
    [~, ipk] = max(abs(h));

    figure; plot(R, 20*log10(abs(h) + eps)); grid on
    xline(R_true, 'r--'); xlim([0 2500]); ylim([-80 5])
    xlabel('Range (m)'); ylabel('|h| (dB)')
    title(sprintf('wl = %d: peak %.2f m (true %.1f m)', wlList(i), R(ipk), R_true))
    
    C = chirp_rom(T);

    e = double(C) - C_ref;

    maxErr(i)      = max(abs(e));
    signal_power   = mean(abs(C_ref).^2);          % = 1
    noise_power    = mean(abs(e).^2);
    sqnr_dB(i)     = 10*log10(signal_power / noise_power);

    fprintf('wl = %2d : max error = %.2e, sqnr = %.2f dB\n', ...
        wlList(i), maxErr(i), sqnr_dB(i));
end

% figure; plot(wlList, sqnr_dB, '-o', 'LineWidth', 2); grid on;
% xlabel('ROM word length (bits per component)', 'FontSize', 14);
% ylabel('SQNR (dB)', 'FontSize', 14);
% 
% figure; semilogy(wlList, maxErr, '-o', 'LineWidth', 2); grid on;
% xlabel('ROM word length (bits per component)', 'FontSize', 14);
% ylabel('Max error', 'FontSize', 14);

% chosen word length -> dump ROM file
% T = chirp_rom_data_types('FxPt', 14);
% write_chirp_rom_hex(T);
