function P = radar_params()
% Single source of truth for the radar RX model parameters (deck step 5:
% "Finalize Design Parameters").
% !! fs, B, Lc, A are ASSUMED values for testing - replace them with the
% !! real project numbers.

P.N      = 4096;            % FFT/IFFT size (matches the FFT model)
P.nbits  = 12;              % log2(N)
P.c      = 299792458;       % speed of light [m/s]
P.fs     = 100e6;           % sample rate [Hz]            (ASSUMED)
P.B      = 40e6;            % chirp bandwidth [Hz]        (ASSUMED)
P.Lc     = 1024;            % chirp length in samples     (ASSUMED)
P.A      = 1.0;             % chirp amplitude (FFT input is Q4.10, +-8)

P.K      = P.c / (2*P.fs);  % metres per output index  (R = idx * K)
P.maxR   = P.N * P.K;       % range at which the index wraps around
P.dRes   = P.c / (2*P.B);   % range resolution (set by bandwidth, not N)
end
