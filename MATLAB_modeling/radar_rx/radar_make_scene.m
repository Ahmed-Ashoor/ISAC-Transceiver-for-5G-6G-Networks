function [R, S, truth] = radar_make_scene(P, d, alpha, sigma, Tfft, seed)
% Test stimulus: one chirp pulse, one target, receiver noise, then the
% teammate's fixed-point FFT. Returns the two spectra as STREAMS in the
% FFT's output (bit-reversed) order.
%   d     : target delay in samples (0 <= d <= N-Lc-1, no wrap-around)
%   alpha : echo amplitude relative to the transmitted chirp
%   sigma : receiver noise std per I and per Q
%   Tfft  : sdf_fft_data_types('FxPt')  -> bit-accurate FFT
%           []                          -> ideal double FFT, same output order
%   R     : echo spectrum        S : reference (transmitted chirp) spectrum

rng(seed);
N = P.N;  Lc = P.Lc;

t     = ((0:Lc-1).' - Lc/2) / P.fs;
kk    = P.B / (Lc / P.fs);                       % chirp rate [Hz/s]
chirp = P.A * exp(1j*pi*kk*t.^2);                % complex baseband LFM

tx = [chirp; zeros(N-Lc, 1)];
rx = alpha*circshift(tx, d) + sigma*(randn(N,1) + 1j*randn(N,1));

S = fft_stream(tx, Tfft);
R = fft_stream(rx, Tfft);

truth.d    = d;
truth.dist = d * P.K;
end

function Xs = fft_stream(x, Tfft)
if isempty(Tfft)
    perm = double(bitrevorder(0:numel(x)-1)) + 1;
    X    = fft(x);
    Xs   = X(perm);                              % bit-reversed stream
else
    xi = cast(x, 'like', Tfft.x);
    if exist('sdf_fft_mex', 'file') == 3
        Xs = sdf_fft_mex(xi, Tfft);
    else
        Xs = sdf_fft(xi, Tfft);                  % slow, but works
    end
end
end
