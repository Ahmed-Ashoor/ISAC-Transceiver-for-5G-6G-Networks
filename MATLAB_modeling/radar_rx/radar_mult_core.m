function Y = radar_mult_core(R, S, T) %#codegen
% Core algorithm 1 (deck step 1): conjugate -> complex multiplier.
%   R : echo spectrum   (N x 1 complex, in the order the FFT outputs it)
%   S : reference chirp spectrum, SAME order as R (the TX FFT output)
%   Y : R .* conj(S), same order, cast to T.Y (IFFT input format)
% Both inputs stay in the FFT's bit-reversed order: the multiplier is
% element-wise, so it needs no reordering as long as R and S match.

Sc = radar_conj(S, T);
Y  = radar_cmult(R, Sc, T);
end
