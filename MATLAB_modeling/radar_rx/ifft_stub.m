function y = ifft_stub(Y, T, out_order) 
% !! PLACEHOLDER for the IFFT your friend is building. Replace with his model.
% It only exists so the chain can be tested end to end.
%   Y         : multiplier output stream, in the FFT's BIT-REVERSED order
%   out_order : 'natural'  -> output stream in natural order
%                             (what a DIT IFFT that takes bit-reversed input gives)
%               'bitrev'   -> output stream in bit-reversed order
%                             (what a DIF IFFT gives, after a bit-reversal
%                              block on its input)
% Behaviour: bit-reverse the input back to natural order, ideal ifft
% (includes the 1/N), quantise to T.yi (Floor + saturate).

N    = numel(Y);
perm = double(bitrevorder(0:N-1)) + 1;      % 1-based bit-reversal permutation

Yn = zeros(N,1);
Yn(perm) = double(Y);                       % arrival position m -> bin bitrev(m)
yn = ifft(Yn);                              % natural order, 1/N scaling

if strcmp(out_order, 'bitrev')
    ys = yn(perm);                          % stream in bit-reversed order
else
    ys = yn;
end
y = cast(ys, 'like', T.yi);
end
