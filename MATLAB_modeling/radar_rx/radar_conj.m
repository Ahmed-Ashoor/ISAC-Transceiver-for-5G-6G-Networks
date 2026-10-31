function Sc = radar_conj(S, T) %#codegen
% "Conjugate" block: (a + jb) -> (a - jb). Only the imaginary part changes.
% In RTL this is a negation of the Q word (or a ROM that already stores the
% conjugate, or just a sign change in the multiplier adders).
% Note: negating the most negative value (-2^(W-1)) saturates to +max.

N  = numel(S);
Sc = cast(complex(zeros(N,1)), 'like', T.Sc);
for i = 1:N
    re = real(S(i));
    im = imag(S(i));
    Sc(i) = complex(cast(re, 'like', T.Sc), cast(-im, 'like', T.Sc));
end
end
