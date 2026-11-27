function C = chirp_rom(T) %#codegen
% Contents of the chirp ROM: one unit-magnitude complex value per sensing
% bin, with quadratic phase  C[m] = exp(j*pi*m^2/M),  m = 0..M-1

M = 1638;                          % sensing bins (proposed, not final)

m  = (0:M-1).';
% phase reduced modulo 2*pi with integers first, so no large-argument error
ph = pi * mod(m.^2, 2*M) / M;
C_double = exp(1i * ph);           % double

C = cast(C_double, 'like', T.C);   % casting

end
