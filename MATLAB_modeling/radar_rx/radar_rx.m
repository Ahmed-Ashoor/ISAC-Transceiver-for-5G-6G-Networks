function [Y, yi, npk, pk_idx, pk_mag, pk_dist, ovf] = radar_rx(R, C, T, P) %#codegen
% Radar RX path after the FFT (see gp_basis.pdf):
%   conjugate -> multiplier -> [IFFT] -> magnitude -> peak detection -> distance
%
%   R : N x 1 complex  echo spectrum, straight from the FFT (BIT-REVERSED stream order)
%   C : N x 1 complex  chirp reference ROM, SAME stream order as R,
%                      unit magnitude on sensing bins and 0 on the comm bins
%   T : types table (radar_rx_data_types)
%   P : parameter struct, passed as a compile-time constant (see test_radar_rx.m)
%
%   Y       : multiplier output (goes to the IFFT, still bit-reversed order)
%   yi      : IFFT output (natural order)  [placeholder here, the real IFFT is the teammate's]
%   npk     : number of peaks found (objects)
%   pk_idx  : delay index of each peak (0-based, uint16)
%   pk_mag  : magnitude (or power) of each peak
%   pk_dist : distance of each peak in metres
%   ovf     : more than P.MAXPK peaks were found

N  = 4096;
sc = 2^-P.shift;      % Y = R*conj(C) * 2^-shift  (binary point move in the rtl)

R = cast(R, 'like', T.R);
C = cast(C, 'like', T.C);

%% 1) conjugate: (a + jb) -> (a - jb)
Sc = cast(complex(zeros(N,1)), 'like', T.Sc);
for i = 1:N
    Sc(i) = complex(cast(real(C(i)), 'like', T.Sc), cast(-imag(C(i)), 'like', T.Sc));
end

%% 2) complex multiplier, bin k meets bin k:  (a+jb)(c+jd) = (ac-bd) + j(ad+bc)
% full precision products and sums, one truncation at the end
Y = cast(complex(zeros(N,1)), 'like', T.Y);
for i = 1:N
    xr = real(R(i));    xq = imag(R(i));
    cr = real(Sc(i));   ci = imag(Sc(i));

    y_re = xr*cr - xq*ci;
    y_im = xr*ci + xq*cr;

    Y(i) = complex(cast(y_re * sc, 'like', T.Y), cast(y_im * sc, 'like', T.Y));
end

%% 3) IFFT (PLACEHOLDER for the teammate's sdf_dit_ifft: bit-reversed in, natural out)
yi = ifft_placeholder(Y, T);

%% 4) magnitude of every IFFT output sample
mag = cast(zeros(N,1), 'like', T.mag);
k38 = cast(0.375, 'like', T.k38);
for n = 1:N
    re = real(yi(n));
    im = imag(yi(n));
    switch P.pmode
        case 'sq'                       % I^2 + Q^2  (2 multipliers)
            mag(n) = cast(re*re + im*im, 'like', T.mag);
        case 'max'                      % max(|I|,|Q|)  (compare + mux only)
            ma = abs(re);
            mb = abs(im);
            if ma > mb
                mag(n) = cast(ma, 'like', T.mag);
            else
                mag(n) = cast(mb, 'like', T.mag);
            end
        otherwise                       % 'ambm': max + 3/8*min  (RTL: mx + (mn>>2) + (mn>>3))
            ma = abs(re);
            mb = abs(im);
            if ma > mb
                mx = ma;  mn = mb;
            else
                mx = mb;  mn = ma;
            end
            mag(n) = cast(mx + mn*k38, 'like', T.mag);
    end
end

%% 5) peak detection: local maximum (sample vs previous and next) + CA-CFAR threshold
%    only inside the search window [P.imin, P.imax] (CP-limited range, leakage blanked)
kcf     = cast(P.kcfar / (2*P.NT), 'like', T.kcf);
pk_idx  = zeros(P.MAXPK, 1, 'uint16');
pk_mag  = cast(zeros(P.MAXPK, 1), 'like', T.mag);
npk     = uint8(0);
ovf     = false;

for n = (P.imin + 1) : (P.imax + 1)            % MATLAB index n  <->  delay index n-1
    m0 = mag(n);
    if m0 > mag(n-1) && m0 >= mag(n+1)         % '>' on one side only: a plateau gives one peak
        % training cells: NT on each side, skipping NG guard cells
        s = cast(0, 'like', T.sum);
        for j = (P.NG + 1) : (P.NG + P.NT)
            s = cast(s + mag(n-j) + mag(n+j), 'like', T.sum);
        end
        thr = cast(s * kcf, 'like', T.mag);    % mean of training cells * K_cfar
        if m0 > thr
            if npk < P.MAXPK
                npk         = npk + 1;
                pk_idx(npk) = uint16(n - 1);
                pk_mag(npk) = m0;
            else
                ovf = true;
            end
        end
    end
end

%% 6) distance = index * c/(2*fs)   (one multiplication by a constant)
Kc      = cast(P.K, 'like', T.K);
pk_dist = cast(zeros(P.MAXPK, 1), 'like', T.dist);
for q = 1:P.MAXPK
    iq         = cast(pk_idx(q), 'like', T.idx);
    pk_dist(q) = cast(iq * Kc, 'like', T.dist);
end

end

%% ------------------------------------------------------------------------
function yi = ifft_placeholder(Y, T)
% PLACEHOLDER: ideal ifft (includes 1/N) so the chain can be tested.
% Same interface as sdf_dit_ifft: bit-reversed order in, natural order out.
N    = numel(Y);
perm = bit_reverse_perm(N);
Yd   = double(Y);
Yn   = complex(zeros(N,1));
for i = 1:N
    Yn(perm(i)) = Yd(i);      % arrival position i-1 carries bin bitrev(i-1)
end
yi = cast(ifft(Yn), 'like', T.yi);
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
