function [idx, pk, idx_raw] = radar_peak_core(y, T, pmode, bitrev) %#codegen
% Core algorithm 2 (deck step 1): power + running-maximum peak search.
% Samples are processed one per "clock", exactly like the RTL would:
%   power -> compare with the best so far -> keep (power, arrival counter).
%
%   y     : IFFT output stream (N x 1 complex), in ARRIVAL order
%   pmode : 'sq'   power = I^2 + Q^2            (two multipliers)
%           'ambm' magnitude = max(|I|,|Q|) + 3/8*min(|I|,|Q|)   (no multiplier)
%           Both give the same argmax in practice; ambm is within
%           -2.8% / +6.8% of the true magnitude.
%   bitrev: false -> arrival position m IS the delay index (natural-order output)
%           true  -> arrival position m carries delay index bitrev(m);
%                    the index is recovered by reversing the bits of the
%                    counter (pure wiring in RTL, no memory, no reorder)
%   idx     : delay index of the strongest sample (uint16)
%   pk      : its power/magnitude
%   idx_raw : arrival position of the strongest sample (uint16)

N  = numel(y);
nb = round(log2(N));

k38  = cast(0.375, 'like', T.k38);
best = cast(0, 'like', T.pwr);
idx_raw = uint16(0);

for m = 0:N-1
    v  = y(m+1);
    re = real(v);
    im = imag(v);

    switch pmode
        case 'sq'
            p = cast(re*re + im*im, 'like', T.pwr);
        otherwise   % 'ambm'
            a = abs(re);
            b = abs(im);
            if a > b
                mx = a;  mn = b;
            else
                mx = b;  mn = a;
            end
            % RTL: mx + (mn>>2) + (mn>>3)
            p = cast(mx + mn*k38, 'like', T.pwr);
    end

    if p > best          % strict '>' : first maximum wins
        best    = p;
        idx_raw = uint16(m);
    end
end

pk = best;
if bitrev
    idx = bit_reverse(idx_raw, nb);
else
    idx = idx_raw;
end
end

function r = bit_reverse(m, nb)
% Reverse the lowest nb bits of m (wiring only in hardware).
r  = uint32(0);
mm = uint32(m);
for k = 0:nb-1
    r = bitor(bitshift(r, 1), bitand(bitshift(mm, -k), uint32(1)));
end
r = uint16(r);
end
