function Y = radar_cmult(R, Sc, T) %#codegen
% "Multiplier" block: element-wise complex multiply, bin k meets bin k.
%   (a + jb) * (c + jd) = (ac - bd) + j(ad + bc)      -> 4 multipliers, 2 adders
% Products and sums stay in full precision; one cast (Floor + saturate) to
% T.Y at the end, like the single truncation stage in the RTL.

N = numel(R);
Y = cast(complex(zeros(N,1)), 'like', T.Y);
for i = 1:N
    a = real(R(i));   b = imag(R(i));
    c = real(Sc(i));  d = imag(Sc(i));

    ac = a * c;
    bd = b * d;
    ad = a * d;
    bc = b * c;

    y_re = ac - bd;
    y_im = ad + bc;

    Y(i) = complex(cast(y_re, 'like', T.Y), cast(y_im, 'like', T.Y));
end
end
