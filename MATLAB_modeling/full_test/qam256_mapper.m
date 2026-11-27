function [I, Q] = qam256_mapper(data, T) %#codegen
%QAM256_MAPPER  256-QAM mapper (3GPP TS 38.211, clause 5.1.5).
%
%   [I, Q] = qam256_mapper(data, T)
%
%   data : 8-bit symbols (0..255), scalar or array. MSB = b0, LSB = b7.
%   T    : struct of type prototypes
%            T.I - data type of the I output
%            T.Q - data type of the Q output
%
%   I, Q : same size as data, scaled by 1/sqrt(170) (unit average power).
%
%   Example:
%     T.I = fi([], 1, 16, 14);
%     T.Q = fi([], 1, 16, 14);
%     [I, Q] = qam256_mapper(uint8(0:255), T);
%     x = complex(I, Q);

    % Per-axis levels, address = {sign, m2, m1, m0}
    levels = [ 5  7  3  1  11  9  13  15, ...
              -5 -7 -3 -1 -11 -9 -13 -15] / sqrt(170);

    LUT_I = cast(levels, 'like', T.I);
    LUT_Q = cast(levels, 'like', T.Q);

    d = uint8(data(:));

    % I address = {b0, b2, b4, b6}
    addrI = bitor(bitor(bitshift(uint8(bitget(d, 8)), 3), bitshift(uint8(bitget(d, 6)), 2)),bitor(bitshift(uint8(bitget(d, 4)), 1), uint8(bitget(d, 2))));

    % Q address = {b1, b3, b5, b7}
    addrQ = bitor(bitor(bitshift(uint8(bitget(d, 7)), 3), bitshift(uint8(bitget(d, 5)), 2)),bitor(bitshift(uint8(bitget(d, 3)), 1), uint8(bitget(d, 1))));

    I = reshape(LUT_I(double(addrI) + 1), size(data));
    Q = reshape(LUT_Q(double(addrQ) + 1), size(data));
end
