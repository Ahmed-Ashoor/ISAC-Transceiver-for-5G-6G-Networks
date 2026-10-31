function T = radar_data_types(dt, pmode)
% Types table for the radar RX blocks (deck step 4 / step 6).
%   dt    : 'double' | 'single' | 'FxPt'
%   pmode : 'sq'   -> power = I^2 + Q^2           (2 multipliers)
%           'ambm' -> magnitude = max + 3/8*min   (no multipliers)
%
% Name  Meaning                                    FxPt (WL,FL)
%  R     echo spectrum from the RX FFT              14,4   (= FFT y type)
%  S     reference chirp spectrum from the TX FFT   14,4
%  Sc    conjugated reference                       14,4
%  Y     multiplier output -> IFFT input            14,1  (13 integer bits)
%  yi    IFFT output (placeholder until the real IFFT type is known)  14,1
%  pwr   power / magnitude used by the comparator
%  k38   constant 3/8 for alpha-max-beta-min
%  idx   peak index   K  metres-per-index   dist  distance output

if nargin < 2, pmode = 'ambm'; end

switch dt
    case 'double'
        T.R = double([]);  T.S = double([]);  T.Sc = double([]);
        T.Y = double([]);  T.yi = double([]); T.pwr = double([]);
        T.k38 = double([]);
        T.idx = double([]); T.K = double([]); T.dist = double([]);

    case 'single'
        T.R = single([]);  T.S = single([]);  T.Sc = single([]);
        T.Y = single([]);  T.yi = single([]); T.pwr = single([]);
        T.k38 = single([]);
        T.idx = single([]); T.K = single([]); T.dist = single([]);

    case 'FxPt'
        % Floor = truncation, to match the RTL (same choice as the FFT model)
        fm = fimath('RoundingMethod','Floor');
        T.R   = fi([], 1, 14, 4,  fm);
        T.S   = fi([], 1, 14, 4,  fm);
        T.Sc  = fi([], 1, 14, 4,  fm);
        T.Y   = fi([], 1, 14, 1,  fm);   % +-4096 range; 13 integer bits
        T.yi  = fi([], 1, 14, 1,  fm);
        switch pmode
            case 'sq'    % 2*(2^12)^2 = 2^25 -> 26 integer bits, 2 frac bits
                T.pwr = fi([], 0, 28, 2,  fm);
            case 'ambm'  % max + 0.375*min <= 1.375*2^12 -> 13 integer bits, 1 frac
                T.pwr = fi([], 0, 14, 1,  fm);
            otherwise
                error('pmode must be ''sq'' or ''ambm''');
        end
        T.k38  = fi([], 0, 3, 3,  fm);    % holds 0.375 exactly
        T.idx  = fi([], 0, 12, 0, fm);
        T.K    = fi([], 0, 20, 18, fm);   % metres per index, range [0,4)
        T.dist = fi([], 0, 16, 3, fm);    % up to 8191.875 m, 0.125 m step
end
end
