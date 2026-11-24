function T = radar_rx_data_types(dt, pmode)
% dt   : 'double' | 'single' | 'FxPt'
% pmode: 'ambm' (default) | 'max' | 'sq'    -> only changes the magnitude/power types

if nargin < 2, pmode = 'ambm'; end

switch dt
    case 'double'
        T.R    = double([]);   % echo spectrum (FFT output)
        T.C    = double([]);   % chirp reference ROM
        T.Sc   = double([]);   % conjugated reference
        T.Y    = double([]);   % multiplier output (IFFT input)
        T.yi   = double([]);   % IFFT output
        T.mag  = double([]);   % magnitude / power
        T.sum  = double([]);   % CFAR training-cell sum
        T.kcf  = double([]);   % CFAR scale constant
        T.k38  = double([]);   % constant 3/8
        T.idx  = double([]);   % peak index
        T.K    = double([]);   % metres per index
        T.dist = double([]);   % distance

    case 'single'
        T.R    = single([]);
        T.C    = single([]);
        T.Sc   = single([]);
        T.Y    = single([]);
        T.yi   = single([]);
        T.mag  = single([]);
        T.sum  = single([]);
        T.kcf  = single([]);
        T.k38  = single([]);
        T.idx  = single([]);
        T.K    = single([]);
        T.dist = single([]);

    case 'FxPt'
        % Floor = truncation, to match the rtl (same as the FFT / IFFT models)
        fm = fimath('RoundingMethod','Floor');

        T.R    = fi([], 1, 10+4,  4,  fm);  % = sdf_fft output type (y)
        T.C    = fi([], 1, 2+14,  14, fm);  % reference ROM, same format as the IFFT twiddles
        T.Sc   = fi([], 1, 2+14,  14, fm);
        T.Y    = fi([], 1, 4+12,  12, fm);  % = sdf_dit_ifft input type (x)
        T.yi   = fi([], 1, -4+16, 14, fm);  % = sdf_dit_ifft output type (y), range +-0.125

        switch pmode
            case 'sq'    % I^2+Q^2 <= 2*0.125^2 -> 24 bits, 28 fraction bits
                T.mag = fi([], 0, 24, 28, fm);
                T.sum = fi([], 0, 28, 28, fm);   % 16 cells -> 4 extra bits
            otherwise    % 'ambm' / 'max' : <= 1.375*0.125 -> 12 bits, 14 fraction bits
                T.mag = fi([], 0, 12, 14, fm);
                T.sum = fi([], 0, 16, 14, fm);   % 16 cells -> 4 extra bits
        end
        T.kcf  = fi([], 0, 8,  6,  fm);     % K_cfar/(2*NT), range [0,4)
        T.k38  = fi([], 0, 3,  3,  fm);     % holds 0.375 exactly
        T.idx  = fi([], 0, 12, 0,  fm);
        T.K    = fi([], 0, 20, 18, fm);     % metres per index, range [0,4)
        T.dist = fi([], 0, 16, 3,  fm);     % up to 8191.875 m, 0.125 m step
end

end
