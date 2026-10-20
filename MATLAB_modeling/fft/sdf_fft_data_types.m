function T = sdf_fft_data_types(dt)
% dt: 'double' or 'single' (extend with fixed-point cases as needed)


switch dt
    case 'double'
        T.x     = double([]);
        T.WN    = double([]);
        T.stg1  = double([]);
        T.stg2  = double([]);
        T.stg3  = double([]);
        T.stg4  = double([]);
        T.stg5  = double([]);
        T.stg6  = double([]);
        T.stg7  = double([]);
        T.stg8  = double([]);
        T.stg9  = double([]);
        T.stg10 = double([]);
        T.stg11 = double([]);
        T.stg12 = double([]);
        T.y     = double([]);

    case 'single'
        T.x     = single([]);
        T.WN    = single([]);
        T.stg1  = single([]);
        T.stg2  = single([]);
        T.stg3  = single([]);
        T.stg4  = single([]);
        T.stg5  = single([]);
        T.stg6  = single([]);
        T.stg7  = single([]);
        T.stg8  = single([]);
        T.stg9  = single([]);
        T.stg10 = single([]);
        T.stg11 = single([]);
        T.stg12 = single([]);
        T.y     = single([]);
        
    case 'FxPt'
        %using a custom fimath to match the rtl truncating of bits instead
        %of default rounding to nearest
        fm = fimath('RoundingMethod','Floor');
        T.x     = fi([], 1, 4+10, 10, fm);
        T.WN    = fi([], 1, 2+12, 12, fm);
        T.stg1  = fi([], 1, 4+10, 10, fm);
        T.stg2  = fi([], 1, 5+9, 9, fm);
        T.stg3  = fi([], 1, 5+9, 9, fm);
        T.stg4  = fi([], 1, 6+8, 8, fm);
        T.stg5  = fi([], 1, 6+8, 8, fm);
        T.stg6  = fi([], 1, 7+7, 7, fm);
        T.stg7  = fi([], 1, 7+7, 7, fm);
        T.stg8  = fi([], 1, 8+6, 6, fm);
        T.stg9  = fi([], 1, 8+6, 6, fm);
        T.stg10 = fi([], 1, 9+5, 5, fm);
        T.stg11 = fi([], 1, 9+5, 5, fm);
        T.stg12 = fi([], 1, 10+4, 4, fm);
        T.y     = fi([], 1, 10+4, 4, fm);
end

end
