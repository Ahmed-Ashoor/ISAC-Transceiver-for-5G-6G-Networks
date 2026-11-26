function T = qam256_mapper_data_types(dt)
% dt: 'double', 'single' or 'FxPt'


switch dt
    case 'double'
        T.I = double([]);
        T.Q = double([]);

    case 'single'
        T.I = single([]);
        T.Q = single([]);

    case 'FxPt'
        %using a custom fimath to match the rtl truncating of bits instead
        %of default rounding to nearest
        fm = fimath('RoundingMethod','Floor');
        T.I = fi([], 1, 2+14, 14, fm);
        T.Q = fi([], 1, 2+14, 14, fm);
end

end
