function T = chirp_rom_data_types(dt, wl)
% dt: 'double' or 'FxPt'
% wl: ROM word length per component (re or im), used by 'FxPt' only

if nargin < 2
    wl = 16;
end

switch dt
    case 'double'
        T.C = double([]);

    case 'FxPt'
        % the ROM contents are computed offline, so round to nearest costs
        % nothing in hardware
        fm = fimath('RoundingMethod','Nearest');
        % cos/sin are in [-1,1]: 1 sign bit, wl-1 fraction bits
        T.C = fi([], 1, wl, wl-1, fm);
end

end
