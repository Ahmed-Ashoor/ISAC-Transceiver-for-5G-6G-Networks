function T = chirp_rom_data_types(dt, wl)
% dt: 'double', 'single' or 'FxPt'
% wl: ROM word length per component (re or im), used by 'FxPt' only

if nargin < 2
    wl = 16;
end

switch dt
    case 'double'
        T.C = double([]);

    case 'single'
        T.C = single([]);

    case 'FxPt'
        % the ROM contents are computed offline, so round to nearest costs
        % nothing in hardware
        fm = fimath('RoundingMethod','Nearest');
        % cos/sin are in [-1,1]: 2 integer bits (incl. sign) so +1 is exact,
        % wl-2 fraction bits. Same format as T.C in radar_rx_data_types
        T.C = fi([], 1, wl, wl-2, fm);
end

end