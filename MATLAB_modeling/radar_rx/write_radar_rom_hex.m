function write_radar_rom_hex(C)
% Dump the chirp reference ROM (visit / stream order, Q2.14) in hex.
%   C : N x 1 fi complex, the T.C-typed reference used by radar_rx
    fid = fopen('chirp_ref_rom.txt','w');
    for n = 1:numel(C)
        fprintf(fid, '%s%s\n', hex(real(C(n))), hex(imag(C(n))));
    end
    fclose(fid);
end
