function write_chirp_rom_hex(T)
% Dump the chirp ROM into a text file in hex, one line per entry.
% Inputs:
%   T : struct with field C (type info for casting)

    C = chirp_rom(T);

    fid = fopen('chirp_rom.txt','w');

    for k = 1:numel(C)
        re_hex = hex(real(C(k)));
        im_hex = hex(imag(C(k)));
        fprintf(fid, '%s%s\n', re_hex, im_hex);
    end

    fclose(fid);
end
