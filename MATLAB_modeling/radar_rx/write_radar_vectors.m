function write_radar_vectors(prefix, R, S, Y, y, idx, pk, dist)
% Dump golden vectors (hex) for the RTL self-checking testbenches.
% One line per sample: "<real hex><imag hex>" (same style as the FFT dump).
%   R, S : multiplier inputs     Y : multiplier output (expected)
%   y    : peak-detector input   idx, pk, dist : expected detector outputs
% Run after test_radar_chain, e.g. (all fi objects from FxPt types):
%   write_radar_vectors('radar_vec', R, S, Y, y, idx, pk, dist)

dump(sprintf('%s_R.txt', prefix), R);
dump(sprintf('%s_S.txt', prefix), S);
dump(sprintf('%s_Y_expected.txt', prefix), Y);
dump(sprintf('%s_ifft_out.txt',  prefix), y);

fid = fopen(sprintf('%s_peak_expected.txt', prefix), 'w');
fprintf(fid, 'idx  %d\n', double(idx));
fprintf(fid, 'pk   %s\n', hex(pk));
fprintf(fid, 'dist %s\n', hex(dist));
fclose(fid);
end

function dump(fname, v)
fid = fopen(fname, 'w');
for k = 1:numel(v)
    fprintf(fid, '%s%s\n', hex(real(v(k))), hex(imag(v(k))));
end
fclose(fid);
end
