function dist = radar_distance(idx, T, K) %#codegen
% Core algorithm 3: delay index -> distance.
%   R = idx * c / (2*fs)     (K = c/(2*fs) metres per index)
% One multiplication by a constant (shift-and-add in RTL; synthesis does it).
%   idx  : peak index (uint16, 0..N-1)
%   K    : metres per index from radar_params (double)
%   dist : metres, cast to T.dist

Kc   = cast(K,   'like', T.K);
i    = cast(idx, 'like', T.idx);
dist = cast(i * Kc, 'like', T.dist);
end
