N = 4096;
nSeeds = 400;

T = qam256_mapper_data_types('FxPt');


error = zeros(nSeeds,1);
sqnr_dB = zeros(nSeeds,1);

% Open files for writing
% fid_in  = fopen('qam_inputs.txt','w');
% fid_out = fopen('qam_outputs.txt','w');

for seed = 1:nSeeds
    rng(seed);

    data = uint8(randi([0 255], N, 1));

    if seed == 1
    buildInstrumentedMex qam256_mapper -args {data, T};
    end

    [I, Q] = qam256_mapper_mex(data, T);

    %reference from 3GPP TS 38.211 equation (column 1 = b0 = MSB)
    s     = 1 - 2*(dec2bin(data, 8) - '0');
    I_ref = s(:,1).*(8 - s(:,3).*(4 - s(:,5).*(2 - s(:,7)))) / sqrt(170);
    Q_ref = s(:,2).*(8 - s(:,4).*(4 - s(:,6).*(2 - s(:,8)))) / sqrt(170);

    y     = complex(double(I), double(Q));
    y_ref = complex(I_ref, Q_ref);

    error(seed) = mean(abs(y - y_ref)./abs(y_ref))*100;
    signal_power   = mean(abs(y_ref).^2);
    noise_power    = mean(abs(y - y_ref).^2);
    sqnr_dB(seed)  = 10*log10(signal_power / noise_power);

    % Write inputs and outputs to files
%     for k = 1:N
%         % Input: 8-bit symbol in hex
%         in_hex = dec2hex(data(k), 2);
%
%         % Output: I and Q in hex (signed fixed point)
%         out_real_hex = hex(I(k));
%         out_imag_hex = hex(Q(k));
%
%         % Write to files
%         fprintf(fid_in,  '%s\n', in_hex);
%         fprintf(fid_out, '%s %s\n', out_real_hex, out_imag_hex);
%     end
end

% fclose(fid_in);
% fclose(fid_out);

avg_error   = mean(error);
avg_sqnr_dB = mean(sqnr_dB);

fprintf('average error percentage = %.4f %%\n', avg_error);
fprintf('average sqnr = %.2f dB\n', avg_sqnr_dB);

figure; plot(1:nSeeds, error, 'LineWidth', 2); grid on;
xlabel('Seed', 'FontSize', 14); ylabel('Error', 'FontSize', 14);

figure; plot(1:nSeeds, sqnr_dB, 'LineWidth', 2); grid on;
xlabel('Seed', 'FontSize', 14); ylabel('SQNR (dB)', 'FontSize', 14);
showInstrumentationResults qam256_mapper_mex -proposeFL -defaultDT numerictype(1,32)
