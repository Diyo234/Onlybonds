% load('');

N = min(2000, size(beamformed_data, 3));
data = beamformed_data(:, :, 1:N);

% SVDparameter
nCut = 4;
nNoiseCut = 0;

%SVD clutter filtering
data_filt = svdClutterFilter(data, nCut, nNoiseCut);

%Power Doppler
PD = sum(abs(data_filt).^2, 3);

% convert to dB
PD_dB = 10*log10(PD + eps);
PD_dB = PD_dB - max(PD_dB(:));


figure;
imagesc(imaging_x_axis, imaging_z_axis, PD_dB, [-40 0]);
axis image;
set(gca, 'YDir', 'normal');
colormap hot;
colorbar;
xlabel('x (m)');
ylabel('z (m)');
title(sprintf('SVD Power Doppler (nCut=%d, N=%d)', nCut, N));