% load('');
% 
% N = min(2000, size(beamformed_data, 3));
% data = beamformed_data(:, :, 1:N);

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



%% If this doesn't work, try and put it through chatgpt once. If it still 
% doesnt work, save your original image
pixelSize = 3.3879*1e-5;


W = 312;
H = 168;

dx = pixelSize;
dz = pixelSize;

xAxis_mm = ((1:W) - 0.5) * dx * 1e3;
zAxis_mm = ((1:H) - 0.5) * dz * 1e3;
figure;
imagesc(xAxis_mm, zAxis_mm, PD_dB, [-40 0]);
set(gca, 'YDir', 'normal');
colormap hot;
colorbar;
axis image off;

hold on;

barLength_mm = 1;

xMin_mm = min(xAxis_mm);
xMax_mm = max(xAxis_mm);
zMin_mm = min(zAxis_mm);
zMax_mm = max(zAxis_mm);

x0 = xMin_mm + 0.05 * (xMax_mm - xMin_mm);
y0 = zMax_mm - 0.05 * (zMax_mm - zMin_mm);

plot([x0, x0 + barLength_mm], [y0, y0], 'w', 'LineWidth', 3);

hold off;





