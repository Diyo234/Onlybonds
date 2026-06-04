N = min(2000, size(beamformed_data, 3));
data = beamformed_data(:, :, 1:N);

nCut = 4;
nNoiseCut = 0;

data_filt = svdClutterFilter(data, nCut, nNoiseCut);

PD = sum(abs(data_filt).^2, 3);

PD_dB = 10*log10(PD + eps);
PD_dB = PD_dB - max(PD_dB(:));

pixelSize = 3.3879*1e-5;
[H, W] = size(PD_dB);

dx = pixelSize;
dz = pixelSize;

xAxis_mm = ((1:W) - 0.5) * dx * 1e3;
zAxis_mm = ((1:H) - 0.5) * dz * 1e3;

figure;
imagesc(xAxis_mm, zAxis_mm, PD_dB, [-40 0]);
impixelinfo;
axis image off;
colormap hot;
colorbar;

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