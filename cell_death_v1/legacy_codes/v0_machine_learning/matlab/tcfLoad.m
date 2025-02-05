function [riTomo, resolutionXY, resolutionZ, deviceInfo] = tcfLoad(tcfFilePath, datasetNo)

% tcfPath = 'D:\Delete later\A2_1 0000NA70.TCF';

riTomo = single(h5read(tcfFilePath, ['/Data/3D/', num2str(datasetNo, '%06.f')]))/10000; %%https://www.mathworks.com/help/matlab/ref/h5read.html
resolutionX = h5readatt(tcfFilePath,'/Data/3D','ResolutionX'); %um per one pixel
resolutionY = h5readatt(tcfFilePath,'/Data/3D','ResolutionY'); %um per one pixel
resolutionXY = [resolutionX  resolutionY];
resolutionZ = h5readatt(tcfFilePath,'/Data/3D','ResolutionZ'); %um per one pixel

deviceInfo = struct();
deviceInfo.na       = h5readatt(tcfFilePath, '/Info/Device', 'NA');
deviceInfo.mag      = h5readatt(tcfFilePath, '/Info/Device', 'Magnification');
% deviceInfo.IterNo   = h5readatt(tcfFilePath, '/Info/Device', 'Iteration');


% figure(); orthosliceViewer(riTomogram);
% figure();imagesc(max(ri,[],3));axis image; title('MIP');
% figure, orthosliceViewer(ri, 'ScaleFactors', [1 1 3]),colorbar;
end