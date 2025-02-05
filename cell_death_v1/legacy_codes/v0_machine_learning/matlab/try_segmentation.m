file_path = "C:\Users\김민욱\Desktop\20230710-Cell death\converted_230512.160429.CD95_TNF_Ctr.003.CD95.A1.T001P20.h5";
ht_data = h5read(file_path, '/ri/000000');
ht_data = rescale(ht_data, 1.33, 1.4);
figure, orthosliceViewer(ht_data);
%%

% Set thresholds
thres = multithresh(ht_data, 4);
mask = ht_data > thres(2);

% Ensure the mask retains the 3D structure
segmented_data = ht_data;
segmented_data(~mask) = 0;

% Visualize the segmented 3D volume
figure, imagesc(squeeze(segmented_data(:,:,30))),axis image off, colormap gray

figure, orthosliceViewer(segmented_data);
%%

% Set target size for cropping
TargetSize = [300 300 59];

% Crop data to target size
win = centerCropWindow3d(size(ht_data), TargetSize);

% Segment the cell using multilevel thresholding
CellThresh = multithresh(ht_data, 3);
temp1 = zeros(size(ht_data));
temp1(ht_data > CellThresh(1)) = 1;
temp2 = imdilate(temp1, strel('sphere', 2));
temp3 = bwconncomp(temp2, 26);
NumPixels = cellfun(@numel, temp3.PixelIdxList);
MaskIdx = find(NumPixels == max(NumPixels(:)));

% Create cell mask
CellMask = zeros(size(ht_data));
% MaskIdx : the index of the largest connected component
CellMask(temp3.PixelIdxList{MaskIdx}) = 1;

% Fill holes in the cell mask
for iter = 1:size(CellMask, 3)
    CellMask(:, :, iter) = imfill(CellMask(:, :, iter), 'holes');
end
%%
% Compute centroid of the cell
CellCenter = regionprops3(CellMask, 'Centroid');
shift = CellCenter.Centroid;
offset_x = shift(1); offset_y = shift(2); offset_z = shift(3);

% Center the cell
offset_x = round(size(ht_data, 1) / 2 - offset_x);
offset_y = round(size(ht_data, 2) / 2 - offset_y);
offset_z = round(size(ht_data, 3) / 2 - offset_z);
CenteredCell = circshift(ht_data, [offset_y offset_x offset_z]);
CenteredMask = circshift(CellMask, [offset_y offset_x offset_z]);

% Crop centered data
CroppedMask = imcrop3(CenteredMask, win);
CroppedData = imcrop3(CenteredCell, win);
data = CroppedData;

%%
% set resx, ri arbitrarily
resx = 1;
ri = 1.33;
% Extract morphological properties
MorphoProps = regionprops3(CellMask, 'all');
CellVolume = MorphoProps.Volume * resx^3;
SurfaceArea = MorphoProps.SurfaceArea * resx^2;
SVRatio = SurfaceArea / CellVolume;
AspectRatio = MorphoProps.PrincipalAxisLength(2) / MorphoProps.PrincipalAxisLength(1);
Solidity = MorphoProps.Solidity;
Sphericity = pi^(1/3) * (6 * CellVolume)^(2/3) / SurfaceArea;
EquivR = MorphoProps.EquivDiameter * resx / 2;
ConvexVolume = MorphoProps.ConvexVolume * resx^3;
Centroid = MorphoProps.Centroid;

CellThreshold = CellThresh(1);
MeanRI = mean(ht_data(find(CellMask)));
RIstd = std(ht_data(find(CellMask)));
Entropy = entropy(double(rescale(ht_data(find(CellMask)))));
ProteinDensity = (MeanRI - ri) / 0.185;
ProteinMass = ProteinDensity * CellVolume;

%%
LDThreshold = 1.39;
LDMask = (ht_data > LDThreshold) .* CellMask;
if isempty(LDMask)
    LDVolume = NaN;
    LDMeanRI = NaN;
    LDDensity = NaN;
    LDMass = NaN;
    LDNumber = NaN;
    MeanLDSize = NaN;
    LDDistMean = NaN;
    LDDistStd = NaN;
else
    LDVolume = sum(LDMask(:)) * resx^3;
    LDMeanRI = mean(data(find(LDMask)));
    LDDensity = (LDMeanRI - ri) / 0.135;
    LDMass = LDDensity * LDVolume;
    LDConncomp = bwconncomp(LDMask, 26);
    NumPixels = cellfun(@numel, LDConncomp.PixelIdxList);
    NumPixels(NumPixels < 10) = [];

    LDNumber = length(NumPixels);

    if LDNumber == 0
        LDVolume = NaN;
        LDMeanRI = NaN;
        LDDensity = NaN;
        LDMass = NaN;
        LDNumber = NaN;
        MeanLDSize = NaN;
        LDDistMean = NaN;
        LDDistStd = NaN;

    else
        MeanLDSize = LDVolume / LDNumber;
        LDProps = regionprops3(LDMask);
        [IndX, IndY, IndZ] = ind2sub(size(LDMask), find(LDMask));
        LDDist = sqrt((IndX - Centroid(1)).^2 + (IndY - Centroid(2)).^2 + (IndZ - Centroid(3)).^2);
        LDDistMean = mean(LDDist) * resx;
        LDDistStd = std(LDDist) * resx;
    end
end

%%
% Initialize feature lists
FeatureList = [];
FeatureName = {};

% Populate feature lists
FeatureList = [CellVolume, SurfaceArea, SVRatio, AspectRatio, Solidity, Sphericity, EquivR, ConvexVolume,...
    CellThreshold, MeanRI, RIstd, Entropy, ProteinDensity, ProteinMass,...
    LDVolume, LDMeanRI, LDDensity, LDMass, LDNumber, MeanLDSize, LDDistMean, LDDistStd];
FeatureName = {'CellVolume', 'SurfaceArea', 'SVRatio', 'AspectRatio', 'Solidity', 'Sphericity', 'EquivR', 'ConvexVolume',...
    'CellThreshold', 'MeanRI', 'RIstd', 'Entropy', 'ProteinDensity', 'ProteinMass',...
    'LDVolume', 'LDMeanRI', 'LDDensity', 'LDMass', 'LDNumber', 'MeanLDSize', 'LDDistMean', 'LDDistStd'};
FeatureIdx = length(FeatureList) + 1;

RI_other_thresh = 1.355:0.005:1.385;

% Iterate through additional thresholds and extract features
for iter = 1:(length(RI_other_thresh) - 1)
    mask2 = (ht_data >= (RI_other_thresh(iter)));
    mask2 = mask2 - (ht_data > (RI_other_thresh(iter + 1)));

    FeatureList(FeatureIdx) = sum(mask2(:)) * resx^3;
    FeatureName{FeatureIdx} = ['Partial Volume_' num2str(RI_other_thresh(iter), '%.3f') '-' num2str(RI_other_thresh(iter + 1), '%.3f')];
    FeatureIdx = FeatureIdx + 1;

    partial_VolumeRatio = sum(mask2(:)) * resx^3 / CellVolume;
    FeatureList(FeatureIdx) = partial_VolumeRatio;
    FeatureName{FeatureIdx} = ['Partial Volume Ratio_' num2str(RI_other_thresh(iter), '%.3f') '-' num2str(RI_other_thresh(iter + 1), '%.3f')];
    FeatureIdx = FeatureIdx + 1;

    if max(mask2(:)) == 0
        mask2 = mk_ellipse_3D(10, 10, 10, size(ht_data, 1), size(ht_data, 2), size(ht_data, 3));
    end

    [IndX, IndY, IndZ] = ind2sub(size(mask2), find(mask2));
    partial_distance = sqrt((IndX - Centroid(1)).^2 + (IndY - Centroid(2)).^2 + (IndZ - Centroid(3)).^2);
    partial_DistMean = mean(partial_distance) * resx;
    partial_DistStd = std(partial_distance) * resx;

    FeatureList(FeatureIdx) = partial_DistMean;
    FeatureName{FeatureIdx} = ['Spatial Distribution Mean_' num2str(RI_other_thresh(iter), '%.3f') '-' num2str(RI_other_thresh(iter + 1), '%.3f')];
    FeatureIdx = FeatureIdx + 1;

    FeatureList(FeatureIdx) = partial_DistStd;
    FeatureName{FeatureIdx} = ['Spatial Distribution Std_' num2str(RI_other_thresh(iter), '%.3f') '-' num2str(RI_other_thresh(iter + 1), '%.3f')];
    FeatureIdx = FeatureIdx + 1;

    FeatureList(FeatureIdx) = (mean(ht_data(find(mask2))) - ri) / 0.185;
    FeatureName{FeatureIdx} = ['Partial Density_' num2str(RI_other_thresh(iter), '%.3f') '-' num2str(RI_other_thresh(iter + 1), '%.3f')];
    FeatureIdx = FeatureIdx + 1;

    FeatureList(FeatureIdx) = (mean(ht_data(find(mask2))) - ri) / 0.185 * sum(mask2(:)) * resx^3;
    FeatureName{FeatureIdx} = ['Partial Mass_' num2str(RI_other_thresh(iter), '%.3f') '-' num2str(RI_other_thresh(iter + 1), '%.3f')];
    FeatureIdx = FeatureIdx + 1;

    FeatureList(FeatureIdx) = (mean(ht_data(find(mask2))) - ri) / 0.185 * sum(mask2(:)) * resx^3 / ProteinMass;
    FeatureName{FeatureIdx} = ['Partial Mass Ratio_' num2str(RI_other_thresh(iter), '%.3f') '-' num2str(RI_other_thresh(iter + 1), '%.3f')];
    FeatureIdx = FeatureIdx + 1;
end

% Save features
save_path = ['C:\Users\김민욱\Desktop\Cell death' '\test' ];
if ~exist(save_path)
    mkdir(save_path)
end
cd(save_path)
save(['Features' '.mat'], 'FeatureList', 'FeatureName')