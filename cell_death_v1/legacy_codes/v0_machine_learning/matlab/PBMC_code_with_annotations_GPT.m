% Clear workspace, command window, and close all figures
clear; clc; close all;

% Add paths to the MATLAB code directory
addpath(genpath('I:\PBMC proj\MATLAB_code'))

% Set target size for cropping
TargetSize = [180 180 90];

% Change directory to the curated data folder
cd('I:\PBMC proj\PBMC_data_curated')

% Get list of groups in the directory
GroupList = dir();
GroupList(1:2) = []; % Remove '.' and '..' entries

% Initialize structure to hold features
Total_Features = struct;
num = 1;

% Loop through each group
for GroupIdx = 1:length(GroupList)
    cd(GroupList(GroupIdx).folder);
    cd(GroupList(GroupIdx).name);
    
    % Get list of patients in the group
    PatientList = dir('HD*');
    for PatientIdx = 1:length(PatientList)
        cd(PatientList(PatientIdx).folder);
        cd(PatientList(PatientIdx).name);
        DataList = dir('*.mat'); % Get list of .mat files
        
        % Loop through each .mat file
        for DataIdx = 1:length(DataList)
            cd(DataList(DataIdx).folder)
            load(DataList(DataIdx).name) % Load the data
            
            % Crop data to target size
            win = centerCropWindow3d(size(data), TargetSize);

            % Segment the cell using multilevel thresholding
            CellThresh = multithresh(data, 2);
            temp1 = zeros(size(data));
            temp1(data > CellThresh(1)) = 1;
            temp2 = imdilate(temp1, strel('sphere', 2));
            temp3 = bwconncomp(temp2, 26);
            NumPixels = cellfun(@numel, temp3.PixelIdxList);
            MaskIdx = find(NumPixels == max(NumPixels(:)));

            % Create cell mask
            CellMask = zeros(size(data));
            CellMask(temp3.PixelIdxList{MaskIdx}) = 1;

            % Fill holes in the cell mask
            for iter = 1:size(CellMask, 3)
                CellMask(:, :, iter) = imfill(CellMask(:, :, iter), 'holes');
            end

            % Compute centroid of the cell
            CellCenter = regionprops3(CellMask, 'Centroid'); 
            shift = CellCenter.Centroid;
            offset_x = shift(1); offset_y = shift(2); offset_z = shift(3);

            % Center the cell
            offset_x = round(size(data, 1) / 2 - offset_x);
            offset_y = round(size(data, 2) / 2 - offset_y);
            offset_z = round(size(data, 3) / 2 - offset_z);
            CenteredCell = circshift(data, [offset_y offset_x offset_z]);
            CenteredMask = circshift(CellMask, [offset_y offset_x offset_z]);

            % Crop centered data
            CroppedMask = imcrop3(CenteredMask, win);
            CroppedData = imcrop3(CenteredCell, win);
            data = CroppedData;

            % Generate Maximum Intensity Projection (MIP) image
            MIP = squeeze(max(data, [], 3));
            MIP = uint16(rescale(MIP, 0, 65280, 'InputMin', 1.34, 'InputMax', 1.41));

            % Save the processed data and MIP image
            save_path = ['I:\PBMC proj\PBMC_data_curated_final' '\' GroupList(GroupIdx).name '\' PatientList(PatientIdx).name];
            if ~exist(save_path)
                mkdir(save_path)
            end
            cd(save_path)
            save(DataList(DataIdx).name, 'data', 'iteration', 'mag', 'na', 'resx', 'resy', 'resz', 'ri');
            imwrite(MIP, [DataList(DataIdx).name(1:end-4) '.png'])
        end
    end
end

% Clear workspace for the next stage
clear; clc; close all;

% Change directory to the final curated data folder
cd('I:\PBMC proj\PBMC_data_curated_final')

% Get list of groups in the directory
GroupList = dir();
GroupList(1:2) = []; % Remove '.' and '..' entries

% Initialize structure to hold features
Total_Features = struct;
num = 1;

% Loop through each group
for GroupIdx = 1:length(GroupList)
    cd(GroupList(GroupIdx).folder);
    cd(GroupList(GroupIdx).name);
    
    % Get list of patients in the group
    PatientList = dir('HD*');
    for PatientIdx = 1:length(PatientList)
        cd(PatientList(PatientIdx).folder);
        cd(PatientList(PatientIdx).name);
        DataList = dir('*.mat'); % Get list of .mat files
        PngList = dir('*.png'); % Get list of .png files
        for DataIdx = 1:length(DataList)
            cd(DataList(DataIdx).folder)
            if ~exist([DataList(DataIdx).name(1:end-4) '.png'])
                continue
            else
                load(DataList(DataIdx).name) % Load the data
                
                % Segment the cell using multilevel thresholding
                CellThresh = multithresh(data, 2);
                temp1 = zeros(size(data));
                temp1(data > CellThresh(1)) = 1;
                temp2 = imdilate(temp1, strel('sphere', 2));
                temp3 = bwconncomp(temp2, 26);
                NumPixels = cellfun(@numel, temp3.PixelIdxList);
                MaskIdx = find(NumPixels == max(NumPixels(:)));

                % Create cell mask
                CellMask = zeros(size(data));
                CellMask(temp3.PixelIdxList{MaskIdx}) = 1;

                % Fill holes in the cell mask
                for iter = 1:size(CellMask, 3)
                    CellMask(:, :, iter) = imfill(CellMask(:, :, iter), 'holes');
                end
                CellMask = imerode(CellMask, strel('sphere', 2));
                MaskConncomp = bwconncomp(CellMask, 26);
                NumPixels = cellfun(@numel, MaskConncomp.PixelIdxList);
                MaskIdx = find(NumPixels == max(NumPixels(:)));

                % Create cell mask
                CellMask = zeros(size(data));
                CellMask(MaskConncomp.PixelIdxList{MaskIdx}) = 1;
                
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
                MeanRI = mean(data(find(CellMask)));
                RIstd = std(data(find(CellMask)));
                Entropy = entropy(double(rescale(data(find(CellMask)))));
                ProteinDensity = (MeanRI - ri) / 0.185;
                ProteinMass = ProteinDensity * CellVolume;
                
                LDThreshold = 1.39;
                LDMask = (data > LDThreshold) .* CellMask;
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
                    mask2 = (data >= (RI_other_thresh(iter)));
                    mask2 = mask2 - (data > (RI_other_thresh(iter + 1)));
                    
                    FeatureList(FeatureIdx) = sum(mask2(:)) * resx^3;
                    FeatureName{FeatureIdx} = ['Partial Volume_' num2str(RI_other_thresh(iter), '%.3f') '-' num2str(RI_other_thresh(iter + 1), '%.3f')];
                    FeatureIdx = FeatureIdx + 1;
                    
                    partial_VolumeRatio = sum(mask2(:)) * resx^3 / CellVolume;
                    FeatureList(FeatureIdx) = partial_VolumeRatio;
                    FeatureName{FeatureIdx} = ['Partial Volume Ratio_' num2str(RI_other_thresh(iter), '%.3f') '-' num2str(RI_other_thresh(iter + 1), '%.3f')];
                    FeatureIdx = FeatureIdx + 1;

                    if max(mask2(:)) == 0
                        mask2 = mk_ellipse_3D(10, 10, 10, size(data, 1), size(data, 2), size(data, 3));
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
                    
                    FeatureList(FeatureIdx) = (mean(data(find(mask2))) - ri) / 0.185;
                    FeatureName{FeatureIdx} = ['Partial Density_' num2str(RI_other_thresh(iter), '%.3f') '-' num2str(RI_other_thresh(iter + 1), '%.3f')];
                    FeatureIdx = FeatureIdx + 1;
                    
                    FeatureList(FeatureIdx) = (mean(data(find(mask2))) - ri) / 0.185 * sum(mask2(:)) * resx^3;
                    FeatureName{FeatureIdx} = ['Partial Mass_' num2str(RI_other_thresh(iter), '%.3f') '-' num2str(RI_other_thresh(iter + 1), '%.3f')];
                    FeatureIdx = FeatureIdx + 1;
                    
                    FeatureList(FeatureIdx) = (mean(data(find(mask2))) - ri) / 0.185 * sum(mask2(:)) * resx^3 / ProteinMass;
                    FeatureName{FeatureIdx} = ['Partial Mass Ratio_' num2str(RI_other_thresh(iter), '%.3f') '-' num2str(RI_other_thresh(iter + 1), '%.3f')];
                    FeatureIdx = FeatureIdx + 1;                    
                end

                % Save features
                save_path = ['I:\PBMC proj\PBMC_features_curated_final' '\' GroupList(GroupIdx).name '\' PatientList(PatientIdx).name];
                if ~exist(save_path)
                    mkdir(save_path)
                end
                cd(save_path)
                save(['Features-' num2str(DataIdx, '%03d') '.mat'], 'FeatureList', 'FeatureName')
                
                disp(['Group : ' GroupList(GroupIdx).name])
                disp(['Patient : ' PatientList(PatientIdx).name])
                disp(['Data : ' num2str(DataIdx) '/' num2str(length(DataList))])
                
            end
        end
    end
end

% Depth-colorization of the data (creating color images)
[Gmag, ~, ~] = imgradient3(data, 'prewitt');
data(data < 1.34) = 1.34;
data_colored = ones([size(data), 3]);
data_scaled = rescale(data, 0, 1);
MyColor = imresize(CustomColormap, [size(data, 3) 3]);
MyColor = rescale(MyColor, 0, 1);
for zz = 1:size(data, 3)
    data_colored(:, :, zz, 1) = data_scaled(:, :, zz) * MyColor(zz, 1);
    data_colored(:, :, zz, 2) = data_scaled(:, :, zz) * MyColor(zz, 2);
    data_colored(:, :, zz, 3) = data_scaled(:, :, zz) * MyColor(zz, 3);
end
data_MCP = squeeze(max(data_colored, [], 3));
figure, imshow(data_MCP)

% Plot a 3D sphere with surface normals
[X, Y, Z] = sphere(10);
[U, V, W] = surfnorm(X, Y, Z);
quiver3(X, Y, Z, U, V, W, 0)

% Clear workspace for the next stage
clear; clc; close all;

% Change directory to the final curated features folder
cd('I:\PBMC proj\PBMC_features_curated_final')

% Get list of groups in the directory
GroupList = dir();
GroupList(1:2) = []; % Remove '.' and '..' entries

% Initialize structure to hold features
Total_Features = struct;
num = 1;

% Loop through each group
for GroupIdx = 1:length(GroupList)
    cd(GroupList(GroupIdx).folder);
    cd(GroupList(GroupIdx).name);
    
    % Get list of patients in the group
    PatientList = dir('HD*');
    for PatientIdx = 1:length(PatientList)
        cd(PatientList(PatientIdx).folder);
        cd(PatientList(PatientIdx).name);
        FileList = dir('*.mat');
        
        % Load and store features for each file
        for FeatureIdx = 1:length(FileList)
            load(FileList(FeatureIdx).name)
            Total_Features(num).Features = FeatureList;
            Total_Features(num).PatientID = PatientList(PatientIdx).name;
            Total_Features(num).Group = GroupList(GroupIdx).name;
            Total_Features(num).FileName = FileList(FeatureIdx).name;
            num = num + 1;
        end
    end
end

% Combine all features into a single array
TotalFeatures = [];
TotalGroup = {};
TotalPatientID = {};
TotalFileName = {};
num = 1;
for idx = 1:length(Total_Features)
    TotalFeatures = [TotalFeatures; Total_Features(idx).Features];
    TotalGroup{num} = Total_Features(idx).Group;
    TotalPatientID{num} = Total_Features(idx).PatientID;
    TotalFileName{num} = [Total_Features(idx).Group '_' Total_Features(idx).PatientID '_' Total_Features(idx).FileName];
    num = num + 1;
end

% Handle missing data
TotalFeatures(isinf(TotalFeatures(:, 20)), 20) = NaN;
NaN_list = isnan(TotalFeatures); % Empty LD mask
Reduction_list = sum(NaN_list, 2);
TotalGroup(find(Reduction_list)) = [];
TotalPatientID(find(Reduction_list)) = [];
TotalFileName(find(Reduction_list)) = [];
TotalFeatures = rmmissing(TotalFeatures);

% Standardize the features
TotalFeatures_zscore = zscore(TotalFeatures);

% Perform dimensionality reduction using UMAP
addpath(genpath('I:\PBMC proj\MATLAB_code'))

umap_model = UMAP(n_components = 2, min_dist = 1);
umap_results = umap_model.fit_transform(double(TotalFeatures_zscore));

% Define color scheme for the scatter plot
ColorMap = cell(length(TotalGroup), 1);
ColorScheme = [37/255 213/255 66/255; 242/255 73/255 73/255; 41/255 49/255 209/255];
xx = unique(TotalGroup);
for ii = 1:length(xx)
    selected = strcmp(xx(ii), TotalGroup);
    ColorMap(selected) = {ColorScheme(ii, :)};
end 
ColorMap = cell2mat(ColorMap);

% Plot UMAP results
figure,
scatter(umap_results(:, 1), umap_results(:, 2), 7, 'filled', 'CData', ColorMap, 'MarkerFaceAlpha', 0.7), axis image;

% Exclude outliers based on UMAP results
ExclusionList = (umap_results(:, 1) > 5);
TotalGroup(ExclusionList) = [];
TotalFeatures_zscore(ExclusionList, :) = [];
TotalFeatures(ExclusionList, :) = [];
TotalPatientID(ExclusionList) = [];
TotalFileName(ExclusionList) = [];

% Perform UMAP again after excluding outliers
umap_results = umap_model.fit_transform(double(TotalFeatures_zscore));

% Plot UMAP results again
close all
ColorMap = cell(length(TotalGroup), 1);
ColorScheme = [37/255 213/255 66/255; 242/255 73/255 73/255; 41/255 49/255 209/255];
xx = unique(TotalGroup);
for ii = 1:length(xx)
    selected = strcmp(xx(ii), TotalGroup);
    ColorMap(selected) = {ColorScheme(ii, :)};
end 
ColorMap = cell2mat(ColorMap);

figure,
scatter(umap_results(:, 1), umap_results(:, 2), 7, 'filled', 'CData', ColorMap, 'MarkerFaceAlpha', 1), axis image off;

% Plot specific features
close all
for iter = 64
    figure(iter),
    scatter(umap_results(~isoutlier(TotalFeatures_zscore(:, iter)), 1), umap_results(~isoutlier(TotalFeatures_zscore(:, iter)), 2), 7, 'filled',...
        'CData', [153 255 51] ./ 255, 'MarkerFaceAlpha', 'flat', 'AlphaData', rescale(TotalFeatures_zscore(~isoutlier(TotalFeatures_zscore(:, iter)), iter), 0, 1)), axis image off
end

% Plot UMAP results with density contours
addpath("I:\PBMC proj\MATLAB_code\kde2d")
[~, density] = kde2d(umap_results, 2^10, min(umap_results) - 2, max(umap_results) + 2);
density = rescale(density, 0, 1);
[x, y] = meshgrid(linspace(min(umap_results(:, 1) - 2), max(umap_results(:, 1) + 2), size(density, 2)), ...
                  linspace(min(umap_results(:, 2) - 2), max(umap_results(:, 2) + 2), size(density, 1)));
sortedDensity = sort(density(:), 'descend');
idx50 = round(0.2 * length(sortedDensity));
level95 = sortedDensity(idx50);

figure;
scatter(umap_results(:, 1), umap_results(:, 2), 5, 'k', 'filled');
hold on;
s = surf(x, y, density, 'LineStyle', 'none');
view([0, -90]);
colormap hot;
axis image off;
hold on;
[C, h] = contour(x, y, density, [level95 level95], 'LineColor', 'b', 'LineWidth', 3);

% Compute and plot correlation of features
Data_Corr = corr(zscore(FeatureSelected));
schemaball(Data_Corr, FeatureName)
correlationCircles(zscore(FeatureSelected)')
circularGraph(rescale(zscore(FeatureSelected)'))

% Plot individual features
for ParamIdx = 18
    DataGroup = [];
    Data = [];

    for idx = 1:length(TotalFeatures)
        DataGroup = [DataGroup; TotalGroup(idx)];
        Data = [Data; TotalFeatures(idx, ParamIdx)];
    end
    Data(Data == 0) = NaN;
    DataName = FeatureName{ParamIdx};
    GroupName = unique(DataGroup);

    x = [];
    y = [];

    for ii = 1:length(GroupName)
        SelectionList = (strcmp(DataGroup, GroupName(ii)));
        OutlierList = isoutlier(Data(SelectionList));

        temp1 = Data(SelectionList);
        temp1 = temp1(~OutlierList);
        temp2 = DataGroup(SelectionList);
        temp2 = temp2(~OutlierList);

        y = [y; temp1];
        x = [x; temp2];
    end

    figure,
    SwarmPlot = swarmchart(categorical(x), y, '.');
    SwarmPlot.XJitterWidth = 0.5;
    title(FeatureName{ParamIdx})
    xticklabels([])
    
    figure,
    BoxPlot = boxchart(categorical(x), y);
    BoxPlot.BoxWidth = 0.2;
    BoxPlot.MarkerStyle = 'none';
    xticklabels([])
end

% Perform MANOVA on features
for iter = 1:size(TotalFeatures_zscore, 2)
    [~, p] = manova1(TotalFeatures_zscore(:, 5), TotalGroup);
    pause(0.1)
end

% Perform t-tests on features
unique_x = unique(x);
data = [];
clc
for j1 = 1:length(unique_x)
    sep_list = strcmp(x, unique_x(j1));
    data{j1} = y(sep_list);
    disp([unique_x{j1} ': ' num2str(mean(y(sep_list), 'omitnan')) '  ' char(177) ' ' num2str(std(y(sep_list), 'omitnan'))])
    if j1 > 1
        [~, p] = ttest2(rmmissing(data{1})', rmmissing(data{j1})');
        disp(['T test: p value: ' num2str(p)])
    end
end

% Plot heatmap of features
PatientSelectionList = unique(TotalPatientID);
GroupSelectionList = unique(TotalGroup);
FeatureSelected = [];
PatientOrderList = {};
Order = zeros(length(Total_Features), 1);
Patient_num = 1;

for GroupIdx = 1:length(GroupSelectionList)
    for PatientIdx = 1:length(PatientSelectionList)
        SelectionList1 = strcmp(TotalPatientID, PatientSelectionList(PatientIdx));
        SelectionList2 = strcmp(TotalGroup, GroupSelectionList(GroupIdx));
        SelectionList = SelectionList1 & SelectionList2;
        FeatureSelected(Patient_num, :) = mean(TotalFeatures(SelectionList, :));
        PatientOrderList(Patient_num) = PatientSelectionList(PatientIdx);
        Patient_num = Patient_num + 1;
    end
end

% Generate clustergram for the features
CG = clustergram(zscore(FeatureSelected), 'Cluster', 'row');
mycmap = get(gcf, 'Colormap');
colormapeditor
load('MyCmap.mat');
CG.Colormap = CustomColormap;
CG.RowLabels = PatientOrderList;

FeaturesOrder = CG.ColumnLabels;
FeatureName_Ordered = {};
for ii = 1:length(FeaturesOrder)
    FeatureName_Ordered{ii} = FeatureName{str2num(FeaturesOrder{ii})};
end
CG.ColumnLabels = FeatureName_Ordered;

% Perform ANOVA on features and plot results
for iter = 1:length(FeatureName)
    Manova(iter) = anovan(TotalFeatures_zscore(:, iter), {TotalGroup'});
end
LogManova = -log10(Manova);
figure, bar(LogManova)

for ii = 1:length(Manova)
    num = str2num(FeaturesOrder{ii});
    LogManova2(ii) = LogManova(num);
end
figure, bar(LogManova2)
