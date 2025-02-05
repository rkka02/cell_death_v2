%%
% timelapse
timelapse = zeros(1009,1009,47);
for i = 0:46
    tidx = sprintf('%06d', i);
    temp = h5read(file_path, strcat('/ri/', tidx));
    temp = temp(:,:,30);
    timelapse(:,:,i+1) = temp;
end
figure, orthosliceViewer(timelapse);

%%
% Create gif
% Assuming 'data' is your 1009x1009x47 array
data = timelapse;  % Example data, replace with your actual data

% Normalize data to the range [0, 255] for 8-bit image representation
data_min = min(data(:));
data_max = max(data(:));
normalized_data = uint8(255 * (data - data_min) / (data_max - data_min));

% Define the output file name
outputFileName = '3D_array_to_gif.gif';

% Set the delay time (increase this to slow down the GIF)
delayTime = 0.5;  % 0.5 seconds per frame

% Create the GIF file
for z = 1:size(normalized_data, 3)
    % Extract the current slice
    frame = normalized_data(:, :, z);
    
    % Convert the frame to an indexed image
    [indexedFrame, cmap] = gray2ind(frame, 256);
    
    % Write the frame to the GIF file
    if z == 1
        % For the first frame, create the GIF file
        imwrite(indexedFrame, cmap, outputFileName, 'gif', 'LoopCount', Inf, 'DelayTime', delayTime);
    else
        % For subsequent frames, append to the GIF file
        imwrite(indexedFrame, cmap, outputFileName, 'gif', 'WriteMode', 'append', 'DelayTime', delayTime);
    end
end
