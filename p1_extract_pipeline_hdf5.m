%% 1P Data
%% Convert DCIMG to HDF5
addpath('D:\Enrico\Alignment\dcimg');
filename = "m3_d241118_s04_1p";
%filename = "TF19_d241114_s06_1p"
%filepath = "D:\Enrico\Alignment\Data\20241118\meas04\";
filepath = "E:\";

dcimg_filepath = filepath + filename + ".dcimg";
hdf5_filepath = filepath + filename + ".h5";
%convert_dcimg(dcimg_filepath, hdf5_filepath)
%[movie,totalframes,summary] = loadDCIMG(dcimg_filepath);
%movie=LoadDCIMG_2(dcimg_filepath);

h5disp(hdf5_filepath)
%% Inspect HDF5 frame size, number, and image
info = h5info(hdf5_filepath);
%disp(info);

datasetName = '/mov'; % Replace with the correct dataset path
dataInfo = h5info(hdf5_filepath, datasetName);
frameHeight = dataInfo.Dataspace.Size(1);
frameWidth = dataInfo.Dataspace.Size(2);
numFrames = dataInfo.Dataspace.Size(3);
fprintf('Frame size: %dx%d\n', frameHeight, frameWidth);
fprintf('Number of frames: %d\n', numFrames);

frameNumber = 67; % Specify which frame to read
frameData = h5read(hdf5_filepath, datasetName, [1, 1, frameNumber], [frameHeight, frameWidth, 1]);
imshow(frameData, []);

%% Load 1P Data and cut frames
frameCut = 70; % Specify at which frame to cut to eliminate spinning mirror artifact

data = h5read(hdf5_filepath, datasetName);
data = data(:, :, frameCut:end);

%% 
p1_data = data;
p1_data_filename = save_filepath + p1_filename + "_p1_data.mat";
save(p1_data_filename, 'p1_data')
%% Plot integration of all planes in 1P
integrated1P = sum(data, 3);

% Plot the integrated 1P data
figure;
imagesc(integrated1P); % Display the integrated image
colormap('gray'); % Use grayscale for better visualization
colorbar; % Add a colorbar for reference
axis image; % Keep aspect ratio
title('Integrated 1P Data Across Frames');
xlabel('X Pixels');
ylabel('Y Pixels');

%% Save Extract Outputs
save_filename = filename;
save_filepath = "";
output = save_extract(data, save_filepath, save_filename);

%% Return and Save T and S arrays
S_p1 = output.spatial_weights; 
T_p1 = output.temporal_weights;
info = output.info;
config = output.config;
filepath = "D:\Enrico\Alignment\Data\20241118\meas04\";
save_filepath = filepath + "Outputs\";


cell_extract_filename = save_filepath + filename + "_cell_extract_output_v2.mat";
save(cell_extract_filename, 'S_p1', 'T_p1');
%% Plot Outputs
% Plot the masks for the current plane
S_p1 = output.spatial_weights; % 512x512xcell_count array
% Create a figure for the plane
figure;
hold on;

% Number of cells in the current plane
cell_count = size(S_p1, 3);

% Generate a colormap for the cells
colormap = lines(cell_count);

% Loop through each cell in the plane
for cell_idx = 1:cell_count
    % Extract the mask for the current cell
    cell_mask = S_p1(:, :, cell_idx);
    
    % Find the non-zero regions of the mask
    [rows, cols] = find(cell_mask);
    
    % Plot the mask for the current cell
    scatter(cols, rows, 10, colormap(cell_idx, :), 'filled');
end

% Finalize the plot
title('Masks of Cells in a Single Plane');
xlabel('X Coordinate');
ylabel('Y Coordinate');
axis equal;
hold off;

% Force MATLAB to render the figure before proceeding
drawnow;



