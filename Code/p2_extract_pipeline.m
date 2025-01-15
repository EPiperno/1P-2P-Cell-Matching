%% 2P Data
%% Load 2P Data and inspect frame size and number
filename = "m3_d241118_s03_2p_00001";
filepath = "D:\Enrico\Alignment\Data\20241118\meas03\";
save_filepath = filepath + "Outputs\";
tif_filepath = filepath + filename + ".tif";
% Get information about the .tif file
info = imfinfo(tif_filepath);

% Number of frames in the TIFF file
numFrames = numel(info);

% Frame size (assumes all frames have the same size)
frameHeight = info(1).Height;
frameWidth = info(1).Width;

% Output results
fprintf('Number of frames: %d\n', numFrames);
fprintf('Frame size: %d x %d pixels\n', frameHeight, frameWidth);

%% Load TIF data_array
% Open the TIFF file
tiffObj = Tiff(tif_filepath, 'r');

% Preallocate a 3D array for the frames
data = zeros(frameHeight, frameWidth, numFrames, 'like', tiffObj.read());

% Read all frames
for k = 1:numFrames
    tiffObj.setDirectory(k);
    data(:, :, k) = tiffObj.read();
end

% Close the TIFF object
tiffObj.close();

%% Split into n_planes
numArrays = 10
planes_arr = 1:numArrays;

separatedArrays = cell(1, numArrays); % Preallocate a cell array to store results

% Loop through each group
for i = 1:numArrays
    separatedArrays{i} = data(:, :, i:10:end); % Take every 10th frame starting at index i
end

%% Frame cut first 10 frames (i.e. 100 frames from entire .tif file)
frameCut = 10+1;

for i = 1:length(separatedArrays)
    separatedArrays{i} = separatedArrays{i}(:, :, frameCut:end); % Keep frames from frameCut onward
end

p2_data_div = separatedArrays;

p2_data_div_filename = save_filepath + filename + "_p2_data_div.mat";
save(p2_data_div_filename, 'p2_data_div')
%% Plot and visualize ith plane
% Extract the 7th z-stack plane
z_stack_plane = separatedArrays{7};

% Check if there are frames in this z-stack plane
if isempty(z_stack_plane)
    error('The 7th z-stack plane is empty.');
end

% Extract the first frame of the 7th z-stack plane
first_frame = z_stack_plane(:, :, 10);

% Plot the first frame
figure;
imagesc(first_frame); % Display the frame as an image
colormap('gray'); % Use grayscale for visualization
colorbar; % Add a colorbar for intensity reference
axis image; % Ensure the aspect ratio is correct
title('First Frame of the 7th Z-Stack Plane');
xlabel('X Pixels');
ylabel('Y Pixels');

%% Save extract output for z-stack
save_filepath = filepath + "Outputs\";

% Loop through the string arrays
for i = 1:numArrays
    save_filename = save_filepath + p2_filename + sprintf('_f_%d_v2', i)
    output_i = save_extract(separatedArrays{i}, save_filename);


    % Plot the masks for the current plane
    S_p2_i = output_i.spatial_weights; % 512x512xcell_count array
    % Create a figure for the plane
    figure;
    hold on;
    
    % Number of cells in the current plane
    cell_count = size(S_p2_i, 3);
    
    % Generate a colormap for the cells
    colormap = lines(cell_count);
    
    % Loop through each cell in the plane
    for cell_idx = 1:cell_count
        % Extract the mask for the current cell
        cell_mask = S_p2_i(:, :, cell_idx);
        
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
end

%% Plot Results
% Create a figure
figure;
hold on;

% Generate a colormap to assign a unique color to each plane
colormap = lines(length(S_p2));

% Iterate through each plane in S_p2
for plane_idx = 1:length(S_p2)
    plane_data = S_p2{plane_idx}; % Extract the 512x512xcell_count array for the current plane
    
    % Sum across the 3rd dimension to combine masks into a 2D representation
    combined_mask = max(plane_data, [], 3); % Logical OR across masks (equivalent to max for binary data)
    
    % Find the non-zero regions of the combined mask
    [rows, cols] = find(combined_mask);
    
    % Scatter plot for this plane's mask regions, using a unique color
    scatter(cols, rows, 10, colormap(plane_idx, :), 'filled'); 
end

% Finalize the plot
title('Masks from Different Planes in S_{p2}');
xlabel('X Coordinate');
ylabel('Y Coordinate');
axis equal;
hold off;










