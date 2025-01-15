scriptpath =  "D:\Enrico\Alignment";
addpath(scriptpath);
filepath = "D:\Enrico\Alignment\Data\20241118\meas04\";
save_filepath = filepath + "Outputs\";

p1_filename = "m3_d241118_s04_1p";
p2_filename = "m3_d241118_s04_2p_00001";
p2_n_planes = 10;

%% Load EXTRACT Results for 1P
cell_extract_filename = save_filepath + p1_filename + "_cell_extract_output_v2.mat";
load(cell_extract_filename, 'S_p1', 'T_p1');

%% Load EXTRACT Results for 2P
S_p2 = cell(1,10);
T_p2 = cell(1,10);

% Loop through each file
for i = 1:10
    % Generate the filename
    filename = save_filepath + p2_filename + sprintf('_f_%d_v2_extract_output.mat', i);
    
    % Load the 'data' variable from the .mat file
    output = load(filename, 'output');
    S_p2_i=output.output.spatial_weights;
    T_p2_i=output.output.temporal_weights;
    
    % Store the 'data' into the cell array
    S_p2{i} = S_p2_i;
    T_p2{i} = T_p2_i;
end 

%% Plotting S_p1 and S_p2 for visualizing alignment
[cellIDs1P, cellIDs2P] = plotAllCells(S_p1, S_p2);

%% Save EXTRACT Results
extract_arrays_filename = save_filepath + p1_filename + "_extract_arrays_v3.mat";
save(extract_arrays_filename, 'S_p1', 'S_p2', 'T_p1', 'T_p2');

%% Load EXTRACT Results for all data arrays
extract_arrays_filename = save_filepath + p1_filename + "_extract_arrays_v3.mat";
load(extract_arrays_filename, 'S_p1', 'S_p2', 'T_p1', 'T_p2');

%% Matching S_p1 and S_p2 resolutions
target_size = [2048, 2048];
S_p2_matched_resolution = matchResolution(S_p2, target_size);

%[cellIDs1P_1, cellIDs1P_1] = plotAllCells(S_p1, S_p2_matched_resolution);
%% Return AFFINE Transform for S_p1 cells, given manually matched cells
csv_filepath = save_filepath + "matching_meas03.csv";
[S_p1_transformed, cellIDs1P_PostAffine, tform] = alignCellsAffine(S_p1, S_p2_matched_resolution, csv_filepath);

%% Plotting S_p1 and S_p2 for visualizing alignment
[cellIDs1P_1, cellIDs1P_1] = plotAllCells(S_p1_transformed, S_p2_matched_resolution);

%% Pre Processing T_p1
% Retain only the rows corresponding to matchedCellIDs
T_p1_filtered = T_p1(cellIDs1P_PostAffine, :);

% Output:
disp('Filtered T_p1 size:');
disp(size(T_p1_filtered)); % Should show [length(matchedCellIDs), 585]

%% Pre Processing T_p2
% Resample T_p2 so that is has same time bins as T_p1
target_time_bins = size(T_p1, 2); % Number of time bins in T_p1
T_p2_resampled = cell(size(T_p2)); % Initialize resampled T_p2

for z = 1:length(T_p2)
    if ~isempty(T_p2{z})
        % Dynamically measure the time bins of the current T_p2 plane
        [cell_count, original_time_bins_count] = size(T_p2{z});
        
        % Generate the original and target time bin vectors
        original_time_bins = 1:original_time_bins_count; % Time bins for the current plane
        new_time_bins = linspace(1, original_time_bins_count, target_time_bins); % Target time bins

        % Preallocate the resampled matrix
        T_p2_resampled{z} = zeros(cell_count, target_time_bins);

        % Interpolate each row in the matrix
        for i = 1:cell_count
            T_p2_resampled{z}(i, :) = interp1(original_time_bins, T_p2{z}(i, :), new_time_bins, 'linear', 'extrap');
        end
    else
        % Keep empty cells as they are
        T_p2_resampled{z} = [];
    end
end

%% Editing variable names

S_p1_final = S_p1_transformed;
S_p2_final = S_p2_matched_resolution;
T_p1_final = T_p1_filtered;
T_p2_final = T_p2_resampled;
cellIDs = {cellIDs1P, cellIDs2P, cellIDs1P_PostAffine};

%% Saving S_p1_final and S_p2_final
transformation_filename = save_filepath + p1_filename + "_pre_processed.mat";
save(transformation_filename, 'S_p1_final', 'S_p2_final', 'T_p1_final', 'T_p2_final', 'cellIDs', '-v7.3');

%% Load Data
transformation_filename = save_filepath + p1_filename + "_pre_processed.mat";
load(transformation_filename, 'S_p1_final', 'S_p2_final', 'T_p1_final', 'T_p2_final', 'cellIDs');