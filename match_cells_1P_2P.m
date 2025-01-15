function matched_cells = match_cells_1P_2P(T_p1, S_p1, T_p2, S_p2)
% MATCH_CELLS_1P_2P Matches cells between 1 Photon (1P) and 2 Photon (2P) data modalities
% 
% Inputs:
%   T_p1: Time data of 1P cells (cell_count x time_bins)
%   S_p1: Spatial data of 1P cells (resolution_x x resolution_y x cell_count)
%   T_p2: Cell array of 2P time data, where each element is (cell_count x time_bins)
%   S_p2: Cell array of 2P spatial data, where each element is (resolution_x x resolution_y x cell_count)
% 
% Outputs:
%   matched_cells: 3 x matched_cell_count array
%       Row 1: Cell IDs from 1P
%       Row 2: Cell IDs from 2P
%       Row 3: Plane indices from 2P

    matched_cells = [];

    % Thresholds for correlations (can be tuned)
    time_corr_threshold = 0.3; % Temporal cross-correlation threshold
    spatial_corr_threshold = 0.3; % Spatial correlation threshold
    max_lag = 100; % Maximum lag for cross-correlation

    % Get dimensions of 1P data
    [~, time_bins] = size(T_p1);
    num_cells_p1 = size(S_p1, 3);

    % Number of planes in 2P data
    n_planes = numel(T_p2);

    % Initialize matched cells list
    matches = [];

    total_iterations = num_cells_p1 * n_planes; % Total progress steps
    iteration_count = 0; % Counter for iterations
    progress_checkpoint = 5; % Progress feedback every 5%
    next_progress_update = progress_checkpoint;

    for cell_1P = 1:num_cells_p1
        % Extract temporal and spatial data for current 1P cell
        T1 = T_p1(cell_1P, :);
        S1 = S_p1(:, :, cell_1P);

        for plane_idx = 1:n_planes
            % Skip empty planes
            if isempty(T_p2{plane_idx}) || isempty(S_p2{plane_idx})
                continue;
            end

            % Get number of cells in current 2P plane
            num_cells_p2 = size(T_p2{plane_idx}, 1);

            for cell_2P = 1:num_cells_p2
                % Extract temporal and spatial data for current 2P cell
                T2 = T_p2{plane_idx}(cell_2P, :);
                S2 = S_p2{plane_idx}(:, :, cell_2P);

                % Compute time correlation with lag adjustment
                [xc, lags] = xcorr(T1, T2, max_lag, 'coeff');
                max_time_corr = max(xc);

                % Compute spatial correlation
                S1_flat = S1(:) / norm(S1(:));
                S2_flat = S2(:) / norm(S2(:));
                spatial_corr = dot(S1_flat, S2_flat);

                % Check thresholds
                if max_time_corr > time_corr_threshold && spatial_corr > spatial_corr_threshold
                    matches = [matches; cell_1P, cell_2P, plane_idx];
                end
            end

            % Update progress
            iteration_count = iteration_count + 1;
            progress_percentage = (iteration_count / total_iterations) * 100;
            if progress_percentage >= next_progress_update
                fprintf('Progress: %.2f%%\n', progress_percentage);
                next_progress_update = next_progress_update + progress_checkpoint;
            end
        end
    end

    % Format output
    matched_cells = matches'; % 3 x matched_cell_count
end