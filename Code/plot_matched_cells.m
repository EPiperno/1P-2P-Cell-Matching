function plot_matched_cells(matched_cells, S_p1, S_p2, plot_only_matched)
% PLOT_MATCHED_CELLS Visualizes the matched cells between 1P and 2P modalities
% 
% Inputs:
%   matched_cells: 3 x matched_cell_count array
%       Row 1: Cell IDs from 1P
%       Row 2: Cell IDs from 2P
%       Row 3: Plane indices from 2P
%   S_p1: Spatial data for 1P cells (resolution_x x resolution_y x cell_count)
%   S_p2: Cell array of 2P spatial data, where each element is (resolution_x x resolution_y x cell_count)
%   plot_only_matched: Boolean flag to plot only matched cells if true

    % Create a new figure
    figure;

    % Get dimensions of spatial data
    [res_x, res_y, num_cells_1P] = size(S_p1);
    num_planes = numel(S_p2);

    % Ensure square aspect ratio by padding if necessary
    max_dim = max(res_x, res_y);
    padded_1P = zeros(max_dim, max_dim, num_cells_1P);
    padded_1P(1:res_x, 1:res_y, :) = S_p1;
    S_p1 = padded_1P;

    for plane_idx = 1:num_planes
        if ~isempty(S_p2{plane_idx})
            [res_x_2P, res_y_2P, num_cells_2P] = size(S_p2{plane_idx});
            max_dim_2P = max(res_x_2P, res_y_2P);
            padded_2P = zeros(max_dim_2P, max_dim_2P, num_cells_2P);
            padded_2P(1:res_x_2P, 1:res_y_2P, :) = S_p2{plane_idx};
            S_p2{plane_idx} = padded_2P;
        end
    end

    % Generate colors for matched cells using a colormap
    num_matches = size(matched_cells, 2);
    colors = lines(num_matches); % Generate distinct colors

    % Compute unique matched cell counts
    num_unique_1P = numel(unique(matched_cells(1, :))); % Unique cell IDs from 1P
    num_unique_2P = numel(unique(matched_cells(2, :))); % Unique cell IDs from 2P

    % Compute average plane indices (z-coordinates) for each 1P cell
    unique_1P_cells = unique(matched_cells(1, :));
    avg_plane_indices = zeros(1, numel(unique_1P_cells));
    for idx = 1:numel(unique_1P_cells)
        cell_1P = unique_1P_cells(idx);
        matched_indices = find(matched_cells(1, :) == cell_1P);
        avg_plane_indices(idx) = mean(matched_cells(3, matched_indices));
    end

    % Initialize combined masks
    combined_1P = zeros(max_dim, max_dim, 3); % RGB image for 1P
    combined_2P = zeros(max_dim, max_dim, 3); % RGB image for 2P

    if ~plot_only_matched
        % Overlay all unmatched cells in white
        for cell_idx = 1:num_cells_1P
            mask_1P = S_p1(:, :, cell_idx) > 0;
            for channel = 1:3
                combined_1P(:, :, channel) = combined_1P(:, :, channel) + mask_1P;
            end
        end

        for plane_idx = 1:num_planes
            if ~isempty(S_p2{plane_idx})
                [~, ~, num_cells_2P] = size(S_p2{plane_idx});
                for cell_idx = 1:num_cells_2P
                    mask_2P = S_p2{plane_idx}(:, :, cell_idx) > 0;
                    for channel = 1:3
                        combined_2P(:, :, channel) = combined_2P(:, :, channel) + mask_2P;
                    end
                end
            end
        end
    end

    % Overlay matched cells with distinct colors
    for match_idx = 1:num_matches
        cell_1P = matched_cells(1, match_idx); % 1P cell ID
        cell_2P = matched_cells(2, match_idx); % 2P cell ID
        plane_idx = matched_cells(3, match_idx); % 2P plane ID

        % Extract color for this match
        color = colors(match_idx, :);

        % Overlay on 1P combined mask
        mask_1P = S_p1(:, :, cell_1P) > 0;
        for channel = 1:3
            combined_1P(:, :, channel) = combined_1P(:, :, channel) + mask_1P * color(channel);
        end

        % Overlay on 2P combined mask
        if ~isempty(S_p2{plane_idx})
            mask_2P = S_p2{plane_idx}(:, :, cell_2P) > 0;
            for channel = 1:3
                combined_2P(:, :, channel) = combined_2P(:, :, channel) + mask_2P * color(channel);
            end
        end
    end

    % Plot side-by-side overlays with updated titles
    subplot(1, 2, 1);
    imagesc(combined_1P);
    title(sprintf('1P Cells (%d Unique Matched)', num_unique_1P)); % Updated title for 1P
    axis square;
    axis off;

    % Annotate each unique 1P cell with the average plane index
    hold on; % Enable adding annotations
    for idx = 1:numel(unique_1P_cells)
        % Get cell mask
        cell_1P = unique_1P_cells(idx);
        mask = S_p1(:, :, cell_1P) > 0;

        % Find the centroid of the cell
        stats = regionprops(mask, 'Centroid');
        if ~isempty(stats)
            centroid = stats.Centroid; % Centroid of the cell

            % Add text annotation
            text(centroid(1), centroid(2), sprintf('%.1f', avg_plane_indices(idx)), ...
                'Color', 'black', 'FontSize', 8, 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', 'FontWeight', 'bold');
        end
    end
    hold off;

    subplot(1, 2, 2);
    imagesc(combined_2P);
    title(sprintf('2P Cells (%d Unique Matched)', num_unique_2P)); % Updated title for 2P
    axis square;
    axis off;

    sgtitle('Matched Cells Visualization: 1P vs 2P');
end
