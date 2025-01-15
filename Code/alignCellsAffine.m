function [S_p1_transformed, matchedCellIDs, tform] = alignCellsAffine(S_p1, S_p2, csv_filepath)
    % alignCellsAffine: Computes affine transformation for spatial alignment
    % of S_p1 to S_p2 based on manual matching cell data from a CSV file.
    %
    % Inputs:
    %   S_p1: 3D matrix of 1P modality data (2048x2048xcell_count)
    %   S_p2: 1x10 cell array of 2048x2048xcell_count matrices for 2P data
    %   csv_filepath: Path to CSV file containing manual matches
    %
    % Outputs:
    %   S_p1_transformed: Transformed 1P data with valid cells only
    %   matchedCellIDs: Array of cell IDs that match the "3rd" plane

    % Load the CSV file
    raw = readcell(csv_filepath);

    % Extract 1P and 2P cell IDs
    cell_ids_p1 = cell2mat(raw(1, 2:end)); % 1P cell IDs (numeric)
    cell_matches_p2 = cell2mat(raw(2, 2:end)); % 2P matches (numeric)

    % Display extracted cell IDs for sanity check
    disp('Extracted 1P Cell IDs:');
    disp(cell_ids_p1);
    disp('Extracted 2P Matches (First Row):');
    disp(cell_matches_p2);

    % Extract centroids of S_p1 and S_p2
    centroids_p1 = extractCentroids(S_p1);

    % Map centroids of S_p2
    centroids_p2 = [];
    currentID = 1;
    centroid_map_2p = containers.Map('KeyType', 'double', 'ValueType', 'any');

    for z = 1:length(S_p2)
        if ~isempty(S_p2{z})
            zCentroids = extractCentroids(S_p2{z});
            for i = 1:size(zCentroids, 1)
                centroid_map_2p(currentID) = zCentroids(i, :);
                currentID = currentID + 1;
            end
        end
    end

    % Initialize matched coordinates
    matched_p1_coords = [];
    matched_p2_coords = [];

    % Iterate through the matches
    for col = 1:length(cell_ids_p1)
        % Get 1P cell ID and corresponding centroid
        p1_cell_id = cell_ids_p1(col);
        if isnan(p1_cell_id) || p1_cell_id < 1 || p1_cell_id > size(S_p1, 3)
            continue; % Skip invalid IDs
        end
        centroid_p1 = centroids_p1(p1_cell_id, :);

        % Get the 2P match and corresponding centroid
        p2_id = cell_matches_p2(col);
        if isnan(p2_id) || ~isKey(centroid_map_2p, p2_id)
            continue; % Skip if 2P ID is invalid or not found
        end
        centroid_p2 = centroid_map_2p(p2_id);

        % Append to matched coordinates
        matched_p1_coords = [matched_p1_coords; centroid_p1]; %#ok<AGROW>
        matched_p2_coords = [matched_p2_coords; centroid_p2]; %#ok<AGROW>
    end

    % Display matched coordinates for sanity check
    disp('Matched P1 Coordinates:');
    disp(matched_p1_coords);
    disp('Matched P2 Coordinates:');
    disp(matched_p2_coords);

    % Validate matched coordinates
    if size(matched_p1_coords, 2) ~= 2 || size(matched_p2_coords, 2) ~= 2 || isempty(matched_p1_coords)
        error('Matched coordinates must be Nx2 arrays. Check the input data.');
    end

    % Compute affine transformation
    tform = fitgeotrans(matched_p1_coords, matched_p2_coords, 'affine');

    % Transform S_p1 to match S_p2
    valid_cells = false(1, size(S_p1, 3)); % Track valid cells
    S_p1_transformed = [];
    matchedCellIDs = []; % Store cell IDs that match the "3rd" plane

    for i = 1:size(S_p1, 3)
        mask = S_p1(:, :, i) > 0;
        [x, y] = find(mask);
        if isempty(x), continue; end

        transformed = transformPointsForward(tform, [x, y]);

        % Remove points outside the FOV
        valid_idx = transformed(:, 1) > 0 & transformed(:, 1) <= 2048 & ...
                    transformed(:, 2) > 0 & transformed(:, 2) <= 2048;
        transformed = transformed(valid_idx, :);

        % Skip cells that fall entirely outside the FOV
        if isempty(transformed)
            continue;
        end
        valid_cells(i) = true;

        % Round and clip to valid integer indices
        transformed = round(transformed);
        transformed(transformed < 1) = 1; % Ensure indices are >= 1
        transformed(transformed(:, 1) > 2048, 1) = 2048; % Clip row indices
        transformed(transformed(:, 2) > 2048, 2) = 2048; % Clip column indices

        % Populate the transformed mask with morphological dilation
        transformed_mask = zeros(2048, 2048);
        for j = 1:size(transformed, 1)
            transformed_mask(transformed(j, 1), transformed(j, 2)) = 1;
        end

        % Apply morphological dilation to fill gaps in the mask
        transformed_mask = imdilate(transformed_mask, strel('disk', 2));

        % Add the transformed mask to the result
        S_p1_transformed = cat(3, S_p1_transformed, transformed_mask);

        % Check if this cell matches the "3rd" plane
        if valid_cells(i)
            matchedCellIDs = [matchedCellIDs, i];
        end
    end

    % Display valid cells
    disp('Valid Cells:');
    disp(find(valid_cells));
end

%% Helper Function: Extract Centroids
function centroids = extractCentroids(binaryStack)
    % Extract centroids of non-zero regions from a binary stack
    % Inputs:
    %   - binaryStack: 3D binary matrix representing cell masks
    % Outputs:
    %   - centroids: [N x 2] matrix of (x, y) centroids

    centroids = [];
    for i = 1:size(binaryStack, 3)
        stats = regionprops(binaryStack(:, :, i), 'Centroid');
        if ~isempty(stats)
            centroids = [centroids; cat(1, stats.Centroid)]; %#ok<AGROW>
        end
    end
end