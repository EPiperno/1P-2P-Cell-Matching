function S_p2_matched_resolution = matchResolution(S_p2, target_size)
    % matchResolution: Scales up S_p2 to match the resolution of S_p1.
    %
    % Inputs:
    %   - S_p2: 1x10 cell array of 512x512xcell_count matrices for 2P data.
    %   - target_size: [rows, cols] target resolution (e.g., [2048, 2048]).
    %
    % Output:
    %   - S_p2_matched_resolution: 1x10 cell array with resolution matched to S_p1.

    % Initialize output cell array
    S_p2_matched_resolution = cell(size(S_p2));

    % Iterate over each plane in S_p2
    for z = 1:length(S_p2)
        if ~isempty(S_p2{z})
            [rows, cols, cell_count] = size(S_p2{z});
            if rows ~= target_size(1) || cols ~= target_size(2)
                % Create an empty array for scaled data
                scaled_data = zeros(target_size(1), target_size(2), cell_count);

                % Scale each plane in S_p2{z} to the target resolution
                for i = 1:cell_count
                    scaled_data(:, :, i) = imresize(S_p2{z}(:, :, i), target_size, 'nearest');
                end

                % Assign the scaled data to the output cell array
                S_p2_matched_resolution{z} = scaled_data;
            else
                % If already at target resolution, copy data as is
                S_p2_matched_resolution{z} = S_p2{z};
            end
        end
    end
end