function output = save_extract(data, filename)
    fprintf('Computing Motion Correction...')
    % %% Compute Motion Correction
    % addpath('D:\Enrico\Alignment\NoRMCorre');
    % 
    % Y = data;
    % % tic; Y = read_file(name); toc; % read the file (optional, you can also pass the path in the function instead of Y)
    % Y = single(Y);                 % convert to single precision 
    % T = size(Y,ndims(Y));
    % Y = Y - min(Y(:));
    % %% Set parameters (first try out rigid motion correction)
    % options_rigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'bin_width',200,'max_shift',15,'us_fac',50,'init_batch',200);
    % 
    % % perform rigid motion correction
    % tic; [M1,shifts1,template1,options_rigid] = normcorre(Y,options_rigid); toc
    % fprintf('Computed!')

    fprintf('Computing Extract...')
    %% Compute Extract
    config=[];
    config = get_defaults(config); 
    config.avg_cell_radius=5;
    config.trace_output_option='no_constraint';
    config.num_partitions_x=1;
    config.num_partitions_y=1; 
    config.use_gpu=1; 
    config.max_iter = 5; 
    config.cellfind_min_snr=0;
    config.thresholds.T_min_snr=3;
    output=extractor(data,config);
    %% Save Extract Results
    extract_filename = filename + "_extract_output.mat";
    save(extract_filename, 'output', '-v7.3');
    fprintf('Computed!')
end
