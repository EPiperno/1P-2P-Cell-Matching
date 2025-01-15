scriptpath =  "D:\Enrico\Alignment";
addpath(scriptpath);
filepath = "D:\Enrico\Alignment\Data\20241118\meas03\";
save_filepath = filepath + "Outputs\";

p1_filename = "m3_d241118_s03_1p";
p2_filename = "m3_d241118_s03_2p_00001";
p2_n_planes = 10;

%% Load Pre Processed Data
transformation_filename = save_filepath + p1_filename + "_pre_processed.mat";
load(transformation_filename, 'S_p1_final', 'S_p2_final', 'T_p1_final', 'T_p2_final', 'cellIDs');

%% Call Match Cells Function
matched_cells = match_cells_1P_2P(T_p1_final, S_p1_final, T_p2_final, S_p2_final);

%% Plot Matched Cells
plot_matched_cells(matched_cells, S_p1_final, S_p2_final, true)

%% Plot Z Plane Distribution
plot_z_plane_distribution(matched_cells)