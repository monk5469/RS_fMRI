%% search_optimal_ROI_centers.m
% This script identifies optimal ROI centers (rOFA and rFFA) by maximizing 
% the number of activated voxels within a spherical search radius.
%
% Requirements:
%   - SPM12 must be installed and on the MATLAB path.
%   - 'workspace_data.mat' must be located in the parent folder of 'dirname'.
%   - Each subject's activation results should be stored under:
%       target_dir{s}/bold/<stat_dir>/SPM.mat
%
% Outputs:
%   - Optimal ROI center coordinates printed to console
%   - Exported coverage tables: 'coverage_OFA.csv', 'coverage_FFA.csv'

clear;

%% Set paths
dirname = uigetdir(pwd, 'Select project directory');
load(fullfile(fileparts(dirname), 'workspace_data.mat'));  % Load necessary workspace data

stat_dir = 'Imm_Del';  % Contrast folder name
contrast_idx = 2;       % Select contrast index (e.g., Imm > Del)

%% Step 1: Extract peak coordinates from SPM result table
coordinates = cell(length(objects_dir), 1);

for s = 1:length(objects_dir)
    data_path = fullfile(target_dir{s}, 'bold', stat_dir);
    load(fullfile(data_path, 'SPM.mat'));

    xSPM = struct();
    xSPM.swd = SPM.swd;
    xSPM.title = SPM.xCon(contrast_idx).name;
    xSPM.Ic = contrast_idx;
    xSPM.u = 0.05;            % Uncorrected p-value threshold
    xSPM.k = 0;              % Cluster extent threshold
    xSPM.thresDesc = 'none'; % No FWE/FDR correction
    xSPM.Im = [];
    xSPM.pm = [];

    [hReg, xSPM, SPM] = spm_results_ui('setup', xSPM);
    TabDat = spm_list('List', xSPM, hReg);
    coordinates{s} = TabDat.dat(:, end);  % Extract MNI coordinates
end

%% Step 2: Define anatomical ranges for rOFA and rFFA (based on prior literature)
rOFA_x = [30, 40]; rOFA_y = [-95, -70]; rOFA_z = [-20, 0];
rFFA_x = [30, 40]; rFFA_y = [-65, -35]; rFFA_z = [-25, -5];

All_OFA = [];
All_FFA = [];

for s = 1:length(objects_dir)
    coords = coordinates{s};
    OFA_coords = [];
    FFA_coords = [];

    for i = 1:length(coords)
        x = coords{i}(1);
        y = coords{i}(2);
        z = coords{i}(3);

        if x >= rOFA_x(1) && x <= rOFA_x(2) && ...
           y >= rOFA_y(1) && y <= rOFA_y(2) && ...
           z >= rOFA_z(1) && z <= rOFA_z(2)
            OFA_coords = [OFA_coords; x, y, z];
        end

        if x >= rFFA_x(1) && x <= rFFA_x(2) && ...
           y >= rFFA_y(1) && y <= rFFA_y(2) && ...
           z >= rFFA_z(1) && z <= rFFA_z(2)
            FFA_coords = [FFA_coords; x, y, z];
        end
    end

    All_OFA = [All_OFA; OFA_coords];
    All_FFA = [All_FFA; FFA_coords];
end

All_OFA = unique(All_OFA, 'rows');
All_FFA = unique(All_FFA, 'rows');

%% Step 3: Search for optimal ROI center using spherical coverage
search_radius = 10;  % mm

[optimal_OFA, coverage_OFA, candidate_OFA] = find_optimal_center(All_OFA, rOFA_x, rOFA_y, rOFA_z, search_radius);
[optimal_FFA, coverage_FFA, candidate_FFA] = find_optimal_center(All_FFA, rFFA_x, rFFA_y, rFFA_z, search_radius);

fprintf('Optimal rOFA Center: (%.2f, %.2f, %.2f)\n', optimal_OFA);
fprintf('Optimal rFFA Center: (%.2f, %.2f, %.2f)\n', optimal_FFA);

%% Step 4: Export candidate centers and their coverage
T_OFA = table(candidate_OFA(:,1), candidate_OFA(:,2), candidate_OFA(:,3), coverage_OFA(:), ...
              'VariableNames', {'x', 'y', 'z', 'coverage'});
T_FFA = table(candidate_FFA(:,1), candidate_FFA(:,2), candidate_FFA(:,3), coverage_FFA(:), ...
              'VariableNames', {'x', 'y', 'z', 'coverage'});

writetable(T_OFA, fullfile(dirname, 'coverage_OFA.csv'));
writetable(T_FFA, fullfile(dirname, 'coverage_FFA.csv'));

%% ========== FUNCTION ==========

function [optimal_center, coverage, candidate_centers] = find_optimal_center(points, x_range, y_range, z_range, radius)
% Find the optimal ROI center that maximizes voxel coverage within a sphere
% Inputs:
%   - points: Nx3 matrix of voxel coordinates
%   - *_range: two-element vector defining search bounds for x/y/z
%   - radius: search radius in mm
% Outputs:
%   - optimal_center: 1x3 coordinate of optimal center
%   - coverage: vector of voxel counts for each candidate center
%   - candidate_centers: list of all candidate center coordinates

    if isempty(points)
        optimal_center = NaN(1,3);
        coverage = [];
        candidate_centers = [];
        return;
    end

    [X, Y, Z] = ndgrid(x_range(1):2:x_range(2), ...
                       y_range(1):2:y_range(2), ...
                       z_range(1):2:z_range(2));
    candidate_centers = [X(:), Y(:), Z(:)];

    coverage = zeros(size(candidate_centers,1),1);
    max_coverage = 0;
    optimal_center = mean(points, 1, 'omitnan');

    for i = 1:size(candidate_centers, 1)
        center = candidate_centers(i, :);
        distances = sqrt(sum((points - center).^2, 2));
        coverage(i) = sum(distances <= radius);

        if coverage(i) > max_coverage
            max_coverage = coverage(i);
            optimal_center = center;
        end
    end
end
