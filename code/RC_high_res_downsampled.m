

%% High-Res Downsampled Model Simulation (Board E-field Stimulation Version)

% Description: This script simulates an RC network based on the high resolution
% downsampled structural connectome, using downsampled connectivity (SC),
% distance (D), and voltage (V) matrices to model TMS-driven signal propagation.

%% Adjustable Parameters

a = 1;                     % Subject ID or session ID for file naming
namestr = 'Downsampled_TIK'; % Base name for output files
time_scale = 50;            % Total simulation time
num_points = 100;          % Number of points for time resolution

%% Load Input Files

% Load downsampled structural connectome and distance matrices
load('SC_downsampled.mat');  % Downsampled SC matrix
load('D_downsampled.mat');   % Downsampled distance matrix
load('E_downsampled.mat');   % Downsampled initial voltage vector
load('group_indices.mat');   % Downsampled grouping information
load('nodes_fmri_values.mat');
% fMRI_data_full = fMRI_values;
% threshold = 59412;  % Set the threshold value
% filtered_groups = {};  % Initialize a cell array to store filtered groups
% filtered_centroids = [];
% fMRI_downsampled = [];
%fMRI_downsampled = zeros(num_groups, 1);
% for group_id = 1:numel(group_indices)
%     % Get the current group
%     current_group = group_indices{group_id};
% 
%     % Check if all indices in the group satisfy the condition
%     if all(current_group <= threshold)
%         % If condition is satisfied, add the group to the filtered list
%         filtered_groups{end+1} = current_group;
%         filtered_centroids{end+1} = coords(group_id,:);
%     end
% end
% num_groups = size(filtered_groups,2);
% for i = 1:num_groups
%     nodes_in_group_i = filtered_groups{i};
%     fMRI_downsampled(i) = mean(abs(fMRI_data_full(nodes_in_group_i))); % Use 'sum' if required
% end
% target_data = fMRI_downsampled;
%target_data(1,filtered_selected_indices) = 0;

% % Normalize to unit sphere
% coords = coords ./ vecnorm(coords, 2, 2);
% cortical_coords = coords(1:998, :);

% Load the 91282x3 position coordinates
position_coordinates = readmatrix('positionc.csv');  % Replace with your file name

% Assume `group_indices` is a cell array where each cell contains node indices
num_groups = numel(group_indices);

% Initialize the downsampled coordinates matrix
downsampled_positions = zeros(num_groups, 3);

% Loop through each group
for group_id = 1:num_groups
    % Get the indices for the current group
    nodes_in_group = group_indices{group_id};

    % Extract the corresponding positions
    group_positions = position_coordinates(nodes_in_group, :);

    % Compute the centroid (mean position) for the group
    downsampled_positions(group_id, :) = mean(group_positions, 1);
end
% % Load MNI-space coordinates
% group_centroids_mat = load('group_centroids.mat');
% group_centroids = group_centroids_mat.group_centroids;
% Normalize to unit sphere
%coords = downsampled_positions;
coords = downsampled_positions ./ vecnorm(downsampled_positions, 2, 2);


%% Data Preparation

% Define the resistance matrix from the downsampled SC matrix
%R_downsampled = D_downsampled;
R_downsampled = D_downsampled ./ SC_downsampled; % Convert SC to resistance (R = D / SC)
R_downsampled(isnan(R_downsampled) | isinf(R_downsampled)) = 0; % Handle NaNs/Infs

% % Normalize the SC matrix (optional, as done in Atlas model)
% SC = SC_downsampled ./ max(SC_downsampled(:)) * 1000;
% %%
% V = zeros(size(E_downsampled,1),1);
% initial_voltage = 1;
% V(766) = initial_voltage;
%%
% Set the initial voltage from the downsampled V matrix
%V = E_downsampled;


n_weight = 1;
weight = 1;
V = zeros(size(E_downsampled,1),1);
% Find the node with the maximum voltage in V
[max_val, max_idx] = max(E_downsampled);

% % Create a new voltage vector with only the largest node retained
% V_new = zeros(size(V));
% V_new(max_idx) = max_val;  % Retain the maximum voltage at its original index
%
% % Update V to V_new
% V = V_new;
% %% Apply 90% Threshold to Voltage (V)
% % Calculate the 90th percentile of the voltage matrix V across all time steps
% threshold_V = prctile(V(:), 98.8);  % 90th percentile of voltage values
% 0
% % Apply the threshold: Set values below the threshold to zero
% V(V < threshold_V) = 0;  % Replace values below the threshold

% % Define the percentage of top values to retain
% retain_percentage = 0.02;  % Set the percentage of values to retain (e.g., top 10%)
%
% % Calculate the threshold corresponding to the top retain_percentage
% threshold_V = prctile(V(:), 100 - retain_percentage);  % Retain the top X%
%weight = 1;
V(max_idx) = 2*max_val;
%V(1468) =  max_val;re
% % Apply the threshold: Set values below the threshold to zero
% V(V < threshold_V) = 0;  % Replace values below the threshold with zero
%V(565) =  weight*max_val;
 %V(1912) = max_val;    %(single: position: TMS: max_idx = 565, L M1:1468 L Visual: 2072 R SGC:5578)
 %V(5323) = max_val;
% Define capacitance matrix (uniform across nodes)
n = numel(V);  % Number of nodes in downsampled model
%C = ones(n, 1) * 1000;
num_groups = numel(group_indices);
C = zeros(num_groups, 1);

% Loop through each group and assign capacitance based on size
for group_id = 1:num_groups
    % Get the size of the current group
    group_size = numel(group_indices{group_id});

    % Assign capacitance (proportional to group size)
    C(group_id) = group_size * 100;  % Example: directly proportional


end

%% Simulation Setup

% Set up and run the ODE solver with the downsampled RC circuit model
sol = ode23(@(t, v) capODE_downsampled(t, v, R_downsampled, C), [0 time_scale], V);

% Define time points for results
t = linspace(0, time_scale, num_points);
[v, vd] = deval(sol, t);

%% Calculate Currents and Power for Each Node

Ip = zeros(n, num_points);
In = zeros(n, num_points);

for i = 1:n
    for j = 1:n
        if i ~= j && R_downsampled(j, i) ~= 0
            for k = 1:num_points
                current_val = (v(j, k) - v(i, k))  ./ R_downsampled(j, i) ./ 2;
                if v(j, k) >= v(i, k)
                    Ip(i, k) = Ip(i, k) + current_val;
                else
                    In(i, k) = In(i, k) + current_val;
                end
            end
        end
    end
end

% Calculate total current, power metrics, and max values over time
I = In + Ip;
Cp = cumtrapz(abs(Ip), 2);
E1 = abs(I .* v);
E1c = cumtrapz(E1, 2);
Cpmax = max(Cp, [], 2);
E1max = max(E1c, [], 2);
Imax = max(I, [], 2);
Ipmax = max(Ip, [], 2);
%Imax(filtered_selected_indices) = 0;
%[r(n_weight),p(n_weight)] = corr(Imax(1:6192),target_data',"Type","Spearman");

%% Save Results to Matrix Files

output_folder = './ODESolver/';
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

save(fullfile(output_folder, 'voltage.mat'), 'v', 't');
save(fullfile(output_folder, 'I.mat'), 'I');
save(fullfile(output_folder, 'Ip.mat'), 'Ip');
save(fullfile(output_folder, 'Imax.mat'), 'Imax');
save(fullfile(output_folder, 'Ipmax.mat'), 'Ipmax');
save(fullfile(output_folder, 'E1max.mat'), 'E1max');
save(fullfile(output_folder, 'Cpmax.mat'), 'Cpmax');

% Total number of original nodes
num_original_nodes = 91282;

% Initialize arrays for recovered data
I_full = zeros(num_original_nodes, size(I, 2));    % Time-dependent current
Ip_full = zeros(num_original_nodes, size(Ip, 2)); % Time-dependent positive current
E1_full = zeros(num_original_nodes, size(E1, 2)); % Time-dependent energy
Cp_full = zeros(num_original_nodes, size(Cp, 2)); % Time-dependent positive charge accumulated
Imax_full = zeros(num_original_nodes, 1);         % Maximum current
Ipmax_full = zeros(num_original_nodes, 1);        % Maximum positive current
E1max_full = zeros(num_original_nodes, 1);        % Maximum energy
Cpmax_full = zeros(num_original_nodes, 1);        % Maximum positive charge accumulated
%% Expand Downsampled Data to Original Resolution
for group_id = 1:numel(group_indices)
    % Get the indices of nodes in this group
    nodes_in_group = group_indices{group_id};

    % Assign downsampled values to all nodes in the group
    I_full(nodes_in_group, :) = repmat(I(group_id, :), numel(nodes_in_group), 1);
    Ip_full(nodes_in_group, :) = repmat(Ip(group_id, :), numel(nodes_in_group), 1);
    E1_full(nodes_in_group, :) = repmat(E1(group_id, :), numel(nodes_in_group), 1);
    Cp_full(nodes_in_group, :) = repmat(Cp(group_id, :), numel(nodes_in_group), 1);
    Imax_full(nodes_in_group) = Imax(group_id);
    Ipmax_full(nodes_in_group) = Ipmax(group_id);
    E1max_full(nodes_in_group) = E1max(group_id);
    Cpmax_full(nodes_in_group) = Cpmax(group_id);
end


save(fullfile(output_folder, 'Imax_full.mat'), 'Imax_full');
save(fullfile(output_folder, 'Ipmax_full.mat'), 'Ipmax_full');
save(fullfile(output_folder, 'E1max_full.mat'), 'E1max_full');
save(fullfile(output_folder, 'Cpmax_full.mat'), 'Cpmax_full');

% Generate figures and save them as dscalar files
c = cifti_read('ones.dscalar.nii');
clabel0 = c.cdata;

% Time-dependent metrics: I, Ip, E1, Cp
metrics_time_dependent = {I_full, Ip_full, E1_full, Cp_full};
metric_names_time_dependent = {'I', 'Ip', 'E', 'Cp'};

% Log-scaled maximum metrics: Imax, Ipmax, Emax, Cpmax
metrics_max = {Imax_full, Ipmax_full, E1max_full, Cpmax_full};
metric_names_max = {'Imax', 'Ipmax', 'Emax', 'Cpmax'};
%
% Plot time-dependent metrics
for i = 1:length(metrics_time_dependent)
    metric_data = metrics_time_dependent{i};
    for t_idx = 1:num_points
        c.cdata = metric_data(:, t_idx);  % Map the data for this time step
        filename_var = sprintf('%s_time_%d.dscalar.nii', metric_names_time_dependent{i}, t_idx);
        ciftisavereset(c, filename_var);
        fprintf('Saved %s at time step %d\n', metric_names_time_dependent{i}, t_idx);
    end
end
% 
% % Plot log-scaled maximum metrics
% for i = 1:length(metrics_max)
%     metric_data = metrics_max{i};
% 
% 
% 
%     log_metric_data = log10(metric_data + 1);  % Apply log scale
%     c.cdata = log_metric_data;  % Map the log-scaled data
%     filename_var = sprintf('%s_log_scaled_%s.dscalar.nii', metric_names_max{i}, namestr);
%     ciftisavereset(c, filename_var);
%     fprintf('Saved log-scaled %s\n', metric_names_max{i});
% end

% 
% c = cifti_read('ones.dscalar.nii');
% tstep = 1;
% clabel0 = c.cdata;
% n = size(clabel0,1);
% clabel1 = zeros(n,tstep);
% V_full = zeros(num_original_nodes, 1);
% for group_id = 1:numel(group_indices)
%     nodes_in_group = group_indices{group_id};
% 
%     V_full(nodes_in_group) = V(group_id);
% end
% c.cdata = V_full(:);
% filename_var = 'Filter_V.dscalar.nii';
% ciftisavereset(c,filename_var);


%% Plotting to dscalar Files for Connectome Workbench

% Load the Schaefer atlas file to map results to brain regions (or a similar atlas file)
c = cifti_read('ones.dscalar.nii');
%c = cifti_read('Schaefer2018_1000Parcels_7Networks_order.dlabel.nii');
clabel0 = c.cdata;

% Helper function to save each metric to a dscalar file
save_dscalar = @(data, label, suffix) ...
    ciftisavereset(setfield(c, 'cdata', data),  sprintf('%s_SimNIBS_%s.dscalar.nii', label, suffix));

% Mapping function to match data points with atlas regions using the reordered indices
function mapped_data = map_to_cifti(metric, clabel0, ReOrder, n)
    mapped_data = NaN(size(clabel0));
    for i = 1:numel(clabel0)
        if clabel0(i) > 0 && clabel0(i) <= n
            mapped_data(i) = metric(ReOrder(clabel0(i)));
        end
    end
end
% 
% % Plot and save results as dscalar files for Cpmax, Imax, Emax, and Ipmax
% save_dscalar(map_to_cifti(log10(Cpmax + 1), clabel0, ReOrder, n), namestr, 'Cpmax');
% save_dscalar(map_to_cifti(log10(Imax + 1), clabel0, ReOrder, n), namestr, 'Imax');
% save_dscalar(map_to_cifti(log10(E1max + 1), clabel0, ReOrder, n), namestr, 'Emax');
% save_dscalar(map_to_cifti(log10(Ipmax + 1), clabel0, ReOrder, n), namestr, 'Ipmax');

disp('Simulation and output generation complete.');

%% Function: capODE_downsampled with SC and C as Inputs

function dvdt = capODE_downsampled(t, v, R, C)
    % Differential equations for RC circuit in downsampled model
    num_nodes = length(v);  
    dvdt = zeros(num_nodes, 1);  

    % Compute rate of voltage change at each node
    for i = 1:num_nodes
        for j = 1:num_nodes
            if i ~= j && R(i, j) ~= 0
                dvdt(i) = dvdt(i) + (v(j) - v(i)) / (R(i, j) * C(i));
            end
        end
    end
end



load('nodes_fmri_values.mat');
load('filtered_selected_indices.mat');
fMRI_data_full = fMRI_values;
threshold = 59412;  % Set the threshold value
filtered_groups = {};  % Initialize a cell array to store filtered groups
filtered_centroids = [];
fMRI_downsampled = [];
%fMRI_downsampled = zeros(num_groups, 1);
for group_id = 1:numel(group_indices)
    % Get the current group
    current_group = group_indices{group_id};

    % Check if all indices in the group satisfy the condition
    if all(current_group <= threshold)
        % If condition is satisfied, add the group to the filtered list
        filtered_groups{end+1} = current_group;
        filtered_centroids{end+1} = coords(group_id,:);
    end
end
num_groups = size(filtered_groups,2);
for i = 1:num_groups
    nodes_in_group_i = filtered_groups{i};
    fMRI_downsampled(i) = mean(abs(fMRI_data_full(nodes_in_group_i))); % Use 'sum' if required
end
target_data = fMRI_downsampled;
target_data(1,filtered_selected_indices) = 0;
    
target_data1 = target_data(1:6192)';


% 
% %% Scatter Plot for Spearman Correlation 
% Imax(filtered_selected_indices) = 0;
% Ipmax(filtered_selected_indices) = 0;
% E1max(filtered_selected_indices) = 0;
% Cpmax(filtered_selected_indices) = 0;
% Imax(max_idx) = 0;
% Imax(1912) = 0;Imax(5323) = 0;
% metrics = {Imax(1:6192), Ipmax(1:6192), Cpmax(1:6192), E1max(1:6192)};
% metric_names = {'Imax', 'Ipmax', 'Cpmax', 'Emax'};
% 
% 
% 
% z_threshold = 2;  % Define the z-score threshold
% 
% for metric_idx = 1:4
%     %subplot(2, 2, metric_idx);
%     figure;
%     % Extract metric and fMRI signal
%     metric_data = metrics{metric_idx};
%     metric_data = log1p(metric_data);  % log(1 + x) avoids log(0)
%     fMRI_signal = target_data1;  % Assuming a2 contains fMRI signal for 998 nodes
% 
%     % Calculate z-scores for the metric and fMRI signal
%     z_metric = zscore(metric_data);
%     z_fMRI = zscore(fMRI_signal);
% 
%     % Identify nodes with z-scores within the threshold
%     valid_idx = abs(z_metric) <= z_threshold & abs(z_fMRI) <= z_threshold;
% 
%     % Filter data based on valid indices
%     metric_filtered = metric_data(valid_idx);
%     fMRI_filtered = fMRI_signal(valid_idx);
% 
%     % Scatter plot
%     scatter(fMRI_filtered, (metric_filtered), 'b', 'filled', 'MarkerFaceAlpha', 0.5);
%     hold on;
% 
%     % Nonparametric LOESS fit for visualization
%     %smooth_fit = smooth(fMRI_filtered, metric_filtered, 0.2, 'lowess');  % Adjust span (0.2) as needed
%     %plot(fMRI_filtered, smooth_fit, 'r-', 'LineWidth', 2);  % Plot the LOESS curve
% 
%     % Spearman correlation
%     [rho, pval] = corr( metric_data,fMRI_signal, 'Type', 'Spearman');
% 
%     % Add title and annotations
%     %title(sprintf('Scatter Plot with Spearman Correlation: %s', metric_names{metric_idx}), 'FontSize', 20);
%     xlabel('fMRI Signal', 'FontSize', 24);
%     ylabel(['Max Change (', metric_names{metric_idx}, ')'], 'FontSize', 24);
%     legend('Data', 'Location', 'best');
%     text(min(fMRI_filtered), max(metric_filtered), ...
%          sprintf('r = %.2f, p = %.3f', rho, pval), ...
%          'VerticalAlignment', 'top', 'FontSize', 24);
%     set(gca, 'FontSize', 24); % Adjust '14' to your preferred size
%     grid on;
%     hold off;
% end
% 
% metrics = {Imax(1:6192), Ipmax(1:6192), Cpmax(1:6192), E1max(1:6192)};
% metric_names = {'Imax', 'Ipmax', 'Cpmax', 'Emax'};
% figure;
% 
% z_threshold = 2;  % Define the z-score threshold
% 
% for metric_idx = 1:1
%     % Extract metric and fMRI signal
%     metric_data = metrics{metric_idx};
%     fMRI_signal = target_data1;  % Assuming target_data1 contains fMRI signal for 998 nodes
% 
%     % Calculate z-scores for the metric and fMRI signal
%     z_metric = zscore(metric_data);
%     z_fMRI = zscore(fMRI_signal);
% 
%     % Identify nodes with z-scores within the threshold
%     valid_idx = abs(z_metric) <= z_threshold & abs(z_fMRI) <= z_threshold;
% 
%     % Filter data based on valid indices
%     metric_filtered = metric_data(valid_idx);
%     fMRI_filtered = fMRI_signal(valid_idx);
% 
%     % Scatter plot
%     scatter(fMRI_filtered, metric_filtered, 'b', 'filled', 'MarkerFaceAlpha', 0.5);
%     hold on;
% 
%     % Linear fit (least-squares regression)
%     p = polyfit(fMRI_filtered, metric_filtered, 1);  % Linear fit: y = p(1)*x + p(2)
%     linear_fit = polyval(p, fMRI_filtered);
% 
%     % Plot the linear fit line
%     plot(fMRI_filtered, linear_fit, 'r-', 'LineWidth', 2);
% 
%     % Spearman correlation
%     [rho, pval] = corr(metric_filtered, fMRI_filtered, 'Type', 'Spearman');
% 
%     % Add title and annotations
%     xlabel('fMRI Signal', 'FontSize', 24);
%     ylabel(['Max Change (', metric_names{metric_idx}, ')'], 'FontSize', 24);
%     legend('Data', 'Linear Fit', 'Location', 'best');
%     text(min(fMRI_filtered), max(metric_filtered), ...
%          sprintf('r = %.2f, p = %.3f', rho, pval), ...
%          'VerticalAlignment', 'top', 'FontSize', 24);
%     set(gca, 'FontSize', 24); % Adjust '14' to your preferred size
%     grid on;
%     hold off;
% end
% % 
% % Overall plot adjustments
% sgtitle('Scatter Plots of Metrics vs. fMRI Signal (Using LOESS Fit)', 'FontSize', 24);


% 
% load('mni_nodes_fmri_values.mat');
% fMRI_data_full = fMRI_values;
% %% Filter Groups with All Indices <= 64984
% threshold = 64984;  % Set the threshold value
% filtered_groups = {};  % Initialize a cell array to store filtered groups
% 
% %fMRI_downsampled = zeros(num_groups, 1);
% for group_id = 1:numel(group_indices)
%     % Get the current group
%     current_group = group_indices{group_id};
% 
%     % Check if all indices in the group satisfy the condition
%     if all(current_group <= threshold)
%         % If condition is satisfied, add the group to the filtered list
%         filtered_groups{end+1} = current_group;
%     end
% end
% fMRI_downsampled = [];
% num_groups = size(group_indices,2);
% for i = 1:num_groups
%     nodes_in_group_i = group_indices{i};
%     fMRI_downsampled(i) = mean(abs(fMRI_data_full(nodes_in_group_i))); % Use 'sum' if required
% end
% target_data = fMRI_downsampled;
% target_data(1,filtered_selected_indices) = 0;
% 
% target_data = fMRI_downsampled(1:size(group_indices,2))';
% cortical_coords = coords(1:size(group_indices,2),:);
% 
% %% Load Metrics
% % Example metrics (replace with actual loaded data)
% metrics = {Imax(1:num_groups), Ipmax(1:num_groups), Cpmax(1:num_groups), E1max(1:num_groups)}; % Replace with your data (e.g., size 998 x 1)
% metric_names = {'Imax', 'Ipmax', 'Cpmax', 'E1max'}; % Names of metrics
% %% Observed Correlations
% observed_correlations = zeros(length(metrics), 1);
% for metric_idx = 1:length(metrics)
%     metric_data = metrics{metric_idx};
%     [observed_correlations(metric_idx), ~] = corr(metric_data, target_data, 'Type', 'Spearman');
% end
% 
% %% Spin Permutations
% num_permutations = 100; % Number of permutations
% null_correlations = zeros(length(metrics), num_permutations);
% rotated_coords = [];
% for perm_idx = 1:num_permutations
%     % Generate a random rotation matrix
%     [U, ~, V] = svd(randn(3)); % Orthogonal random rotation
%     rotation_matrix = U * diag([1, 1, det(U * V')]) * V';
% 
%     % Apply rotation to cortical coordinates
%     rotated_coords = (rotation_matrix * cortical_coords')';
% 
%     % Map rotated coordinates back to original node indices
%     permuted_indices = assign_closest_vertices(cortical_coords, rotated_coords);
% 
%     % Compute null correlations for each metric
%     for metric_idx = 1:length(metrics)
%         metric_data = metrics{metric_idx};
%         permuted_metric = metric_data(permuted_indices);
%         null_correlations(metric_idx, perm_idx) = corr(permuted_metric, target_data, 'Type', 'Spearman');
%     end
% end
% 
% %% Compute Empirical P-Values
% p_values_spin = zeros(length(metrics), 1);
% for metric_idx = 1:length(metrics)
%     observed_corr = observed_correlations(metric_idx);
%     null_dist = null_correlations(metric_idx, :);
%     % Two-tailed p-value
%     p_values_spin(metric_idx) = mean((null_dist) >= (observed_corr));
% end
% 
% %% Save Results
% save('spin_test_results.mat', 'observed_correlations', 'null_correlations', 'p_values_spin');
% 
% %% Plot Results
% figure;
% for metric_idx = 1:length(metrics)
%     % Null distribution and observed correlation
%     null_dist = null_correlations(metric_idx, :);
%     observed_corr = observed_correlations(metric_idx);
%     spin_p_value = p_values_spin(metric_idx);
% 
%     % Plot histogram of null distribution
%     subplot(2, 2, metric_idx);
%     histogram(null_dist, 30, 'Normalization', 'probability');
%     hold on;
% 
%     % Plot observed correlation
%     xline(observed_corr, 'r', 'LineWidth', 2, 'DisplayName', 'Observed Correlation');
% 
%     % Add title with r and p-value
%     title(sprintf('%s\nr = %.4f, p = %.4f', metric_names{metric_idx}, observed_corr, spin_p_value), 'FontSize', 12);
%     xlabel('Correlation Coefficient', 'FontSize', 20);
%     ylabel('Probability', 'FontSize', 20);
%     grid on;
%     hold off;
% end
% 
% % Add overall title
% sgtitle('Spin Test Results: Null Distributions and Observed Correlations', 'FontSize', 20);
% 
% %% Helper Function: Assign Closest Vertices
% function assigned = assign_closest_vertices(original_coords, rotated_coords)
%     % Map rotated coordinates back to the closest original vertices
%     num_vertices = size(original_coords, 1);
%     assigned = zeros(num_vertices, 1);
%     for v = 1:num_vertices
%         distances = sqrt(sum((original_coords - rotated_coords(v, :)).^2, 2));
%         [~, assigned(v)] = min(distances);
%     end
% end
% 
% 
% 
% %% Inputs
% c3 = cifti_read('Schaefer2018_1000Parcels_7Networks_order.dscalar.nii');
% load('networks.mat');
% c3_data = c3.cdata;  % 91,282 nodes (each assigned to a Schaefer index)
% network_names = fieldnames(networks);  % Names of the 7 networks
% num_networks = length(network_names);
% 
% %% Initialize Storage
% nodes_by_network = cell(num_networks, 1);  % Store 91,282 node indices for each network
% 
% %% Map Nodes to Networks
% for n = 1:num_networks
%     % Get the Schaefer indices for the current network
%     network_indices = networks.(network_names{n});
% 
%     % Find all nodes in c3_data that match the Schaefer indices for this network
%     nodes_by_network{n} = find(ismember(c3_data, network_indices));
% end
% 
% %% Optional: Convert to a Structure for Easier Access
% nodes_by_network_struct = struct();
% for n = 1:num_networks
%     nodes_by_network_struct.(network_names{n}) = nodes_by_network{n};
% end
% %% Inputs
% % Example metrics (replace with your actual loaded data)
% metrics = {Imax_full, Ipmax_full, Cpmax_full, E1max_full}; % Size: 91,282 x 1
% metric_names = {'Imax', 'Ipmax', 'Cpmax', 'E1max'};        % Metric names
% 
% 
% 
% %% Initialize Storage
% num_networks = length(network_names);
% num_metrics = length(metrics);
% network_sums = zeros(num_networks, num_metrics); % Store sums for each metric by network
% 
% %% Sum Metrics for Each Network
% for n = 1:num_networks
%     % Get the indices of nodes belonging to the current network
%     node_indices = nodes_by_network{n};
% 
%     % Sum each metric for the current network
%     for m = 1:num_metrics
%         metric_data = metrics{m};
%         network_sums(n, m) = sum(metric_data(node_indices)); % Sum values for these nodes
%     end
% end
% 
% %% Save Results
% % Convert to table for easier visualization
% network_sums_table = array2table(network_sums, ...
%     'VariableNames', metric_names, ...
%     'RowNames', network_names);
% 
% % Save as .mat and .csv
% save('network_sums.mat', 'network_sums', 'network_sums_table');
% writetable(network_sums_table, 'network_sums.csv', 'WriteRowNames', true);
% 
% 
% %% Visualization: Separate Figures for Each Metric
% for m = 1:num_metrics
%     % Create a new figure for each metric
%     figure;
% 
%     % Plot the summed values for the current metric
%     bar(network_sums(:, m));
%     xticks(1:num_networks);
%     xticklabels(network_names);
%     xlabel('Networks', 'FontSize', 20);
%     ylabel(['Sum of ', metric_names{m}], 'FontSize', 20);
%     title(['Sum of ', metric_names{m}, ' Across Networks'], 'FontSize', 20);
%     grid on;
% 
%     % Save each figure
%     saveas(gcf, sprintf('Sum_%s_Across_Networks.fig', metric_names{m}));
%     saveas(gcf, sprintf('Sum_%s_Across_Networks.png', metric_names{m}));
% end
% 
% 
% 
% 
% 
% 
% %% Inputs
% c3 = cifti_read('Schaefer2018_1000Parcels_7Networks_order.dscalar.nii');
% load('Ro_networks.mat');
% c3_data = c3.cdata;  % 91,282 nodes (each assigned to a Schaefer index)
% network_names = fieldnames(Ro_networks);  % Names of the 7 networks
% num_networks = length(network_names);
% 
% %% Initialize Storage
% nodes_by_network = cell(num_networks, 1);  % Store 91,282 node indices for each network
% 
% %% Map Nodes to Networks
% for n = 1:num_networks
%     % Get the Schaefer indices for the current network
%     network_indices = Ro_networks.(network_names{n});
% 
%     % Find all nodes in c3_data that match the Schaefer indices for this network
%     nodes_by_network{n} = find(ismember(c3_data, network_indices));
% end
% 
% %% Optional: Convert to a Structure for Easier Access
% nodes_by_network_struct = struct();
% for n = 1:num_networks
%     nodes_by_network_struct.(network_names{n}) = nodes_by_network{n};
% end
% %% Inputs
% % Example metrics (replace with your actual loaded data)
% metrics = {Imax_full, Ipmax_full, Cpmax_full, E1max_full}; % Size: 91,282 x 1
% metric_names = {'Imax', 'Ipmax', 'Cpmax', 'E1max'};        % Metric names
% 
% 
% 
% %% Initialize Storage
% num_networks = length(network_names);
% num_metrics = length(metrics);
% network_sums = zeros(num_networks, num_metrics); % Store sums for each metric by network
% 
% %% Sum Metrics for Each Network
% for n = 1:num_networks
%     % Get the indices of nodes belonging to the current network
%     node_indices = nodes_by_network{n};
% 
%     % Sum each metric for the current network
%     for m = 1:num_metrics
%         metric_data = metrics{m};
%         network_sums(n, m) = sum(metric_data(node_indices)); % Sum values for these nodes
%     end
% end
% 
% %% Save Results
% % Convert to table for easier visualization
% network_sums_table = array2table(network_sums, ...
%     'VariableNames', metric_names, ...
%     'RowNames', network_names);
% 
% % Save as .mat and .csv
% save('network_sums.mat', 'network_sums', 'network_sums_table');
% writetable(network_sums_table, 'network_sums.csv', 'WriteRowNames', true);
% 
% %% Visualization: Separate Figures for Each Metric
% for m = 1:num_metrics
%     % Create a new figure for each metric
%     figure;
% 
%     % Plot the summed values for the current metric
%     bar(network_sums(:, m));
%     xticks(1:num_networks);
%     xticklabels(network_names);
%     xlabel('Networks', 'FontSize', 20);
%     ylabel(['Sum of ', metric_names{m}], 'FontSize', 20);
%     title(['Sum of ', metric_names{m}, ' Across Networks'], 'FontSize', 20);
%     grid on;
% 
%     % Save each figure
%     saveas(gcf, sprintf('Sum_%s_Across_Networks.fig', metric_names{m}));
%     saveas(gcf, sprintf('Sum_%s_Across_Networks.png', metric_names{m}));
% end
% 
% 
% 
% 

permno = 2000;
wsname = 'rotationFS_test.mat';
load('vertex_indices_L.mat');
load('vertex_indices_R.mat');
vertex_indices_L = vertex_indices_L + 1;
vertex_indices_R = vertex_indices_R + 1;

SL_num = 29696;
SR_num = 29716;
down_SL_num = 3091;
down_SR_num = 3101;

% Initialize file names and labels
data_files = {'Imax.mat', 'Ipmax.mat', 'E1max.mat', 'Cpmax.mat'};
data_labels = {'Imax', 'Ipmax', 'E1max', 'Cpmax'};

% Loop through each data file and plot results
for data_idx = 1:length(data_files)
    % Load data
    data_file = ['./ODESolver/', data_files{data_idx}];
    data_struct = load(data_file);
    data_field = data_labels{data_idx};
    data_values = data_struct.(data_field);


    data_values(filtered_selected_indices) = 0;

    % Split into left and right hemisphere data
    datal = data_values(1:down_SL_num);
    datar = data_values(down_SL_num+1:down_SL_num+down_SR_num);

    load('left_sphere_32k.mat');
    load('right_sphere_32k.mat');
    load('group_indices.mat');
    verticesl = verticesl(vertex_indices_L,:);
    verticesr = verticesr(vertex_indices_R,:);
    % %%extract the corresponding sphere surface coordinates for rotation
    % [verticesl, ~] = freesurfer_read_surf(fullfile(fshome,'subjects/fsaverage5/surf/lh.sphere'));
    % [verticesr, ~] = freesurfer_read_surf(fullfile(fshome,'subjects/fsaverage5/surf/rh.sphere'));

    load('group_indices.mat');   % Downsampled grouping information
    start_index_l = 1;
    end_index_l =  SL_num; % Equivalent to 59413
    start_index_r = 1+SL_num;
    end_index_r =  SR_num+SL_num; % Equivalent to 59413

    % Initialize an empty array to store the matching group indices
    matching_groups_l = {};
    matching_groups_r = {};
    % Loop through each group and check if it falls within the range
    for i = 1:length(group_indices)
        group = group_indices{i};

        % Check if all elements of the group are within the range
        if all(group >= start_index_l & group <= end_index_l)
            matching_groups_l{end+1} = group; % Append the group
        elseif all(group >= start_index_r & group <= end_index_r)
            matching_groups_r{end+1} = group; % Append the group
        end
    end



    num_groups = size(matching_groups_l,2);
    verticesld = zeros(size(matching_groups_l, 2), 3); % Preallocate for speed

    for i = 1:num_groups
        nodes_in_group_i = matching_groups_l{i};
        verticesld(i,:) = mean(verticesl(nodes_in_group_i,:),1); % Use 'sum' if required
    end
    num_groups = size(matching_groups_r,2);
    verticesrd = zeros(size(matching_groups_r, 2), 3); % Preallocate for speed
    for i = 1:num_groups
        nodes_in_group_i = matching_groups_r{i}-SL_num;
        verticesrd(i,:) = mean(verticesr(nodes_in_group_i,:),1); % Use 'sum' if required
    end

    load('nodes_fmri_values.mat');
    fMRI_data_full = fMRI_values;
    threshold = 59412;  % Set the threshold value
    filtered_groups = {};  % Initialize a cell array to store filtered groups
    filtered_centroids = [];
    fMRI_downsampled = [];
    %fMRI_downsampled = zeros(num_groups, 1);
    for group_id = 1:numel(group_indices)
        % Get the current group
        current_group = group_indices{group_id};

        % Check if all indices in the group satisfy the condition
        if all(current_group <= threshold)
            % If condition is satisfied, add the group to the filtered list
            filtered_groups{end+1} = current_group;
            filtered_centroids{end+1} = coords(group_id,:);
        end
    end
    num_groups = size(filtered_groups,2);
    for i = 1:num_groups
        nodes_in_group_i = filtered_groups{i};
        fMRI_downsampled(i) = mean(abs(fMRI_data_full(nodes_in_group_i))); % Use 'sum' if required
    end
    target_data = fMRI_downsampled;
    target_data(1,filtered_selected_indices) = 0;
    target_data1 = target_data(1:6192);
    % Calculate correlations
    rng(0);
    bigrotl = [];
    bigrotr = [];
    I1 = eye(3, 3);
    I1(1, 1) = -1;
    bl = verticesld;
    br = verticesrd;

    for j = 1:permno
        A = normrnd(0, 1, 3, 3);
        [TL, temp] = qr(A);
        TL = TL * diag(sign(diag(temp)));
        if (det(TL) < 0)
            TL(:, 1) = -TL(:, 1);
        end
        TR = I1 * TL * I1;
        bl = bl * TL;
        br = br * TR;

        Il = nearestneighbour(verticesld', bl');
        Ir = nearestneighbour(verticesrd', br');

        bigrotl = [bigrotl; datal(Il)'];
        bigrotr = [bigrotr; datar(Ir)'];
    end

    save(wsname, 'bigrotl', 'bigrotr');
    [r(data_idx), p(data_idx)] = corr(cat(1, datal, datar), target_data1', 'Type', 'Spearman');
    r_spin = zeros(permno, 4);
    for i = 1:permno
        [r_spin(i,data_idx), ~] = corr(cat(2, bigrotl(i, :), bigrotr(i, :))', target_data1', 'Type', 'Spearman');
    end

    observed_corr(data_idx) = r(data_idx);
    null_dist = r_spin(:,data_idx);
    p_values_spin(data_idx) = mean((null_dist) >= (observed_corr(data_idx)));

    % Plot null distribution
    figure;
    histogram(null_dist, 30, 'Normalization', 'probability');
    hold on;

    % Plot observed correlation
    xline(observed_corr(data_idx), 'r', 'LineWidth', 2, 'DisplayName', 'Observed Correlation');
    %title(sprintf('%s Null Distribution and Observed Correlation', data_field), 'FontSize', 24);
    xlabel('Correlation Coefficient', 'FontSize', 24);
    ylabel('Probability', 'FontSize', 24);
    subtitle(sprintf('r = %.3f, spin test p-value = %.3f', observed_corr(data_idx), p_values_spin(data_idx)),"FontSize",18);
    set(gca, 'FontSize', 24); % Adjust '14' to your preferred size
    grid on;
    hold off;
end
