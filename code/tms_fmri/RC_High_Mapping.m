

%% High-Res Downsampled Model Simulation (Board E-field Stimulation Version)

% Description: This script simulates an RC network based on the high resolution
% downsampled structural connectome, using downsampled connectivity (SC),
% distance (D), and voltage (V) matrices to model TMS-driven signal propagation.

%% Adjustable Parameters

a = 1;                     % Subject ID or session ID for file naming
namestr = 'Downsampled_TIK'; % Base name for output files
time_scale = 100;            % Total simulation time
num_points = 100;          % Number of points for time resolution
stim_freq = 100;                % Stimulation frequency (Hz)
num_pulses = 3;                % Number of pulses
pulse_interval = 1000 / stim_freq; % Interval between pulses (ms)
stim_times = 0:pulse_interval:(pulse_interval*(num_pulses-1)); % [0,100,200] ms
stim_voltage = 1;              % Stimulation voltage (1V)
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
stim_voltage = 2*max_val;             % Stimulation voltage
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
%sol = ode23(@(t, v) capODE_downsampled(t, v, R_downsampled, C), [0 time_scale], V);
% %% Simulation Setup with repeated stimulation
% sol = ode23(@(t, v) capODE_downsampled_rTMS(t, v, R_downsampled, C, stim_times, stim_voltage, max_idx), [0 time_scale], V);
%stim_times = [0, 100, 200]; % ms, times of pulses

% Define options with events to trigger pulses
options = odeset('Events', @(t,v)stim_events(t,v,stim_times));

% Initialize variables
tspan = [0 time_scale];
t_all = [];
v_all = [];
V_current = V; % initial condition


for pulse_idx = 1:length(stim_times)

    % Solve ODE up to next stimulation time
    if pulse_idx < length(stim_times)
        t_end = stim_times(pulse_idx+1);
    else
        t_end = time_scale;
    end

    % Run ODE solver from current time to next stimulation time
    [t_sol, v_sol] = ode23(@(t, v) capODE_downsampled(t, v, R_downsampled, C), [stim_times(pulse_idx), t_end], V_current, options);

    % Append solution
    t_all = [t_all; t_sol];
    v_all = [v_all, v_sol'];

    % Apply stimulation at end of current interval (except last)
    if pulse_idx < length(stim_times)
        V_current = v_sol(end,:)';
        V_current(max_idx) = V_current(max_idx) + stim_voltage;
    end
end
% Remove duplicate time points
[t_all_unique, ia] = unique(t_all, 'stable');
v_all_unique = v_all(:,ia);

% Interpolate to uniform time grid
%v_uniform = interp1(t_all_unique, v_all_unique', t_uniform)';  % → [num_nodes x num_points]
v_uniform = interp1(t_all_unique, v_all_unique', t_uniform, 'linear', 'extrap')';

% % Define uniform time points for analysis
% t_uniform = linspace(0, time_scale, num_points);
% v_uniform = interp1(t_all, v_all', t_uniform)'; % interpolate along columns
t = t_uniform;  % uniform time vector
v = v_uniform;  % interpolated voltage matrix

% 
% 
% % Define time points for results
% t = linspace(0, time_scale, num_points);
% [v, vd] = deval(sol, t);
%%
figure;

% Select nodes
stim_node = max_idx;
distant_node = randi([1 n]);

% Plot voltage at stimulated and distant node
plot(t_all, v_all(stim_node,:), 'LineWidth', 2); hold on;
plot(t_all, v_all(distant_node,:), '--', 'LineWidth', 2);

% Mark stimulation times
for pulse_t = stim_times
    xline(pulse_t,'k--','LineWidth',1.5,'Label',['Pulse ',num2str(pulse_t),' ms'],'LabelOrientation','horizontal');
end

xlabel('Time (ms)');
ylabel('Voltage (V)');
title('Corrected Voltage-Time Plot with rTMS Pulses');
legend(['Stimulated Node (',num2str(stim_node),')'], ...
       ['Distant Node (',num2str(distant_node),')']);
grid on;



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
% % Plot time-dependent metrics
% for i = 1:length(metrics_time_dependent)
%     metric_data = metrics_time_dependent{i};
%     for t_idx = 1:num_points
%         c.cdata = metric_data(:, t_idx);  % Map the data for this time step
%         filename_var = sprintf('%s_time_%d.dscalar.nii', metric_names_time_dependent{i}, t_idx);
%         ciftisavereset(c, filename_var);
%         fprintf('Saved %s at time step %d\n', metric_names_time_dependent{i}, t_idx);
%     end
% end
% % 
% Plot log-scaled maximum metrics
for i = 1:length(metrics_max)
    metric_data = metrics_max{i};



    log_metric_data = metric_data ;  % Apply log scale
    c.cdata = log_metric_data;  % Map the log-scaled data
    filename_var = sprintf('%s_log_scaled_%s.dscalar.nii', metric_names_max{i}, namestr);
    ciftisavereset(c, filename_var);
    fprintf('Saved log-scaled %s\n', metric_names_max{i});
end

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
load('filtered_selected_indices_high_res.mat');
%load('filtered_selected_indices.mat');
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


% Load or use existing:
% filtered_selected_indices: the node indices you want to remove from visualization

% Step 1: Mask filtered nodes
Imax(filtered_selected_indices) = 0;
Ipmax(filtered_selected_indices) = 0;
E1max(filtered_selected_indices) = 0;
Cpmax(filtered_selected_indices) = 0;
% Imax_full(filtered_selected_indices)  = NaN;
% Ipmax_full(filtered_selected_indices) = NaN;
% E1max_full(filtered_selected_indices) = NaN;
% Cpmax_full(filtered_selected_indices) = NaN;
Imax_full = zeros(num_original_nodes, 1);         % Maximum current
Ipmax_full = zeros(num_original_nodes, 1);        % Maximum positive current
E1max_full = zeros(num_original_nodes, 1);        % Maximum energy
Cpmax_full = zeros(num_original_nodes, 1);        % Maximum positive charge accumulated

% % Step 2: Create log-scaled versions
% log_Imax  = log10(Imax_full  + 1);
% log_Ipmax = log10(Ipmax_full + 1);
% log_Emax  = log10(E1max_full + 1);
% log_Cpmax = log10(Cpmax_full + 1);
%E1max_log = log10(E1max);
E1max_log = 1000000*E1max';
N_Down = 1:9656;
target_full = zeros(91282,1);
for group_id = 1:numel(group_indices)
    % Get the indices of nodes in this group
    nodes_in_group = group_indices{group_id};
    if group_id<=6192
        Imax_full(nodes_in_group) = Imax(group_id);
        Ipmax_full(nodes_in_group) = Ipmax(group_id);
        E1max_full(nodes_in_group) = E1max_log(group_id);
        Cpmax_full(nodes_in_group) = Cpmax(group_id);
    end
end

c = cifti_read('ones.dscalar.nii');



tstep = 1;
clabel0 = c.cdata;
clabel1  = zeros(size(clabel0,1),tstep);
for i = 1:size(clabel0,1)

  
        %clabel1(i)=SC(ReOrder(clabel0(i)),ReOrder(268));
        %clabel1(i)=D(ReOrder(clabel0(i)),ReOrder(268));
        clabel1(i)=Imax_full(i);
    
end


%%
% Save the new cifti file
c.cdata = clabel1;
filename_var = strcat('Imax_filtered_aud','.dscalar.nii');
ciftisavereset(c,filename_var);



c = cifti_read('ones.dscalar.nii');



tstep = 1;
clabel0 = c.cdata;
clabel1  = zeros(size(clabel0,1),tstep);
for i = 1:size(clabel0,1)

  
        %clabel1(i)=SC(ReOrder(clabel0(i)),ReOrder(268));
        %clabel1(i)=D(ReOrder(clabel0(i)),ReOrder(268));
        clabel1(i)=Ipmax_full(i);
    
end


%%
% Save the new cifti file
c.cdata = clabel1;
filename_var = strcat('Ipmax_filtered_aud','.dscalar.nii');
ciftisavereset(c,filename_var);


c = cifti_read('ones.dscalar.nii');



tstep = 1;
clabel0 = c.cdata;
clabel1  = zeros(size(clabel0,1),tstep);
for i = 1:size(clabel0,1)

  
        %clabel1(i)=SC(ReOrder(clabel0(i)),ReOrder(268));
        %clabel1(i)=D(ReOrder(clabel0(i)),ReOrder(268));
        clabel1(i)=Cpmax_full(i);
    
end


%%
% Save the new cifti file
c.cdata = clabel1;
filename_var = strcat('Cpmax_filtered_aud','.dscalar.nii');
ciftisavereset(c,filename_var);


c = cifti_read('ones.dscalar.nii');



tstep = 1;
clabel0 = c.cdata;
clabel1  = zeros(size(clabel0,1),tstep);
for i = 1:size(clabel0,1)

  
        %clabel1(i)=SC(ReOrder(clabel0(i)),ReOrder(268));
        %clabel1(i)=D(ReOrder(clabel0(i)),ReOrder(268));
        clabel1(i)=E1max_full(i);
    
end


%%
% Save the new cifti file
c.cdata = clabel1;
filename_var = strcat('Emax_filtered_aud','.dscalar.nii');
ciftisavereset(c,filename_var);



% 
% % Step 3: Load CIFTI template
% c = cifti_read('ones.dscalar.nii');  % Assume same length as 91k+ map
% clabel0 = c.cdata;
% 
% % Step 4: Save each as separate .dscalar.nii
% c.cdata = log_Imax;
% ciftisavereset(c, sprintf('%s_Imax_filtered.dscalar.nii', namestr));
% 
% c.cdata = log_Ipmax;
% ciftisavereset(c, sprintf('%s_Ipmax_filtered.dscalar.nii', namestr));
% 
% c.cdata = log_Emax;
% ciftisavereset(c, sprintf('%s_Emax_filtered.dscalar.nii', namestr));
% 
% c.cdata = log_Cpmax;
% ciftisavereset(c, sprintf('%s_Cpmax_filtered.dscalar.nii', namestr));
% 
% fprintf('Filtered dscalar files saved.\n');
% 
