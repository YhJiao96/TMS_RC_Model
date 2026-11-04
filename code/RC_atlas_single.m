%% Atlas Segmentation Model Simulation (Single Node Stimulation Version)

% Description: This script simulates an RC network based on a brain atlas, 
% using structural connectivity and distance data to model TMS-driven signal 
% propagation. It reads connectivity and distance matrices, processes them, 
% and solves an ODE for voltage propagation across brain regions. This code 
% version supports only single-node stimulation.

%% Adjustable Parameters

a = 1;                     % Subject ID or session ID for file naming
namestr = 'Atlas_Sch1000_dlpfc'; % Base name for output files
time_scale = 400;          % Total simulation time
num_points = 100;          % Number of points for time resolution
target_node = 268;         % Target node (TMS target for one-node stimulation)
initial_voltage = 1;       % Initial voltage at TMS target node (for one-node stimulation)

%% Load Input Files

% Load the Schaefer atlas reorder table to map node indices
SinaSchOrder = table2array(readtable('Schreorder.csv'));
ReOrder = SinaSchOrder(:, 1);  % Node reordering based on atlas

% Load connectivity matrix and distance matrix
load('SC_ds_hcp_tract_pro_ctx_sch7n1000p_sctx_fs14_n_1000_m_1000_group_mean_slcthr_5_nf_na.mat');  % SC matrix
load('ED_ctx_sch7n1000p_sctx_na.mat');  % Distance matrix
reverse_reorder = zeros(size(ReOrder));
reverse_reorder(ReOrder) = 1:numel(ReOrder);
reverse_reorder = reverse_reorder(1:998);
D1 = cat(2,D,zeros(1000,13));
Dn = cat(1,D1,zeros(13,1013));
Dror = Dn(reverse_reorder,:);
Dro = Dror(:,reverse_reorder);
D = Dro(1:998,1:998);

%% Data Preparation

% Reorder stimulation location based on atlas order instead of reordering SC matrix
target_node = ReOrder(target_node);  % Adjust target node based on reordering


% % Process structural connectivity matrix based on distance matrix
% SCnew = D(1:998,1:998) ;
% SCnew(isnan(SCnew)) = 0;
% SCnew(SCnew == Inf) = 0;
% SC = SCnew ./ max(SCnew(:)) ;  % Normalize connectivity to 1000 scale

%% Simulation Setup

n = size(SC, 1);           % Number of nodes
C = ones(n, 1) * 1000;     % Capacitance for each node
R = 1./SC;
R(isinf(a2aver))=0;
R(isnan(a2aver))=0;
% Initialize voltage vector with target node voltage for single-node stimulation
V = zeros(n, 1);
V(target_node) = initial_voltage;

% Run the ODE solver with capODESC_Sch function that takes SC and C
sol = ode23(@(t, v) capODESC_Sch(t, v, R, C), [0 time_scale], V);

% Define time points for results
t = linspace(0, time_scale, num_points);
[v, vd] = deval(sol, t);

%% Calculate Currents and Power for Each Node

Ip = zeros(n, num_points);
In = zeros(n, num_points);

for i = 1:n
    for j = 1:n
        if i ~= j && SC(j, i) ~= 0
            for k = 1:num_points
                current_val = (v(j, k) - v(i, k)) * 1000 ./R(j, i) ./ 2;
                if v(j, k) >= v(i, k)
                    Ip(i, k) = Ip(i, k) + current_val;
                else
                    In(i, k) = In(i, k) + current_val;
                end
            end
        end
    end
end

% Calculate total current, power-related metrics, and max values over time
I = In + Ip;
Cp = cumtrapz(abs(Ip), 2);
E1 = abs(I .* v);
E1c = cumtrapz(E1, 2);
Cpmax = max(Cp, [], 2);
E1max = max(E1c, [], 2);
Imax = max(I, [], 2);
Ipmax = max(Ip, [], 2);

%% Save Time-Series Data to dscalar Files

% Load the Schaefer atlas file to map results to brain regions
c = cifti_read('Schaefer2018_1000Parcels_7Networks_order.dlabel.nii');
clabel0 = c.cdata;

% Define the number of nodes and time points
n = size(SC, 1);
time_frames = num_points;

% Helper function to save each metric to a dscalar file
save_dscalar = @(data, label, suffix, frame) ...
    ciftisavereset(setfield(c, 'cdata', data), sprintf('%s_%s_t%d.dscalar.nii', label, suffix, frame));

% Mapping function to match data points with Schaefer atlas regions using the reordered indices
function mapped_data = map_to_cifti(metric, clabel0, ReOrder, n)
    mapped_data = NaN(size(clabel0));
    for i = 1:numel(clabel0)
        if clabel0(i) > 0 && clabel0(i) <= n
            mapped_data(i) = metric(ReOrder(clabel0(i)));
            
        end
    end
end

% Iterate over each time point to save voltage, current, positive current, and charge
for frame = 1:time_frames
    % Get the metrics for the current time frame
    V_frame = v(:, frame); 
    I_frame = I(:, frame);
    Ip_frame = Ip(:, frame);
    Cp_frame = Cp(:, frame);
    %c = cifti_read('Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR_Tian_Subcortex_S1.dlabel.nii');
    c = cifti_read('Schaefer2018_1000Parcels_7Networks_order.dlabel.nii'); % This is the Schaefer dlabel map we used for our atlas level Model.
    clabel0 = c.cdata;
    clabel1  = zeros(size(clabel0,1),1);
    for i = 1:size(clabel0,1)
        if clabel0(i)==0
            continue

        elseif (clabel0(i))==0
            continue
        elseif (clabel0(i))>size(V_frame,1)
            continue
        else

            clabel1(i) = Ip_frame(ReOrder(clabel0(i)));

        end
    end


    c.cdata = clabel1;
    filename_var = strcat('Atlas_Ip_',num2str(frame),'.dscalar.nii');
    ciftisavereset(c,filename_var);



    % Map each metric to the dscalar format for this time frame
    V_dscalar = map_to_cifti(V_frame, clabel0, ReOrder, n);
    I_dscalar = map_to_cifti(I_frame, clabel0, ReOrder, n);
    Ip_dscalar = map_to_cifti(Ip_frame, clabel0, ReOrder, n);
    Cp_dscalar = map_to_cifti(Cp_frame, clabel0, ReOrder, n);

    % Save each metric as a dscalar file
    save_dscalar(V_dscalar, namestr, 'Voltage', frame);
    save_dscalar(I_dscalar, namestr, 'Current', frame);
    save_dscalar(Ip_dscalar, namestr, 'PosCurrent', frame);
    save_dscalar(Cp_dscalar, namestr, 'Charge', frame);
end

disp('Time-series data saved as dscalar files.');



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

%% Plotting to dscalar Files for Connectome Workbench

% Load the Schaefer atlas file to map results to brain regions
c = cifti_read('Schaefer2018_1000Parcels_7Networks_order.dlabel.nii');
clabel0 = c.cdata;

% Helper function to save each metric to a dscalar file
save_dscalar = @(data, label, suffix) ...
    ciftisavereset(setfield(c, 'cdata', data),  sprintf('%s_%s_%d.dscalar.nii', label, suffix, target_node));

% % Mapping function to match data points with Schaefer atlas regions using the reordered indices
% function mapped_data = map_to_cifti(metric, clabel0, ReOrder, n)
%     mapped_data = NaN(size(clabel0));
%     for i = 1:numel(clabel0)
%         if clabel0(i) > 0 && clabel0(i) <= n
%             mapped_data(i) = metric(ReOrder(clabel0(i)));
%         end
%     end
% end

% Plot and save results as dscalar files for Cpmax, Imax, Emax, and Ipmax
% with a name start with the name_str (in this case
% "Atlas_Sch1000_xx_316.dscalar.nii")
save_dscalar(map_to_cifti(log10(Cpmax + 1), clabel0, ReOrder, n), namestr, 'Cpmax');
save_dscalar(map_to_cifti(log10(Imax + 1), clabel0, ReOrder, n), namestr, 'Imax');
save_dscalar(map_to_cifti(log10(E1max + 1), clabel0, ReOrder, n), namestr, 'Emax');
save_dscalar(map_to_cifti(log10(Ipmax + 1), clabel0, ReOrder, n), namestr, 'Ipmax');

disp('Simulation and output generation complete.');

%%
noise_loc = [];
%noise_loc  = [96 101 92 108 93 109 90 99 238 234 253 256 106 95 89 97 111 98 105 86 259 263 102 265 103 87 83 607 601 597 606 754 613 612 600 603 611 589 598 596 590 608 787 776 756 614 609 610 ];%already reordered auditory regions
c2 = cifti_read('FMRIsurface.dscalar.nii');
Surf1  = (c2.cdata);
a2 = zeros(1012,1);
%vertices = zeros(1012,3);
a2num = zeros(1012,1);
for i = 1:1012
    for j = 1:size(c1.cdata,1)
        if c1.cdata(j) == i
            a2(i) = a2(i)+abs(Surf1(j)); %abs optional
            a2num(i) = a2num(i)+1;
            
            %Atlas_normE(i)=Atlas_normE(i)+normE(j);
        else 
            continue
        end
    end
end
a2aver = a2./a2num;
a2aver(isinf(a2aver))=0;
a2aver(isnan(a2aver))=0;
a2n = a2aver(1:1000);
a2n_ro = zeros(max(ReOrder),1);
for i = 1:size(a2aver)
    if i>=1000
        continue
    elseif ReOrder(i) == 0 
        continue
    else
        a2n_ro(ReOrder(i)) = a2aver(i);
    end
end
a2n_ro(ReOrder(noise_loc)) = 0;
a2 = a2n_ro(1:998);
% max_half = max(find(ReOrder<=1000));
% a2n = a2n(1:max_half);
%a2n(noise_loc) = 0;

%%

metric_names = {'Imax', 'Ipmax', 'Cpmax', 'E1max'};
metrics_storage = zeros(1, length(metric_names));  % Store correlation values
p_values_storage = zeros(1, length(metric_names));  % Store p-values
metrics = {Imax, Ipmax, Cpmax, E1max};
for j = 1:length(metrics)
    a1t = zeros(size(a2));
    a1t = metrics{j};
    [s1,p1] = corr(a1t,a2,'Type', 'Spearman');
    fprintf('Correlation results for correlation and p-value\n');
    fprintf('  %s Correlation: %.4f, p-value: %.4f\n', metric_names{j}, s1, p1);
end

for j = 1:length(metrics)
    a1n = zeros(size(a2));  % Preallocate

    a1_metric = metrics{j};  % Current metric data
    %a1n(valid_indices) = a1_metric(ReOrder(valid_indices));  % Reorder the data based on ReOrder
    %a1n((noise_loc)) = 0;  % Set noise locations to 0
    a1n = a1_metric;
    %a1n(valid_indices) = a1_metric((valid_indices));
    %a1n(ReOrder(noise_loc)) = 0;  % Set noise locations to 0
    % Filter out invalid values in a2 and the current metric data (a1n)
    valid_mask = ~(a1n == 0 | isnan(a1n) | isinf(a1n) | a2 == 0 | isnan(a2) | isinf(a2));
    a1n_filtered = a1n(valid_mask);  % Filtered metric data
    a2_filtered = a2(valid_mask);    % Filtered fMRI data

    % Calculate Spearman correlation and p-value between filtered a1n and a2
    [spearman_corr, p_value] = corr(a1n_filtered, a2_filtered, 'Type', 'Spearman');

    % Store results for each power and metric
    metrics_storage(1, j) = spearman_corr;
    p_values_storage(1, j) = p_value;
end

fprintf('Correlation results for correlation and p-value\n');

    
for j = 1:length(metric_names)
    fprintf('  %s Correlation: %.4f, p-value: %.4f\n', metric_names{j}, metrics_storage(1, j), p_values_storage(1, j));
end

%% Function: capODESC_Sch with R and C as Inputs

function dvdt = capODESC_Sch(t, v, R, C)
    % Differential equations for RC circuit in atlas model
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



time = 1:100;  % Define the time vector
figure;
plot(time, v');  % Plot each row as a curve (transpose to align rows as curves)
title('Node Status Over Time');
xlabel('Time Frame');
ylabel('Node Status');
grid on;

%% Compute Correlations for Each Time Frame
time_frames = 20;  % Total number of time frames
correlations = zeros(4, time_frames);  % Store correlations for each metric
metric_names = {'I', 'Ip', 'E', 'Cp'};  % Metric names for plotting
metrics = {I, Ip, E1, Cp};  % Corresponding metric matrices

% Loop over time frames and compute correlations
for frame = 1:time_frames
    for metric_idx = 1:length(metrics)
        metric_data = metrics{metric_idx}(:, frame);  % Data for the current metric at this time frame

        % Filter out invalid values in metric_data and a2
        valid_mask = ~(metric_data == 0 | isnan(metric_data) | isinf(metric_data) | ...
                       a2 == 0 | isnan(a2) | isinf(a2));
        metric_filtered = metric_data(valid_mask);
        a2_filtered = a2(valid_mask);

        % Compute Spearman correlation
        if ~isempty(metric_filtered) && ~isempty(a2_filtered)
            correlations(metric_idx, frame) = corr(metric_filtered, a2_filtered, 'Type', 'Spearman');
        else
            correlations(metric_idx, frame) = NaN;  % Assign NaN if no valid data exists
        end
    end
end

%% Plot Correlation Change Over Time
figure;
hold on;
time = 1:time_frames;  % Time vector

% Plot correlation curves for all metrics
for metric_idx = 1:length(metrics)
    plot(time, correlations(metric_idx, :), '-o', 'LineWidth', 2, ...
         'DisplayName', [metric_names{metric_idx} ' Correlation']);
end

% Customize the plot
title('Correlation of Metrics with a2 Over Time Frames');
xlabel('Time Frame');
ylabel('Spearman Correlation');
legend('Location', 'best');
grid on;
hold off;

%% Save the Plot
saveas(gcf, 'Correlation_Over_Time.png');
