%% Atlas Segmentation Model Simulation (Single Node Stimulation Version)

% Description: This script simulates an RC network based on a brain atlas, 
% using structural connectivity and distance data to model TMS-driven signal 
% propagation. It reads connectivity and distance matrices, processes them, 
% and solves an ODE for voltage propagation across brain regions. This code 
% version supports only single-node stimulation.

%% Adjustable Parameters

a = 1;                     % Subject ID or session ID for file naming
namestr = 'Atlas_Sch1000'; % Base name for output files
time_scale = 0.5;          % Total simulation time
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


%% Data Preparation

% Reorder stimulation location based on atlas order instead of reordering SC matrix
target_node = ReOrder(target_node);  % Adjust target node based on reordering

% Process structural connectivity matrix based on distance matrix
SCnew = SC(1:1000,1:1000) ./ D;
SCnew(isnan(SCnew)) = 0;
SCnew(SCnew == Inf) = 0;
SC = SCnew ./ max(SCnew(:)) * 1000;  % Normalize connectivity to 1000 scale

%% Simulation Setup

n = size(SC, 1);           % Number of nodes
C = ones(n, 1) * 1000;     % Capacitance for each node

% Initialize voltage vector with target node voltage for single-node stimulation
V = zeros(n, 1);
V(target_node) = initial_voltage;

% Run the ODE solver with capODESC_Sch function that takes SC and C
sol = ode23(@(t, v) capODESC_Sch(t, v, SC, C), [0 time_scale], V);

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
                current_val = (v(j, k) - v(i, k)) * 1000 * SC(j, i) / 2;
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

% Mapping function to match data points with Schaefer atlas regions using the reordered indices
function mapped_data = map_to_cifti(metric, clabel0, ReOrder, n)
    mapped_data = NaN(size(clabel0));
    for i = 1:numel(clabel0)
        if clabel0(i) > 0 && clabel0(i) <= n
            mapped_data(i) = metric(ReOrder(clabel0(i)));
        end
    end
end

% Plot and save results as dscalar files for Cpmax, Imax, Emax, and Ipmax
% with a name start with the name_str (in this case
% "Atlas_Sch1000_xx_316.dscalar.nii")
save_dscalar(map_to_cifti(log10(Cpmax + 1), clabel0, ReOrder, n), namestr, 'Cpmax');
save_dscalar(map_to_cifti(log10(Imax + 1), clabel0, ReOrder, n), namestr, 'Imax');
save_dscalar(map_to_cifti(log10(E1max + 1), clabel0, ReOrder, n), namestr, 'Emax');
save_dscalar(map_to_cifti(log10(Ipmax + 1), clabel0, ReOrder, n), namestr, 'Ipmax');

disp('Simulation and output generation complete.');


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
