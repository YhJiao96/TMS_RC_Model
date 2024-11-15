

%% High-Res Downsampled Model Simulation (Board E-field Stimulation Version)

% Description: This script simulates an RC network based on the high resolution 
% downsampled structural connectome, using downsampled connectivity (SC), 
% distance (D), and voltage (V) matrices to model TMS-driven signal propagation.

%% Adjustable Parameters

a = 1;                     % Subject ID or session ID for file naming
namestr = 'Downsampled_SC'; % Base name for output files
time_scale = 1;            % Total simulation time
num_points = 100;          % Number of points for time resolution

%% Load Input Files

% Load downsampled structural connectome and distance matrices
load('SC_downsampled.mat');  % Downsampled SC matrix
load('D_downsampled.mat');   % Downsampled distance matrix
load('V_downsampled.mat');   % Downsampled initial voltage vector

%% Data Preparation

% Define the resistance matrix from the downsampled SC matrix
R_downsampled = 1 ./ SC_downsampled; % Convert SC to resistance (R = 1 / SC)
R_downsampled(isnan(R_downsampled) | isinf(R_downsampled)) = 0; % Handle NaNs/Infs

% Normalize the SC matrix (optional, as done in Atlas model)
SC = SC_downsampled ./ max(SC_downsampled(:)) * 1000;

% Set the initial voltage from the downsampled V matrix
V = V_downsampled;

% Define capacitance matrix (uniform across nodes)
n = numel(V);  % Number of nodes in downsampled model
C = ones(n, 1) * 1000;

%% Simulation Setup

% Set up and run the ODE solver with the downsampled RC circuit model
sol = ode23(@(t, v) capODE_downsampled(t, v, SC, C), [0 time_scale], V);

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

% Calculate total current, power metrics, and max values over time
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

% Load the Schaefer atlas file to map results to brain regions (or a similar atlas file)
c = cifti_read('Schaefer2018_1000Parcels_7Networks_order.dlabel.nii');
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

% Plot and save results as dscalar files for Cpmax, Imax, Emax, and Ipmax
save_dscalar(map_to_cifti(log10(Cpmax + 1), clabel0, ReOrder, n), namestr, 'Cpmax');
save_dscalar(map_to_cifti(log10(Imax + 1), clabel0, ReOrder, n), namestr, 'Imax');
save_dscalar(map_to_cifti(log10(E1max + 1), clabel0, ReOrder, n), namestr, 'Emax');
save_dscalar(map_to_cifti(log10(Ipmax + 1), clabel0, ReOrder, n), namestr, 'Ipmax');

disp('Simulation and output generation complete.');

%% Function: capODE_downsampled with SC and C as Inputs

function dvdt = capODE_downsampled(t, v, SC, C)
    % Differential equations for RC circuit in downsampled model
    num_nodes = length(v);  
    dvdt = zeros(num_nodes, 1);  

    % Compute rate of voltage change at each node
    for i = 1:num_nodes
        for j = 1:num_nodes
            if i ~= j && SC(i, j) ~= 0
                dvdt(i) = dvdt(i) + (v(j) - v(i)) / (SC(i, j) * C(i));
            end
        end
    end
end
