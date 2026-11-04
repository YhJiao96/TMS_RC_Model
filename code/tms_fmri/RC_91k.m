

%% High-Res Model Simulation (Board E-field Stimulation Version)

% Description: This script simulates an RC network based on the high resolution 
% structural connectome, using full high resolution group level connectivity (SC), 
% distance (D), and voltage (V) matrices (node size as 91282) to model 
% TMS-driven signal propagation.
% Warning: This code needs a memeory of computer at least 64 GB.

%% Adjustable Parameters

a = 1;                     % Subject ID or session ID for file naming
namestr = 'High_Res_SC'; % Base name for output files
time_scale = 100;            % Total simulation time
num_points = 100;          % Number of points for time resolution

%% Load Input Files

% Load downsampled structural connectome and distance matrices
load /data/gpfs/projects/punim2041/SCCount.mat; % load the SC matrix 
GroupSC_Count = spconvert(SCCount);
DiagGSC = diag(diag(GroupSC_Count));
SC = GroupSC_Count-DiagGSC;

load /data/gpfs/projects/punim2041/SCLength.mat; % load the Distance matrix 
GroupSC_Length = spconvert(SCLength);
DiagGSL = diag(diag(GroupSC_Length));
D = GroupSC_Length - DiagGSL;

filename = strcat('/data/gpfs/projects/punim2041/subject_coord_91k_normE.csv');
E = importdata(filename);
V = E;

% Calculate the 99th percentile of the voltage matrix V across all time steps
threshold_V = prctile(V(:), 99);  % 90th percentile of voltage values

% Apply the threshold: Set values below the threshold to zero
V(V < threshold_V) = 0;  % Replace values below the threshold


%% Data Preparation

% Define the resistance matrix from the downsampled SC matrix
R = D ./ SC; % Convert SC to resistance (R = 1 / SC)
R(isnan(R) | isinf(R)) = 0; % Handle NaNs/Infs



% Define capacitance matrix (uniform across nodes)
n = numel(V);  % Number of nodes in downsampled model
C = ones(n, 1) * 1000;

%% Simulation Setup

% Set up and run the ODE solver with the downsampled RC circuit model
sol = ode23(@(t, v) capODE_high_Res(t, v, R, C), [0 time_scale], V);

% Define time points for results
t = linspace(0, time_scale, num_points);
[v, vd] = deval(sol, t);

%% Calculate Currents and Power for Each Node

Ip = zeros(n, num_points);
In = zeros(n, num_points);

for i = 1:n
    for j = 1:n
        if i ~= j && R(j, i) ~= 0
            for k = 1:num_points
                current_val = (v(j, k) - v(i, k)) * 1000 ./ R(j, i) ./ 2;
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


%% Function: capODE_downsampled with SC and C as Inputs

function dvdt = capODE_high_Res(t, v, R, C)
    % Differential equations for RC circuit in downsampled model
    num_nodes = length(v);  
    dvdt = zeros(num_nodes, 1);  

    % Compute rate of voltage change at each node
    for i = 1:num_nodes
        for j = 1:num_nodes
            if i ~= j && R(i, j) ~= 0
                dvdt(i) = dvdt(i) + (v(j) - v(i)) ./ (R(i, j) .* C(i));
            end
        end
    end
end
