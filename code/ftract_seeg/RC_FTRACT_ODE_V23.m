load('SC_data_neuron_2023_yihang.mat');
load('idx_sc.mat');
load('ReOreder_F_tract.mat');
load('ReOrder_Caio_SC.mat');
load('SC_data_neuron_2023_yihang.mat');
%% There is a ranking mismatch in the code, need to be checked!

atlasDir = '/Users/yihjiao/Documents/MATLAB/F-tract-2';
%age_group = '0_15';
parcellation = 'MNI-HCP-MMP1';
%parcellation = 'Lausanne2008-500';

N_file = fullfile(atlasDir,parcellation,'N_implantations.txt.gz');
output_dir = fullfile(atlasDir,parcellation);  % You can also use atlasDir if you prefer
gunzip(N_file, output_dir);

txt_file = fullfile(output_dir, 'N_implantations.txt');
N_data_matrix = readtable(txt_file); 
N_data_matrix_num = N_data_matrix(8:end,1); 

% Convert table (360x1) of strings to cell array of strings
N_row_strings = table2cell(N_data_matrix_num);  % now 360x1 cell

% Parse each row string into a numeric row
row_values = cellfun(@(s) sscanf(s, '%f')', N_row_strings, 'UniformOutput', false);

% Combine into final matrix
N_final_matrix = vertcat(row_values{:});  % should now be 360 x 360

%%

%amp_file = fullfile(atlasDir,parcellation,'amplitude__mad.txt.gz');
amp_file = fullfile(atlasDir,parcellation,'amplitude__median.txt.gz');
output_dir = fullfile(atlasDir,parcellation);  % You can also use atlasDir if you prefer
gunzip(amp_file, output_dir);

%txt_file = fullfile(output_dir, 'amplitude__mad.txt');
txt_file = fullfile(output_dir, 'amplitude__median.txt');
amp_data_matrix = readtable(txt_file); 
amp_data_matrix_num = amp_data_matrix(8:end,1); 

% Convert table (360x1) of strings to cell array of strings
amp_row_strings = table2cell(amp_data_matrix_num);  % now 360x1 cell

% Parse each row string into a numeric row
row_values = cellfun(@(s) sscanf(s, '%f')', amp_row_strings, 'UniformOutput', false);

% Combine into final matrix
amp_final_matrix = vertcat(row_values{:});  % should now be 360 x 360



prob_file = fullfile(atlasDir,parcellation,'probability.txt.gz');
output_dir = fullfile(atlasDir,parcellation);  % You can also use atlasDir if you prefer
gunzip(prob_file, output_dir);

txt_file = fullfile(output_dir, 'probability.txt');
prob_data_matrix = readtable(txt_file); 
prob_data_matrix_num = prob_data_matrix(8:end,1); 

% Convert table (360x1) of strings to cell array of strings
prob_row_strings = table2cell(prob_data_matrix_num);  % now 360x1 cell

% Parse each row string into a numeric row
row_values = cellfun(@(s) sscanf(s, '%f')', prob_row_strings, 'UniformOutput', false);

% Combine into final matrix
prob_final_matrix = vertcat(row_values{:});  % should now be 360 x 360


N_file = fullfile(atlasDir,parcellation,'N_stimulations.txt.gz');
output_dir = fullfile(atlasDir,parcellation);  % You can also use atlasDir if you prefer
gunzip(N_file, output_dir);

txt_file = fullfile(output_dir, 'N_stimulations.txt');
N_data_matrix = readtable(txt_file); 
N_data_matrix_num = N_data_matrix(8:end,1); 

% Convert table (360x1) of strings to cell array of strings
N_row_strings = table2cell(N_data_matrix_num);  % now 360x1 cell

% Parse each row string into a numeric row
row_values = cellfun(@(s) sscanf(s, '%f')', N_row_strings, 'UniformOutput', false);

% Combine into final matrix
N_sti_final_matrix = vertcat(row_values{:});  % should now be 360 x 360

%amp_file = fullfile(atlasDir,parcellation,'amplitude__mad.txt.gz');
peak_delay_file = fullfile(atlasDir,parcellation,'amplitude__median.txt.gz');
output_dir = fullfile(atlasDir,parcellation);  % You can also use atlasDir if you prefer
gunzip(peak_delay_file, output_dir);

%txt_file = fullfile(output_dir, 'amplitude__mad.txt');
txt_file = fullfile(output_dir, 'amplitude__median.txt');
peak_delay_data_matrix = readtable(txt_file); 
peak_delay_data_matrix_num = peak_delay_data_matrix(8:end,1); 

% Convert table (360x1) of strings to cell array of strings
peak_delay_row_strings = table2cell(peak_delay_data_matrix_num);  % now 360x1 cell

% Parse each row string into a numeric row
row_values = cellfun(@(s) sscanf(s, '%f')', peak_delay_row_strings, 'UniformOutput', false);

% Combine into final matrix
peak_delay_final_matrix = vertcat(row_values{:});  % should now be 360 x 360


% Extract the unified order
%[~, sort_idx_col1] = sort(ReOrder_F_tract(:,1));  % for Amp rows
%[~, sort_idx_col2] = sort(ReOrder_F_tract(:,2));  % for SC columns
sc_order_1 = ReOrder_Caio_SC; 
d_order_1 = ReOrder_Caio_SC;
f_order_1 = ReOrder_Caio_SC;
amp_order_1  = ReOrder_F_tract(:, 2);

% Get the matching order from original indices
% sc_order = ReOrder_F_tract(sort_idx_col1, 1);  % should be 1:360
% amp_order  = ReOrder_F_tract(sort_idx_col2, 2);

% Reorder Amp (row-wise by column 1)
%SC_sorted = wb_sc_avg(sort_idx_col1, sort_idx_col1);  % symmetric
Amp_sorted = amp_final_matrix(idx_sc(1:360), idx_sc(1:360));
Prob_sorted = prob_final_matrix(idx_sc(1:360), idx_sc(1:360));
Amp_sorted = amp_final_matrix(amp_order_1, amp_order_1);
Prob_sorted = prob_final_matrix(amp_order_1, amp_order_1);
%Edist = euclidian_distance_value(sort_idx_col1, sort_idx_col1);
N_imp_prob = N_final_matrix(amp_order_1, amp_order_1);
N_sti = N_sti_final_matrix(amp_order_1, amp_order_1);
N_imp = N_sti.*Prob_sorted;

amplitude_value_cor = Amp_sorted;
prob_value_cor = Prob_sorted;

SC_sorted = SC(sc_order_1, sc_order_1);
D_sorted = D(d_order_1,d_order_1);
F_sorted = F(f_order_1,f_order_1);
SC = SC_sorted; % strength
%SC(F_sorted<0.25)=0;

SC(isnan(SC) | isinf(SC)) = 0;
D_sorted(isnan(D_sorted) | isinf(D_sorted)) = 0;

% Avoid division by zero
%R = 1./(D.*SC);
R =  D_sorted./ SC;
R(isnan(R) | isinf(R)) = 0;

n = size(R, 1);
C = 1 * ones(n, 1);  % uniform capacitance

time_scale = 200;     % ms
num_points = 5000;
dt = time_scale / num_points;
t = linspace(0, time_scale, num_points);

function dVdt = rc_ode_rhs(t, V, R, C, I_input)
    n = length(C);
    dVdt = zeros(n, 1);
    for i = 1:n
        inflow = 0;
        for j = 1:n
            if R(i,j) > 0
                inflow = inflow + (V(j) - V(i)) / R(i,j);
            end
        end
        dVdt(i) = (inflow + I_input(i)) / C(i);
    end
end
function I_matrix = compute_current_matrix(R, V_out, I_input)
    n = size(V_out,1);
    T = size(V_out,2);
    I_matrix = zeros(n, T);
    for t = 1:T
        V = V_out(:,t);
        for i = 1:n
            inflow = 0;
            for j = 1:n
                if R(i,j) > 0
                    inflow = inflow + (V(j) - V(i)) / R(i,j);
                end
            end
            I_matrix(i,t) = inflow + I_input(i);
        end
    end
end

t_span = linspace(0, time_scale, num_points);
V_results = zeros(n, num_points, n);  % [target_node, time, stim_node]
I_max = zeros(n, n, n);  
% --- Loop over each node as the stim node ---
for stim_idx = 1:n
    I_input = zeros(n, 1);
    %I_input(stim_idx) = 1;  % Stimulate this node with unit input

    V0 = zeros(n, 1);  % Initial voltage is zero everywhere
    V0(stim_idx) = 1;
    % Solve the ODE using ode45
    ode_rhs = @(t, V) rc_ode_rhs(t, V, R, C, I_input);
    [~, V_traj] = ode45(ode_rhs, t_span, V0);
    V_traj = V_traj';  % Now V_traj is [n x num_points]
    % Store result: transpose to match shape [n x num_points]
    V_results(:,:,stim_idx) = V_traj;
    I_temp_max = zeros(n, n);  % Per stimulation max current
    for t_idx = 1:num_points
        V_t = V_traj(:, t_idx);
        I_mat = zeros(n, n);
        for i = 1:n
            for j = 1:n
                if R(i,j) > 0
                    I_mat(i,j) = (V_t(i) - V_t(j)) / R(i,j);
                end
            end
        end
        I_temp_max = max(I_temp_max, abs(I_mat));
    end
    I_max(:,:,stim_idx) = I_temp_max;
    fprintf('Stim node %d done\n', stim_idx);
end

I_node_max = zeros(n, n);  % [node, stim]

for stim = 1:n
    V_traj = V_results(:,:,stim);  % [n x time]
    node_current = zeros(n, 1);

    for t = 1:num_points
        V = V_traj(:, t);
        total_I = zeros(n, 1);
        for i = 1:n
            for j = 1:n
                if R(i,j) > 0
                    I_ij = (V(i) - V(j)) / R(i,j);
                    total_I(i) = total_I(i) + abs(I_ij);  % outflow from node i
                end
            end
        end
        node_current = max(node_current, total_I);  % max across time
    end
    I_node_max(:, stim) = node_current;
end

% 
% %% --- Precompute conductance graph and Laplacian ---
% SC(isnan(SC) | isinf(SC)) = 0;
% D(isnan(D)   | isinf(D))  = 0;
% % Avoid division by zero
% R =  D_sorted./ SC;
% R(isnan(R) | isinf(R)) = 0;
% 
% n = size(R, 1);
% C = 1 * ones(n, 1);  % uniform capacitance
% 
% G = zeros(n,n);
% maskG = (SC > 0) & (D > 0);
% G(maskG) = SC(maskG) ./ D(maskG);   % = 1./R
% 
% % Optional: sparsify for speed if graph is large
% % G = sparse(G);
% 
% deg = sum(G,2);
% L = diag(deg) - G;                  % Graph Laplacian (n×n)
% 
% % Column vector capacitances, safeguard zeros
% C = C(:); 
% C(C==0) = eps;
% Cinv = 1 ./ C;
% 
% t_span = linspace(0, time_scale, num_points);
% 
% % Outputs [node × stim]
% I_node_max       = zeros(n, n);  % (you already had this)
% I_pos_node_max   = zeros(n, n);  % NEW: max positive outgoing current
% E_cap_node_accum = zeros(n, n);  % NEW: accumulated capacitor energy ∫(1/2 C V^2) dt
% Q_node_max       = zeros(n, n);  % NEW: max charge |C V|
% 
% % Helper for ODE (vectorized)
% ode_rhs = @(t, V) ( -L*V + zeros(n,1) ) .* Cinv;  % I_input = 0 (since you excite via V0)
% 
% for stim_idx = 1:n
%     % Initial condition: unit voltage at stim node (impulse-like)
%     V0 = zeros(n,1);
%     V0(stim_idx) = 1;
% 
%     % Solve
%     [~, V_traj] = ode45(ode_rhs, t_span, V0);
%     V_traj = V_traj.';                % [n × T]
%     T = size(V_traj,2);
% 
%     % ---------- Currents (vectorized per frame) ----------
%     % I_ij(t) = G_ij * (V_i - V_j). We'll accumulate per-node totals.
%     node_I_max      = zeros(n,1);     % total |I| out of node i (max over t)
%     node_Ipos_max   = zeros(n,1);     % total positive outflow only (max over t)
% 
%     % ---------- Capacitor energy & charge ----------
%     % E_i(t) = 1/2 * C_i * V_i(t)^2, Q_i(t) = C_i * V_i(t)
%     E_i_t = 0.5 .* (C .* (V_traj.^2));             % [n × T]
%     E_cap_node_accum(:, stim_idx) = trapz(t_span, E_i_t, 2);  % ∫ E_i(t) dt per node
%     Q_node_max(:, stim_idx)       = max(abs(C .* V_traj), [], 2);
% 
%     % Sweep over time to get current-based maxima
%     for t_idx = 1:T
%         v = V_traj(:, t_idx);                      % [n×1]
% 
%         % Outer difference ΔV_ij = V_i - V_j, then multiply by conductance
%         % NOTE: for dense G this is O(n^2) per frame; with sparse G it’s fast.
%         dV = v - v.';                              % [n×n]
%         I_mat = G .* dV;                           % [n×n], current from i->j
% 
%         % Per-node outgoing totals
%         I_out_abs = sum(abs(I_mat), 2);            % total magnitude out of node i
%         I_out_pos = sum(max(I_mat, 0), 2);         % only positive (i->j where V_i>V_j)
% 
%         node_I_max    = max(node_I_max,    I_out_abs);
%         node_Ipos_max = max(node_Ipos_max, I_out_pos);
%     end
% 
%     I_node_max(:,     stim_idx) = node_I_max;
%     I_pos_node_max(:, stim_idx) = node_Ipos_max;
% 
%     if mod(stim_idx,25)==0 || stim_idx==n
%         fprintf('Stim node %d/%d done\n', stim_idx, n);
%     end
% end
% 
%%
% Graph Theory
MFPT = mean_first_passage_time(SC);
diff_eff = 1./MFPT;
diff_eff(isinf(diff_eff))=0;

SI = search_information(SC, 1./SC);
SI(isinf(SI))=0;

s = diag(sum(SC,1));
CMY = expm(s^(-1/2)*SC*s^(-1/2));

L = 1./SC;
L(L == Inf) = 0;
[sr, PL_bin, PL_wei, PL_dis, paths] = navigation_wu(L, D_sorted);
Nav_eff = 1./PL_wei;
Nav_eff(isinf(Nav_eff))=0;

[SPL,hops,Pmat] = distance_wei_floyd(1./SC);
%[SPL,hops,Pmat] = distance_wei_floyd(-log(SC));
SPE = 1./SPL;
SPE(isinf(SPL))=0;

