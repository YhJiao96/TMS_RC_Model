load('SC_data_neuron_2023_yihang.mat');
load('idx_sc.mat');
load("Y_struct_v1-v17_PAL_res.mat");
load("X_struct_v1-19_FCO.mat");
amp_final_matrix = Y_struct.v5.A;
prob_final_matrix = Y_struct.v5.P;
emp_idx = Y_struct.v5.parcel_list_sidx;
%%

idx_sc = 1:374;
sc_order_1 = idx_sc; 
d_order_1 = idx_sc;
f_order_1 = idx_sc;
SC =  X_struct.v8.W;
D = X_struct.v8.ED;
F = X_struct.v8.FCO;
Amp_sorted = amp_final_matrix;
Prob_sorted = prob_final_matrix;

N_prob = Y_struct.v5.Ni;
N_sti = Y_struct.v5.Ns;
N_imp = Y_struct.v5.Ni;
N_resp = Y_struct.v5.Nr; 


amplitude_value_cor = Amp_sorted;
prob_value_cor = Prob_sorted;

SC_sorted = SC(sc_order_1, sc_order_1);
D_sorted = D(d_order_1,d_order_1);
F_sorted = F(f_order_1,f_order_1);
SC = SC_sorted; % strength
D = D_sorted;
%SC(F_sorted<0.25)=0;

%

SC(isnan(SC) | isinf(SC)) = 0;
D_sorted(isnan(D_sorted) | isinf(D_sorted)) = 0;

% Avoid division by zero
R =  D_sorted./ SC;
R(isnan(R) | isinf(R)) = 0;

n = size(R, 1);
C = 0.1 * ones(n, 1);  % uniform capacitance
ratio = 5;
%C(361:end) = ratio;
time_scale = 200;     % ms
num_points = 2000;
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
I_node_t = zeros(n,n,num_points);
E_node_t      = zeros(n, num_points, n);  % [node, time, stim]
E_node_accum  = zeros(n, n);              % time-integrated energy per node & stim
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
    total_I_t = zeros(n,num_points);
    I_temp_max = zeros(n, n);  % Per stimulation max current
        % ---- Capacitor energy per node & time ----
    % V_traj: [n x num_points], C: [n x 1]
    % E_i(t) = 0.5 * C_i * V_i(t)^2
    E_t = 0.5 * (C .* (V_traj.^2));   % [n x num_points]
    
    % Store full time course
    E_node_t(:,:,stim_idx) = E_t;
    
    % Time-integrated (accumulated) energy using trapezoidal rule
    E_node_accum(:, stim_idx) = trapz(t_span, E_t, 2);   % integrate along time (dim 2)

    for t_idx = 1:num_points
        V_t = V_traj(:, t_idx);
        I_mat = zeros(n, n);
        for i = 1:n
            for j = 1:n
                if R(i,j) > 0
                    I_mat(i,j) = (V_t(i) - V_t(j)) / R(i,j);
                    total_I_t(i,t_idx) = total_I_t(i,t_idx)+I_mat(i,j);
                end
            end
        end
        
        I_temp_max = max(I_temp_max, abs(I_mat));
    end
    I_max(:,:,stim_idx) = I_temp_max;
    I_node_t(stim_idx,:,:) =  total_I_t;
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

%%
% Correlation for Prob
N_thres = N_sti;
emp_vec = Prob_sorted;
T = 100;

metrics = {
    'I_node_max', I_node_max(1:360,1:360)
    'Diff_eff',   diff_eff(1:360,1:360)
    'SI',         SI(1:360,1:360)
    %'CMY',        CMY(1:360,1:360)
    'Nav_eff',    Nav_eff(1:360,1:360)
    'SPE',        SPE(1:360,1:360)
    'ED',         D_sorted(1:360,1:360)      % Euclidean distance
    'SC',         SC(1:360,1:360)     % structural connectivity
    'Energy_acc'  E_node_accum(1:360,1:360)
};

nMetrics = size(metrics,1);
names   = cell(nMetrics,1);
r_vals  = nan(nMetrics,1);
p_vals  = nan(nMetrics,1);
pr_vals = nan(nMetrics,1);
pp_vals = nan(nMetrics,1);
count_vals = nan(nMetrics,1);

D_cor = D_sorted(1:360,1:360);
F_cor = F_sorted(1:360,1:360);

for m = 1:size(metrics,1)
    name = metrics{m,1};
    sim_vec = metrics{m,2};
    mask = (N_thres > T);
    %mask = triu(true(360),1);
    mask = mask & (blkdiag(ones(180), ones(180)));
    %mask = mask & (N_thres > T);
    mask = mask & (D_cor > 20);
    %mask = mask & (F_cor < 0.25);
    %mask = mask & (F_cor > 0.75);
    mask = mask & isfinite(emp_vec) ;
    %mask = mask & isfinite(emp_vec) & isfinite(sim_vec);

    emp_vals = emp_vec(mask);
    sim_vals = sim_vec(mask);
    D_vals = D_cor(mask);

    if isempty(emp_vals)
        error('Mask selected zero pairs. Check T and D thresholds.');
    end

    [r_spear, p_spear]   = corr(emp_vals, sim_vals, 'Type','Spearman');

    [pr_spear, pp_spear] = partialcorr(emp_vals, sim_vals, D_vals, ...
        'Type','Spearman', 'Rows','complete');


    names{m}  = name;
    r_vals(m) = r_spear;
    p_vals(m) = p_spear;
    pr_vals(m)= pr_spear;
    pp_vals(m)= pp_spear;
    count_vals(m) = nnz(mask);
end

results_table = table(names, count_vals, r_vals, p_vals, pr_vals, pp_vals, ...
    'VariableNames', {'Metric','Counts','Spearman_r','Spearman_p','Partial_r','Partial_p'});

disp(results_table);

results_matrix = [r_vals p_vals pr_vals pp_vals];


% Compute absolute values
abs_spear   = abs(results_table.Spearman_r);
abs_partial = abs(results_table.Partial_r);

% Sort by Spearman absolute correlation
[~, sort_idx] = sort(abs_spear, 'descend');

% Reorder everything
X = categorical(results_table.Metric(sort_idx));
X = reordercats(X, results_table.Metric(sort_idx));
Y = [abs_spear(sort_idx), abs_partial(sort_idx)];

% Plot grouped bars
figure; hold on;
bar(X, Y, 'grouped');

% Labels and title
legend({'|Spearman|','|Partial|'}, 'Location','best');
ylabel('Absolute Correlation (|r|)');
xlabel('Metric');
title(sprintf('Absolute Correlation with Empirical Data:Prob (T > %g, D > 20)', T));

grid on;



%%
% Correlation for Amp
N_thres = N_resp;
emp_vec = Amp_sorted;
T = 100;

metrics = {
    'I_node_max', I_node_max(1:360,1:360)
    'Diff_eff',   diff_eff(1:360,1:360)
    'SI',         SI(1:360,1:360)
    %'CMY',        CMY(1:360,1:360)
    'Nav_eff',    Nav_eff(1:360,1:360)
    'SPE',        SPE(1:360,1:360)
    'ED',         D_sorted(1:360,1:360)      % Euclidean distance
    'SC',         SC(1:360,1:360)     % structural connectivity
    'Energy_acc'  E_node_accum(1:360,1:360)
};

nMetrics = size(metrics,1);
names   = cell(nMetrics,1);
r_vals  = nan(nMetrics,1);
p_vals  = nan(nMetrics,1);
pr_vals = nan(nMetrics,1);
pp_vals = nan(nMetrics,1);
count_vals = nan(nMetrics,1);

D_cor = D_sorted(1:360,1:360);
F_cor = F_sorted(1:360,1:360);

for m = 1:size(metrics,1)
    name = metrics{m,1};
    sim_vec = metrics{m,2};
    mask = (N_thres > T);
    %mask = mask & triu(true(360),1);
    %mask = mask & (blkdiag(ones(180), ones(180)));
    %mask = mask & (N_thres > T);
    mask = mask & (D_cor > 20);
    mask = mask & (F_cor < 0.25);
    %mask = mask & (F_cor > 0.75);
    mask = mask & isfinite(emp_vec) ;
    %mask = mask & isfinite(emp_vec) & isfinite(sim_vec);

    emp_vals = emp_vec(mask);
    sim_vals = sim_vec(mask);
    D_vals = D_cor(mask);

    if isempty(emp_vals)
        error('Mask selected zero pairs. Check T and D thresholds.');
    end

    [r_spear, p_spear]   = corr(emp_vals, sim_vals, 'Type','Spearman');

    [pr_spear, pp_spear] = partialcorr(emp_vals, sim_vals, D_vals, ...
        'Type','Spearman', 'Rows','complete');


    names{m}  = name;
    r_vals(m) = r_spear;
    p_vals(m) = p_spear;
    pr_vals(m)= pr_spear;
    pp_vals(m)= pp_spear;
    count_vals(m) = nnz(mask);
end

results_table = table(names, count_vals, r_vals, p_vals, pr_vals, pp_vals, ...
    'VariableNames', {'Metric','Counts','Spearman_r','Spearman_p','Partial_r','Partial_p'});

disp(results_table);

results_matrix = [r_vals p_vals pr_vals pp_vals];


% Compute absolute values
abs_spear   = abs(results_table.Spearman_r);
abs_partial = abs(results_table.Partial_r);

% Sort by Spearman absolute correlation
[~, sort_idx] = sort(abs_spear, 'descend');

% Reorder everything
X = categorical(results_table.Metric(sort_idx));
X = reordercats(X, results_table.Metric(sort_idx));
Y = [abs_spear(sort_idx), abs_partial(sort_idx)];

% Plot grouped bars
figure; hold on;
bar(X, Y, 'grouped');

% Labels and title
legend({'|Spearman|','|Partial|'}, 'Location','best');
ylabel('Absolute Correlation (|r|)');
xlabel('Metric');
title(sprintf('Absolute Correlation with Empirical Data (T > %g, D > 20)', T));

grid on;
