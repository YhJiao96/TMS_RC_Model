%% ============================================================
%  ED-based RC control analysis for F-TRACT comparison
%
%  Models:
%    RC_full_EDSC   : G = SC / ED        (your current full ED-based RC)
%    RC_SC_only     : G = SC             (strength only, no ED)
%    RC_ED_edges    : G = 1 / ED         (distance only, same SC edge skeleton)
%    RC_topology    : G = 1              (binary topology only)
%    RC_ED_fullgeom : G = 1 / ED         (fully connected geometry-only model)
%
%  Main outputs:
%    Probability barplot
%    Amplitude barplot
%
%  NOTE:
%    RC simulation output is [target x stim].
%    This script transposes it to [stim x target] for F-TRACT comparison.
%% ============================================================

clear; clc;

%% =================== Load data ===================

%load('SC_data_neuron_2023_yihang.mat');
load('idx_sc.mat');
load('ReOreder_F_tract.mat');
load('ReOrder_Caio_SC.mat');
% NEW:
load("X_struct_v1-19_FCO.mat");

SC = X_struct.v8.W;
D  = X_struct.v8.ED;
F  = X_struct.v8.FCO;

atlasDir = '/Users/yihjiao/Documents/MATLAB/F-tract-2';
parcellation = 'MNI-HCP-MMP1';

%% =================== Read F-TRACT matrices ===================

% ---------- N implantations ----------
N_file = fullfile(atlasDir, parcellation, 'N_implantations.txt.gz');
output_dir = fullfile(atlasDir, parcellation);
gunzip(N_file, output_dir);

txt_file = fullfile(output_dir, 'N_implantations.txt');
N_data_matrix = readtable(txt_file);
N_data_matrix_num = N_data_matrix(8:end,1);

N_row_strings = table2cell(N_data_matrix_num);
row_values = cellfun(@(s) sscanf(s, '%f')', N_row_strings, 'UniformOutput', false);
N_final_matrix = vertcat(row_values{:});

% ---------- Amplitude median ----------
amp_file = fullfile(atlasDir, parcellation, 'amplitude__median.txt.gz');
gunzip(amp_file, output_dir);

txt_file = fullfile(output_dir, 'amplitude__median.txt');
amp_data_matrix = readtable(txt_file);
amp_data_matrix_num = amp_data_matrix(8:end,1);

amp_row_strings = table2cell(amp_data_matrix_num);
row_values = cellfun(@(s) sscanf(s, '%f')', amp_row_strings, 'UniformOutput', false);
amp_final_matrix = vertcat(row_values{:});

% ---------- Probability ----------
prob_file = fullfile(atlasDir, parcellation, 'probability.txt.gz');
gunzip(prob_file, output_dir);

txt_file = fullfile(output_dir, 'probability.txt');
prob_data_matrix = readtable(txt_file);
prob_data_matrix_num = prob_data_matrix(8:end,1);

prob_row_strings = table2cell(prob_data_matrix_num);
row_values = cellfun(@(s) sscanf(s, '%f')', prob_row_strings, 'UniformOutput', false);
prob_final_matrix = vertcat(row_values{:});

% ---------- N stimulations ----------
N_file = fullfile(atlasDir, parcellation, 'N_stimulations.txt.gz');
gunzip(N_file, output_dir);

txt_file = fullfile(output_dir, 'N_stimulations.txt');
N_data_matrix = readtable(txt_file);
N_data_matrix_num = N_data_matrix(8:end,1);

N_row_strings = table2cell(N_data_matrix_num);
row_values = cellfun(@(s) sscanf(s, '%f')', N_row_strings, 'UniformOutput', false);
N_sti_final_matrix = vertcat(row_values{:});

%% =================== Reorder matrices ===================

sc_order_1 = ReOrder_Caio_SC;
d_order_1  = ReOrder_Caio_SC;
f_order_1  = ReOrder_Caio_SC;

amp_order_1 = ReOrder_F_tract(:,2);

Amp_sorted  = amp_final_matrix(amp_order_1, amp_order_1);
Prob_sorted = prob_final_matrix(amp_order_1, amp_order_1);

N_imp_prob = N_final_matrix(amp_order_1, amp_order_1);
N_sti      = N_sti_final_matrix(amp_order_1, amp_order_1);

% This is response count approximation
N_resp = N_sti .* Prob_sorted;
N_imp  = N_imp_prob;

SC_sorted = SC(sc_order_1, sc_order_1);
D_sorted  = D(d_order_1, d_order_1);
%F_sorted  = F(f_order_1, f_order_1);
F_sorted  = F;
SC = SC_sorted;

SC(~isfinite(SC)) = 0;
D_sorted(~isfinite(D_sorted)) = 0;
F_sorted(~isfinite(F_sorted)) = 0;

%% =================== Define model and F-TRACT sizes ===================

N_model  = size(SC_sorted,1);   % usually 374: 360 cortex + 14 subcortex
N_ftract = 360;                 % F-TRACT only has 360 cortical regions
cort_idx = 1:N_ftract;

fprintf('SC model size: %d nodes\n', N_model);
fprintf('F-TRACT comparison size: %d cortical nodes\n', N_ftract);

% Keep full 374-node matrices for RC simulation
SC_model = SC_sorted;
D_model  = D_sorted;
F_model  = F_sorted;

SC_model(~isfinite(SC_model)) = 0;
D_model(~isfinite(D_model))   = 0;
F_model(~isfinite(F_model))   = 0;

% F-TRACT empirical matrices are 360 x 360
Amp_sorted  = Amp_sorted(1:N_ftract,1:N_ftract);
Prob_sorted = Prob_sorted(1:N_ftract,1:N_ftract);
N_sti       = N_sti(1:N_ftract,1:N_ftract);
N_resp      = N_resp(1:N_ftract,1:N_ftract);
N_imp       = N_imp(1:N_ftract,1:N_ftract);

% Cortical-cortical blocks for validation masks
D_cor = D_model(cort_idx,cort_idx);
F_cor = F_model(cort_idx,cort_idx);

%% ============================================================
%  RC ED-control models
%% ============================================================

fprintf('\n===== Running ED-based RC control models =====\n');

N = N_model;

W  = SC_model;
ED = D_model;

W(~isfinite(W)) = 0;
ED(~isfinite(ED)) = 0;

W(1:N+1:end) = 0;
ED(1:N+1:end) = 0;

edge_mask = (W > 0) & (ED > 0);

Cnode = ones(N,1);

time_scale = 200;      % ms
num_points = 2000;     % for final, can increase to 2000 or 5000

%% ---------- Define conductance matrices ----------

% Full ED/SC RC model: R = ED / SC, so G = SC / ED
G_RC_full_ED = zeros(N,N);
G_RC_full_ED(edge_mask) = W(edge_mask) ./ ED(edge_mask);

% SC-only RC model: R = 1 / SC, so G = SC
G_RC_SC_only = zeros(N,N);
G_RC_SC_only(edge_mask) = W(edge_mask);

% ED-only RC on real SC edges: R = ED, so G = 1 / ED
G_RC_ED_edges = zeros(N,N);
G_RC_ED_edges(edge_mask) = 1 ./ ED(edge_mask);

% Binary topology-only RC: R = 1 on real SC edges, so G = 1
G_RC_topology = zeros(N,N);
G_RC_topology(edge_mask) = 1;

% Fully connected geometry-only model: G = 1 / ED for all off-diagonal pairs
G_RC_ED_fullgeom = zeros(N,N);
geom_mask = (ED > 0) & ~eye(N);
G_RC_ED_fullgeom(geom_mask) = 1 ./ ED(geom_mask);

%% ---------- Prepare conductance matrices ----------

G_RC_full_ED     = prep_G(G_RC_full_ED);
G_RC_SC_only     = prep_G(G_RC_SC_only);
G_RC_ED_edges    = prep_G(G_RC_ED_edges);
G_RC_topology    = prep_G(G_RC_topology);
G_RC_ED_fullgeom = prep_G(G_RC_ED_fullgeom);

%% ---------- Run RC simulations ----------

fprintf('\nRunning RC full ED/SC model...\n');
I_node_max_RC_full_ED_raw = run_RC_Imax_from_G(G_RC_full_ED, Cnode, time_scale, num_points);

fprintf('\nRunning RC SC-only model...\n');
I_node_max_RC_SC_only_raw = run_RC_Imax_from_G(G_RC_SC_only, Cnode, time_scale, num_points);

fprintf('\nRunning RC ED-only on SC edges model...\n');
I_node_max_RC_ED_edges_raw = run_RC_Imax_from_G(G_RC_ED_edges, Cnode, time_scale, num_points);

fprintf('\nRunning RC topology-only model...\n');
I_node_max_RC_topology_raw = run_RC_Imax_from_G(G_RC_topology, Cnode, time_scale, num_points);

fprintf('\nRunning RC full-geometry ED-only model...\n');
I_node_max_RC_ED_fullgeom_raw = run_RC_Imax_from_G(G_RC_ED_fullgeom, Cnode, time_scale, num_points);

%% ---------- Orientation check ----------
%% ---------- Orientation check using cortical-cortical block only ----------

tmp374 = I_node_max_RC_full_ED_raw;

% RC raw is [target x stim], size 374 x 374
% F-TRACT is [stim x target], size 360 x 360
% Therefore compare only cortical-cortical block: 1:360, 1:360

tmp_noT_360 = tmp374(cort_idx, cort_idx);      % still [target x stim]
tmp_T_360   = tmp374(cort_idx, cort_idx).';    % now [stim x target]

r_noT = corr(Prob_sorted(:), tmp_noT_360(:), ...
    'type','Spearman', 'rows','complete');

r_T = corr(Prob_sorted(:), tmp_T_360(:), ...
    'type','Spearman', 'rows','complete');

fprintf('\nOrientation check using Prob_sorted cortical block:\n');
fprintf('  RC raw cortical block       [target x stim] vs Prob: rho = %.4f\n', r_noT);
fprintf('  RC transpose cortical block [stim x target] vs Prob: rho = %.4f\n', r_T);


% Usually F-TRACT is [stim x target], while RC raw is [target x stim].
% Therefore we use transpose by default.
TRANSPOSE_RC_FOR_FTRACT = true;

if TRANSPOSE_RC_FOR_FTRACT
    I_node_max_RC_full_ED_374     = I_node_max_RC_full_ED_raw.';
    I_node_max_RC_SC_only_374     = I_node_max_RC_SC_only_raw.';
    I_node_max_RC_ED_edges_374    = I_node_max_RC_ED_edges_raw.';
    I_node_max_RC_topology_374    = I_node_max_RC_topology_raw.';
    I_node_max_RC_ED_fullgeom_374 = I_node_max_RC_ED_fullgeom_raw.';
else
    I_node_max_RC_full_ED_374     = I_node_max_RC_full_ED_raw;
    I_node_max_RC_SC_only_374     = I_node_max_RC_SC_only_raw;
    I_node_max_RC_ED_edges_374    = I_node_max_RC_ED_edges_raw;
    I_node_max_RC_topology_374    = I_node_max_RC_topology_raw;
    I_node_max_RC_ED_fullgeom_374 = I_node_max_RC_ED_fullgeom_raw;
end

% Keep only cortical stim x cortical target block for F-TRACT comparison
I_node_max_RC_full_ED     = I_node_max_RC_full_ED_374(cort_idx,cort_idx);
I_node_max_RC_SC_only     = I_node_max_RC_SC_only_374(cort_idx,cort_idx);
I_node_max_RC_ED_edges    = I_node_max_RC_ED_edges_374(cort_idx,cort_idx);
I_node_max_RC_topology    = I_node_max_RC_topology_374(cort_idx,cort_idx);
I_node_max_RC_ED_fullgeom = I_node_max_RC_ED_fullgeom_374(cort_idx,cort_idx);

save('RC_ED_control_models.mat', ...
    'I_node_max_RC_full_ED', ...
    'I_node_max_RC_SC_only', ...
    'I_node_max_RC_ED_edges', ...
    'I_node_max_RC_topology', ...
    'I_node_max_RC_ED_fullgeom', ...
    'G_RC_full_ED', ...
    'G_RC_SC_only', ...
    'G_RC_ED_edges', ...
    'G_RC_topology', ...
    'G_RC_ED_fullgeom', ...
    'time_scale', ...
    'num_points', ...
    'TRANSPOSE_RC_FOR_FTRACT');

fprintf('\n===== RC ED-control models completed and saved. =====\n');

%% ============================================================
%  Graph-theoretical comparison models
%% ============================================================

fprintf('\n===== Computing graph-theoretical models =====\n');

SC_graph = SC_model;
SC_graph(~isfinite(SC_graph)) = 0;
SC_graph(SC_graph < 0) = 0;
SC_graph(1:N+1:end) = 0;

% Diffusion efficiency
try
    MFPT = mean_first_passage_time(SC_graph);
    diff_eff = 1 ./ MFPT;
    diff_eff(~isfinite(diff_eff)) = 0;
catch
    warning('mean_first_passage_time failed. Setting diff_eff to NaN.');
    diff_eff = nan(N,N);
end

% Search information
try
    L_for_SI = 1 ./ SC_graph;
    L_for_SI(~isfinite(L_for_SI)) = 0;
    SI = search_information(SC_graph, L_for_SI);
    SI(~isfinite(SI)) = 0;
catch
    warning('search_information failed. Setting SI to NaN.');
    SI = nan(N,N);
end

% Communicability
try
    s = diag(sum(SC_graph,1));
    CMY = expm(s^(-1/2) * SC_graph * s^(-1/2));
    CMY(~isfinite(CMY)) = 0;
catch
    warning('communicability failed. Setting CMY to NaN.');
    CMY = nan(N,N);
end

% Navigation efficiency
try
    L_nav = 1 ./ SC_graph;
    L_nav(~isfinite(L_nav)) = 0;
    [sr, PL_bin, PL_wei, PL_dis, paths] = navigation_wu(L_nav, D_sorted);
    Nav_eff = 1 ./ PL_wei;
    Nav_eff(~isfinite(Nav_eff)) = 0;
catch
    warning('navigation_wu failed. Setting Nav_eff to NaN.');
    Nav_eff = nan(N,N);
end

% Shortest path efficiency
try
    L_spl = 1 ./ SC_graph;
    L_spl(~isfinite(L_spl)) = 0;
    [SPL,hops,Pmat] = distance_wei_floyd(L_spl);
    SPE = 1 ./ SPL;
    SPE(~isfinite(SPE)) = 0;
catch
    warning('distance_wei_floyd failed. Setting SPE to NaN.');
    SPE = nan(N,N);
end

SC_cort  = SC_model(cort_idx,cort_idx);
ED_cort  = D_model(cort_idx,cort_idx);

SI       = SI(cort_idx,cort_idx);
SPE      = SPE(cort_idx,cort_idx);
CMY      = CMY(cort_idx,cort_idx);
Nav_eff  = Nav_eff(cort_idx,cort_idx);
diff_eff = diff_eff(cort_idx,cort_idx);

fprintf('===== Graph-theoretical models completed =====\n');

%% ============================================================
%  Plot configuration
%% ============================================================

T_pairs = 100;
D_min   = 20;
same_hemi = blkdiag(ones(180), ones(180)) > 0;

B = 1000;
rng(17);

STYLE = struct( ...
    'figPixels',   [100 100 2600 700], ...
    'titleFS',     20, ...
    'suptitleFS',  22, ...
    'axisFS',      16, ...
    'labelFS',     18, ...
    'legendFS',    15, ...
    'tickAngle',   25, ...
    'barWidth',    0.9, ...
    'errDX',       0.08, ...
    'catGap',      0.10 ...
    );

Cmap = containers.Map;

Cmap('RC_full_EDSC')   = [0.8392 0.1529 0.1569];
Cmap('RC_SC_only')     = [1.0000 0.4980 0.0549];
Cmap('RC_ED_edges')    = [0.1216 0.4667 0.7059];
Cmap('RC_topology')    = [0.5804 0.4039 0.7412];
Cmap('RC_ED_fullgeom') = [0.6000 0.6000 0.6000];

Cmap('SC')             = [0.5490 0.3373 0.2941];
Cmap('ED')             = [0.9294 0.6941 0.1255];
Cmap('SI')             = [0.1216 0.4667 0.7059];
Cmap('SPE')            = [0.0902 0.7451 0.8118];
Cmap('CMY')            = [0.5804 0.4039 0.7412];
Cmap('Nav_eff')        = [1.0000 0.4980 0.0549];
Cmap('Diff_eff')       = [0.1725 0.6275 0.1725];

%% ============================================================
%  Model list
%% ============================================================

METRICS = {
    'RC_full_EDSC',    I_node_max_RC_full_ED
    'RC_SC_only',      I_node_max_RC_SC_only
    'RC_ED_edges',     I_node_max_RC_ED_edges
    'RC_topology',     I_node_max_RC_topology
    'RC_ED_fullgeom',  I_node_max_RC_ED_fullgeom
    'SC',              SC_cort
    'ED',              ED_cort
    'SI',              SI
    'SPE',             SPE
    };

N = N_ftract;
%same_hemi = blkdiag(ones(180), ones(180)) > 0;

D_cor = D_model(cort_idx,cort_idx);
F_cor = F_model(cort_idx,cort_idx);

%% =================== Mask builders ===================

build_base_mask = @(N_thres) ...
    (N_thres >= T_pairs) & ...
    (N_imp >= 5) & ...
    ~eye(N_ftract) & ...
    (D_cor >= D_min);

mask_all = @(N_thres) build_base_mask(N_thres);

mask_connected = @(N_thres) ...
    build_base_mask(N_thres) & (F_cor >= 0.75);

mask_unconnected = @(N_thres) ...
    build_base_mask(N_thres) & (F_cor <= 0.25);

maskSet = {mask_all, mask_connected, mask_unconnected};

%% =================== Probability figure ===================

N_thres_prob = N_sti;
emp_prob = Prob_sorted;

plot_triptych_signed( ...
    'Probability prediction: full vs ED-partial Spearman correlation', ...
    emp_prob, ...
    N_thres_prob, ...
    maskSet, ...
    METRICS, ...
    D_cor, ...
    Cmap, ...
    B, ...
    STYLE);

%% =================== Amplitude figure ===================

N_thres_amp = N_resp;
emp_amp = Amp_sorted;

plot_triptych_signed( ...
    'Amplitude prediction: full vs ED-partial Spearman correlation', ...
    emp_amp, ...
    N_thres_amp, ...
    maskSet, ...
    METRICS, ...
    D_cor, ...
    Cmap, ...
    B, ...
    STYLE);

%% ============================================================
%  Local functions
%% ============================================================

function G = prep_G(G)

    G(~isfinite(G)) = 0;
    G(G < 0) = 0;

    n = size(G,1);
    G(1:n+1:end) = 0;

    G = (G + G.') ./ 2;

    nz = G(G > 0);
    if ~isempty(nz)
        G = G ./ mean(nz);
    end
end


function I_node_max = run_RC_Imax_from_G(G, Cnode, time_scale, num_points)

    n = size(G,1);

    G(~isfinite(G)) = 0;
    G(G < 0) = 0;
    G(1:n+1:end) = 0;

    % Laplacian
    deg = sum(G,2);
    Lg = diag(deg) - G;

    % Capacitance
    Cnode = Cnode(:);
    Cnode(~isfinite(Cnode) | Cnode <= 0) = 1;
    Cinv = 1 ./ Cnode;

    % System matrix: dV/dt = A V
    A = -diag(Cinv) * Lg;

    % Time vector
    t_span = linspace(0, time_scale, num_points);

    % Output: [target x stim]
    I_node_max = zeros(n,n);

    % Initial conditions for all stimulation nodes at once
    % column s = stimulation at node s
    V0_all = eye(n);

    for tt = 1:num_points

        % Voltage at this time for all stimulations
        % V_all(:,s) = voltage vector for stimulation s
        V_all = expm(A * t_span(tt)) * V0_all;

        % Compute node-wise current magnitude for each stimulation
        % Loop over stim columns only
        I_abs_all = zeros(n,n);

        for s = 1:n
            v = V_all(:,s);

            % Current from i to j
            dV = v - v.';
            I_mat = G .* dV;

            % Total absolute current for each target node
            I_abs_all(:,s) = sum(abs(I_mat),2);
        end

        I_node_max = max(I_node_max, I_abs_all);

        if mod(tt,100)==0 || tt==num_points
            fprintf('  time point %d/%d done\n', tt, num_points);
        end
    end

    I_node_max(1:n+1:end) = NaN;
end

function plot_triptych_signed(fig_title, emp_mat, N_thres, maskFns, METRICS, D_cor, Cmap, B, STYLE)

    titles = {'All pairs', 'Connected pairs', 'Unconnected pairs'};

    figure('Color','w','Position',STYLE.figPixels);
    tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

    for k = 1:3
        nexttile;
        ax = gca;

        M = maskFns{k}(N_thres) & isfinite(emp_mat);

        [names, counts, rF, rP, ciF, ciP] = compute_stats(emp_mat, METRICS, D_cor, M, B);

        pretty_dualbars_signed(ax, names, rF, rP, ciF, ciP, counts, Cmap, STYLE);

        title(sprintf('%s (n = %d)', titles{k}, nnz(M)), ...
            'FontWeight','bold', ...
            'FontSize',STYLE.titleFS);

        ax.FontSize = STYLE.axisFS;
        ax.LineWidth = 1.2;
        ylabel('Spearman \rho', 'FontSize', STYLE.labelFS);
    end

    title(tl, fig_title, 'FontSize', STYLE.suptitleFS, 'FontWeight','bold');
end


function [names, counts, r_full, r_part, ci_full, ci_part] = compute_stats(emp_mat, METRICS, D_cor, mask_all, B)

    m = size(METRICS,1);

    names   = strings(m,1);
    r_full  = nan(m,1);
    r_part  = nan(m,1);
    ci_full = nan(m,2);
    ci_part = nan(m,2);
    counts  = nan(m,1);

    for i = 1:m

        names(i) = string(METRICS{i,1});
        sim = METRICS{i,2};

        M = mask_all & isfinite(sim);

        emp  = emp_mat(M);
        simv = sim(M);
        Dv   = D_cor(M);

        counts(i) = numel(emp);

        if counts(i) < 5
            continue;
        end

        oF = boot_corr(emp, simv, [], ...
            'Partial', false, ...
            'B', B);

        oP = boot_corr(emp, simv, Dv, ...
            'Partial', true, ...
            'B', B);

        r_full(i) = oF.r_hat;
        r_part(i) = oP.r_hat;

        ci_full(i,:) = oF.ci;
        ci_part(i,:) = oP.ci;
    end

    % Sort by absolute full correlation
    [~, idx] = sort(abs(r_full), 'descend', 'MissingPlacement','last');

    names   = names(idx);
    counts  = counts(idx);
    r_full  = r_full(idx);
    r_part  = r_part(idx);
    ci_full = ci_full(idx,:);
    ci_part = ci_part(idx,:);
end


function pretty_dualbars_signed(ax, names, rF, rP, ciF, ciP, counts, Cmap, STYLE)

    y1 = rF;
    y2 = rP;

    n = numel(y1);

    raw = cellstr(names);
    dispLabels = raw;

    for j = 1:n
        if strcmp(raw{j}, 'RC_full_EDSC')
            dispLabels{j} = 'RC full';
        elseif strcmp(raw{j}, 'RC_SC_only')
            dispLabels{j} = 'RC SC';
        elseif strcmp(raw{j}, 'RC_ED_edges')
            dispLabels{j} = 'RC ED-edge';
        elseif strcmp(raw{j}, 'RC_topology')
            dispLabels{j} = 'RC topology';
        elseif strcmp(raw{j}, 'RC_ED_fullgeom')
            dispLabels{j} = 'RC geometry';
        elseif strcmp(raw{j}, 'I_node_max')
            dispLabels{j} = 'RC model';
        end
    end

    col = repmat([0.5 0.5 0.5], n, 1);

    if isa(Cmap, 'containers.Map')
        for j = 1:n
            key = raw{j};
            if isKey(Cmap, key)
                val = Cmap(key);
                if isnumeric(val) && numel(val) == 3
                    col(j,:) = val(:).';
                end
            end
        end
    end

    x = 1:n;
    width = max(0.05, 1 - STYLE.catGap);

    hold(ax,'on');

    b1 = bar(ax, x, y1, width, ...
        'FaceColor','flat', ...
        'EdgeColor','none', ...
        'FaceAlpha',0.30);

    b1.CData = col;

    b2 = bar(ax, x, y2, width, ...
        'FaceColor','flat', ...
        'EdgeColor','none', ...
        'FaceAlpha',0.95);

    b2.CData = col;

    uistack(b2,'top');

    loF = rF - ciF(:,1);
    hiF = ciF(:,2) - rF;

    loP = rP - ciP(:,1);
    hiP = ciP(:,2) - rP;

    x1 = b1.XEndPoints;
    x2 = b2.XEndPoints;
    dx = STYLE.errDX;

    errorbar(x1-dx, y1, loF, hiF, ...
        'k', ...
        'LineStyle','none', ...
        'LineWidth',1);

    errorbar(x2+dx, y2, loP, hiP, ...
        'k', ...
        'LineStyle','none', ...
        'LineWidth',1);

    yline(ax, 0, 'k-', 'LineWidth', 1);

    set(ax, 'TickLabelInterpreter','none');
    xticks(ax, x);
    xticklabels(ax, dispLabels);
    xtickangle(ax, STYLE.tickAngle);

    grid(ax, 'on');
    ax.GridAlpha = 0.2;
    box(ax, 'on');

    xlim(ax, [0.5, n+0.5]);

    finiteVals = [y1(:); y2(:); ciF(:); ciP(:)];
    finiteVals = finiteVals(isfinite(finiteVals));

    if isempty(finiteVals)
        ylim(ax, [-0.1 0.1]);
    else
        ymax = max(abs(finiteVals));
        ylim(ax, [-1.15*ymax, 1.15*ymax]);
    end

    FULL_SW = [1.0 0.6 0.6];
    PART_SW = [0.8 0.0 0.0];

    hFull = patch(nan, nan, FULL_SW, ...
        'EdgeColor','none', ...
        'FaceAlpha',0.30);

    hPart = patch(nan, nan, PART_SW, ...
        'EdgeColor','none', ...
        'FaceAlpha',0.95);

    lg = legend([hFull, hPart], ...
        {'Full', 'ED-partial'}, ...
        'Location','northeast');

    lg.Box = 'off';
    lg.FontSize = STYLE.legendFS;
end


function out = boot_corr(emp, sim, dist, varargin)

    p = inputParser;

    addParameter(p, 'Partial', false, @(x)islogical(x)&&isscalar(x));
    addParameter(p, 'B', 1000, @(x)isnumeric(x)&&isscalar(x)&&x>10);

    parse(p, varargin{:});

    isPartial = p.Results.Partial;
    B = p.Results.B;

    emp = emp(:);
    sim = sim(:);

    if isPartial
        dist = dist(:);
        valid = isfinite(emp) & isfinite(sim) & isfinite(dist);
    else
        valid = isfinite(emp) & isfinite(sim);
    end

    emp = emp(valid);
    sim = sim(valid);

    if isPartial
        dist = dist(valid);
    end

    n = numel(emp);

    if n < 5
        out = struct( ...
            'r_hat', NaN, ...
            'ci', [NaN NaN], ...
            'se', NaN, ...
            'p_boot', NaN, ...
            'r_boot', []);
        return;
    end

    if isPartial
        r_hat = partialcorr(emp, sim, dist, 'type','Spearman');
    else
        r_hat = corr(emp, sim, 'type','Spearman', 'rows','complete');
    end

    idxMat = randi(n, [n, B]);
    r_boot = nan(B,1);

    for b = 1:B
        I = idxMat(:,b);

        if isPartial
            r_boot(b) = partialcorr(emp(I), sim(I), dist(I), 'type','Spearman');
        else
            r_boot(b) = corr(emp(I), sim(I), 'type','Spearman', 'rows','complete');
        end
    end

    out.ci = quantile(r_boot, [0.025 0.975]);
    out.se = std(r_boot, 1);
    out.p_boot = mean(abs(r_boot) >= abs(r_hat));
    out.r_hat = r_hat;
    out.r_boot = r_boot;
end

