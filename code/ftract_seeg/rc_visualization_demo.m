% RC network simulation + dynamic visualization (static plots & MP4)
% Patched version with transform options, auto SC/D detection, τ estimation
% ------------------------------------------------------------------------
% Before running, you typically did:
  load('SC_data_neuron_2023_yihang.mat');  % provides SC-like, D/ED, etc.
  load('idx_sc.mat');                      % optional index
%   load('Y_struct_v1-v17_PAL_res.mat');     % optional, for masks/labels
%
% Then just call:
%   rc_visualization_demo
%
% Outputs:
%   - PNG figures and MP4 movie in current folder
%   - rc_outputs_minimal.mat with I_node_max, TTP, Tau, HalfLife, features

%% ------------------ User parameters ------------------
time_scale     = 200;           % ms
num_points     = 2000;          % time steps
use_initial_pulse = true;       % true: V0(stim)=1; false: I_input(stim)=1
C_value        = 1.0;           % uniform capacitance (scalar -> vector internally)
STIM_SET       = [];            % [] = all nodes, else e.g., 1:100 or [5 17 123]
MAKE_MOVIE     = true;
%MOVIE_NODES     = [279];
MOVIE_STRIDE   = 10;
%%

% if you have certain nodes to plot together:
% want = {
%  'L_p9-46v_ROI','L_46_ROI','L_8Av_ROI','L_8C_ROI', ...
%  'L_FEF_ROI','L_LIPd_ROI','L_LIPv_ROI','L_IPS1_ROI', ...
%  'L_4_ROI','L_6d_ROI','L_6a_ROI','L_i6-8_ROI','L_SCEF_ROI', ...
%  'L_AAIC_ROI','L_AVI_ROI','L_24dd_ROI','L_24dv_ROI','L_d32_ROI', ...
%  'L_RSC_ROI','L_PreS_ROI','L_POS1_ROI','L_31a_ROI','L_v23ab_ROI', ...
%  'L_OFC_ROI','L_47m_ROI','L_25_ROI','L_pOFC_ROI', ...
%  'L_A1_ROI','L_STGa_ROI','L_TE1p_ROI','L_45_ROI','L_44_ROI', ...
%  'L_V1_ROI','L_V3A_ROI','L_V4_ROI','L_MT_ROI', ...
%  'Left-Hippocampus','Left-Amygdala'
% };
want = {
 'L_p9-46v_ROI', ...
 'L_FEF_ROI', ...
 'L_4_ROI', ...
 'L_AAIC_ROI', ...
 'L_RSC_ROI', ...
 'L_OFC_ROI',...
 'L_A1_ROI', ...
 'L_V1_ROI', ...
 'Left-Hippocampus',
};

idx = cellfun(@(s) find(strcmp(NODE_NAMES, s), 1, 'first'), want, 'uni', 1);
idx = idx(~cellfun('isempty', num2cell(idx)));       % drop missing (if any)
MOVIE_NODES = unique(idx);                               % plug into your script
%%



% Visualization transform
TRANSFORM.name = 'log1p';       % 'linear' | 'log1p' | 'gamma' | 'zscore' | 'symlog'
TRANSFORM.gamma = 0.5;          % used when name='gamma'
TRANSFORM.symlog_c = [];        % if [], auto-set from data (90th %ile)

SMOOTH.win     = 0;             % 0=no smoothing; else integer window for movmean/movmedian
SMOOTH.kind    = 'mean';        % 'mean' or 'median'

HEAT_PCT       = 99;            % robust upper caxis percentile
TOP_K_TRACES   = 10;             % for traces panel
SNAP_EDGES     = 150;           % #edges in snapshot panel

TOP_K_LABELS_DYNAMIC = 10;   % how many names to show per time step (movie & snapshots)
EXCLUDE_STIM_IN_TOPK  = true; % do not label the stim node among “top receivers”

rng(17);                        % reproducible


%% ------------------ Time vector ------------------
t = linspace(0, time_scale, num_points).';  % column, in ms

%% ------------------ Auto-detect SC and D ------------------
% Expect either SC or Test_X_W / SC_mean; D or ED / D_matrix, etc.
% Try common names in the base workspace (or current function workspace)
vars = whos;
have = @(s) any(strcmp({vars.name}, s));

if ~have('SC')
    if have('Test_X_W'), SC = evalin('base','Test_X_W'); else, try SC = Test_X_W; catch, end
    end
end
if ~exist('SC','var') && have('SC_mean'), SC = evalin('base','SC_mean'); end

if ~have('D')
    if have('ED'),       D = evalin('base','ED');
    elseif have('D_matrix'), D = evalin('base','D_matrix');
    else, try D = ED; catch, end
    end
end

% If still missing, attempt to pull from base workspace explicitly
if ~exist('SC','var')
    try SC = evalin('base','SC'); catch, end
end
if ~exist('D','var')
    try D = evalin('base','D'); catch, end
end

assert(exist('SC','var')==1 && exist('D','var')==1, ...
    ['SC and D must be available. Try:\n' ...
    '  load(''SC_data_neuron_2023_yihang.mat'');\n' ...
    '  SC = Test_X_W;  %% or SC_mean\n' ...
    '  D  = ED;        %% or D_matrix\n']);


%% ------------------ Optional ordering/index ------------------
n = size(SC,1);
if have('idx_sc')
    try
        idx_sc = evalin('base','idx_sc');
    catch
        idx_sc = (1:n).';
    end
else
    idx_sc = (1:n).';
end
SC = SC(idx_sc, idx_sc);
D  = D(idx_sc,  idx_sc);

% Your 360 cortical centroids + names in Amp order:
centroid_mni     = mni_cort;     % [360 x 3]
NODE_NAMES_reord = NODE_NAMES;    % {360 x 1} same order as mni_cort
% centroid_mni: 360x3  (cortex only, same order as NODE_NAMES_reord_cort)
% NODE_NAMES_reord_cort: 360x1 names (cortex only, Amp order)
% NODE_NAMES: n x 1 names in the current script order (n=360 or 374)

coords_cur = nan(numel(NODE_NAMES), 3);

% Use your 360-name cortical list here:
names_cort = NODE_NAMES_reord;         % <-- ensure this is the 360-only list

[tf, loc] = ismember(NODE_NAMES, names_cort);

% keep only indices that land within centroid_mni
idx_valid = tf & loc >= 1 & loc <= size(centroid_mni,1);

coords_cur(idx_valid, :) = centroid_mni(loc(idx_valid), :);

coords_ok = sum(idx_valid) >= 350;     % sanity check (most cortex mapped)


% 
% %% ------------------ MNI coordinates: align to current NODE_NAMES order ------------------
% % Inputs you already have:
% %   centroid_mni        % (360 x 3) in Amp list order (NODE_NAMES_reord)
% %   NODE_NAMES_reord    % (360 x 1) names, same order as centroid_mni
% %   NODE_NAMES          % (n x 1)  names after idx_sc reordering in THIS script
% %                        (n can be 360 or 374 if subcortex present)
% 
% USE_REAL_COORDS = true;      % master toggle
% COORD_PROJ      = 'axial';   % 'axial' (X-Y), 'sagittal' (Y-Z), 'coronal' (X-Z), or '3d'
% 
% coords_ok = false;
% coords_cur = [];   % will become (n x 3) aligned to current NODE_NAMES
% 
% try
%     % Map Amp-order names -> current order (by exact name match)
%     [tf_map, loc_map] = ismember(NODE_NAMES, NODE_NAMES_reord);
%     % Build coords array aligned to current order
%     % coords_cur = nan(numel(NODE_NAMES), 3);
%     % coords_cur(1:size(centroid_mni,1),:) = centroid_mni;
% 
%     % If you have >360 nodes (e.g., subcortex), the extra names won’t map; they remain NaN
%     % That’s fine: we’ll just skip those in spatial plots.
% 
%     % Sanity: require at least all cortical to be present
%     have_360 = sum(tf_map(1:min(360,numel(tf_map)))) >= min(360, numel(tf_map));
%     coords_ok = have_360;
% catch
%     coords_ok = false;
% end




% Try F if present (not required)
if have('F')
    F = evalin('base','F');
    F = F(idx_sc, idx_sc);
else
    F = [];
end

%% ------------------ Clean numerics & build R/G/L ------------------
SC(~isfinite(SC)) = 0;
D(~isfinite(D))   = 0;

mask = (SC>0) & (D>0);
R = zeros(n);
R(mask) = D(mask) ./ SC(mask);
R(~isfinite(R)) = 0;

G = zeros(n);
G(mask) = 1 ./ R(mask);              % conductance

L = diag(sum(G,2)) - G;              % graph Laplacian (conductance form)
C = C_value * ones(n,1) .* 0.5;

%% ------------------ Optional Y_struct pickup ------------------
Y_struct = [];
if have('Y_struct')
    Y_struct = evalin('base','Y_struct');
end
if ~isempty(Y_struct) && isfield(Y_struct,'v5')
    Amp_sorted  = Y_struct.v5.A; %#ok<NASGU>
    Prob_sorted = Y_struct.v5.P; %#ok<NASGU>
    % Ns, Ni, Nr etc. available if needed
end

%% ------------------ Stim selection ------------------
if isempty(STIM_SET)
    STIM_SET = 1:n;
end
n_stims = numel(STIM_SET);

%% ------------------ Storage ------------------
V_results  = zeros(n, num_points, n_stims, 'single');
I_node_t   = zeros(n, num_points, n_stims, 'single');  % total outgoing |current| per node over time
I_node_max = zeros(n, n_stims, 'single');
TTP        = nan(n, n_stims, 'single');                % time-to-peak (ms)

% Decay features (linear space)
Tau       = nan(n, n_stims, 'single');                 % tau (ms)
HalfLife  = nan(n, n_stims, 'single');                 % t1/2 (ms)

%% ------------------ Integrate per stim ------------------
fprintf('Simulating %d stim nodes (n=%d)...\n', n_stims, n);
for s_idx = 1:n_stims
    stim = STIM_SET(s_idx);

    I_input = zeros(n,1);
    V0      = zeros(n,1);
    if use_initial_pulse
        V0(stim) = 1;
    else
        I_input(stim) = 1;
    end

    ode_rhs = @(tt,V) (-L*V + I_input) ./ C;

    % Integrate
    [~, V_traj] = ode45(ode_rhs, t, V0);
    V_traj = V_traj.';                 % n x T
    V_results(:,:,s_idx) = single(V_traj);

    % Per-time node total outgoing |current| (linear)
    for k = 1:num_points
        V = V_traj(:, k);
        I_mat = G .* (V - V.');        % directed edge currents
        I_node_t(:, k, s_idx) = single(sum(abs(I_mat), 2));  % total outgoing |current|
    end

    % Max over time & TTP
    [peakVals, peakIdx]  = max(I_node_t(:,:,s_idx), [], 2);
    I_node_max(:, s_idx) = peakVals;
    TTP(:, s_idx)        = t(peakIdx);

    % ---- Decay tau / half-life (linear space, log-tail fit) ----
    for i = 1:n
        % y as row -> then make column for regression
        y = double(I_node_t(i,:,s_idx));          % 1 x T
        [~, pk] = max(y);
        if pk >= num_points-5 || y(pk) <= 0
            continue;
        end

        % Tail segments and force column vectors
        tail_t = t(pk:end);          tail_t = tail_t(:);      % [m x 1]
        tail_y = y(pk:end);          tail_y = tail_y(:);      % [m x 1]

        % Keep meaningful tail portion (avoid near-zero / NaN)
        thr = max(1e-8, 0.01*max(y));
        msk = isfinite(tail_y) & (tail_y > thr);
        tail_t = tail_t(msk);
        tail_y = tail_y(msk);

        if numel(tail_t) < 6
            continue;
        end

        % Design matrix and log response (avoid log(0))
        X = [ones(numel(tail_t),1), tail_t];      % [m x 2]
        ylog = log(max(tail_y, eps));             % [m x 1]

        % Solve: ylog ≈ a + b*t   => b ≈ -1/τ
        beta = X \ ylog;                           % 2 x 1
        b = beta(2);
        if b < 0
            Tau(i, s_idx)      = -1 / b;          % ms
            HalfLife(i, s_idx) = Tau(i, s_idx) * log(2);
        end
    end

    if mod(s_idx, max(1, round(n_stims/10)))==0
        fprintf('  ... %d/%d stims done\n', s_idx, n_stims);
    end
end
fprintf('Simulation complete.\n');

%%
for nm = 1:numel(MOVIE_NODES)
    MOVIE_STIM = MOVIE_NODES(nm);
    MOVIE_STIM_NAME = NODE_NAMES{MOVIE_STIM};
    MOVIE_STIM_NAME_SAFE = regexprep(MOVIE_STIM_NAME, '[^a-zA-Z0-9_-]', '');
%% ------------------ Static Plot 1: Max |current| (receiver × stim) ------------------
figure('Color','w','Position',[80 80 900 720]);
imagesc(I_node_max);
axis image; box off;
xlabel('Stim node'); ylabel('Receiver node');
title('Max outgoing |current| per node over time');
clim = [0, prctile(double(I_node_max(:)), HEAT_PCT)];
if isempty(clim) || clim(2)<=0, clim(2)=max(double(I_node_max(:))); end
caxis(clim); colormap(parula); colorbar;
set(gca,'FontName','Helvetica','FontSize',11);
fname = sprintf('fig_max_node_current_%03d_%s.png', ...
    MOVIE_STIM, MOVIE_STIM_NAME_SAFE);

saveas(gcf, fname);

%%
% Correlation for Prob
N_thres = N_sti;
emp_vec = Prob_sorted;
T = 100;

metrics = {
    'I_node_max', I_node_max(1:360,1:360)
    'Diff_eff',   diff_eff(1:360,1:360)
    'SI',         SI(1:360,1:360)
    'CMY',        CMY(1:360,1:360)
    'Nav_eff',    Nav_eff(1:360,1:360)
    'SPE',        SPE(1:360,1:360)
    'ED',         D_sorted(1:360,1:360)      % Euclidean distance
    'SC',         SC(1:360,1:360)     % structural connectivity
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
    mask = mask & (F_cor > 0.75);
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
N_thres = N_sti.*Prob_sorted;
emp_vec = Amp_sorted;
T = 100;

metrics = {
    'I_node_max', I_node_max(1:360,1:360)
    'Diff_eff',   diff_eff(1:360,1:360)
    'SI',         SI(1:360,1:360)
    'CMY',        CMY(1:360,1:360)
    'Nav_eff',    Nav_eff(1:360,1:360)
    'SPE',        SPE(1:360,1:360)
    'ED',         D_sorted(1:360,1:360)      % Euclidean distance
    'SC',         SC(1:360,1:360)     % structural connectivity
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
    mask = mask & (blkdiag(ones(180), ones(180)));
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
%% ------------------ Static Plot 3: Node×Time heatmap + top-K traces ------------------
s_view = clamp_index(MOVIE_STIM, STIM_SET);
It = double(I_node_t(:,:,s_view));               % linear currents
It_viz = apply_transform(It, TRANSFORM);         % transformed for visualization
It_viz = maybe_smooth(It_viz, SMOOTH);

[peakVals, ~] = max(It, [], 2);                  % choose top nodes by linear peak
[~, topIdx]    = maxk(peakVals, min(TOP_K_TRACES, n));

figure('Color','w','Position',[60 60 1200 460]);

subplot(1,2,1);
imagesc(t, 1:n, It_viz);
axis tight; set(gca,'YDir','normal');
xlabel('Time (ms)'); ylabel('Node');
ttl = sprintf('Node currents over time (stim %d) [%s]', STIM_SET(s_view), TRANSFORM.name);
title(ttl);
% ---- Highlight the stimulation node on heatmap ----
stim_node = STIM_SET(s_view);
hold on;
yl = ylim; xl = xlim;
plot([xl(1) xl(1)+0.03*(xl(2)-xl(1))], [stim_node stim_node], 'r-', 'LineWidth', 2.5); % short red tick at y
text(xl(1)-0.05*(xl(2)-xl(1)), stim_node, sprintf('→ Stim %d', stim_node), ...
    'Color',[0.85 0 0], 'FontWeight','bold', ...
    'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',10);
hold off;
cHi = prctile(It_viz(:), HEAT_PCT); if isempty(cHi) || cHi<=0, cHi = max(It_viz(:)); end
if ~isfinite(cHi) || cHi==0, cHi = 1; end
caxis([min(It_viz(:)), cHi]); colormap(parula); colorbar;
set(gca,'FontName','Helvetica','FontSize',11);

% ---- Label top-K receivers by name on the heatmap's right edge ----
hold on;
stim_node = STIM_SET(s_view);
xR = max(t);                 % right edge in time units
for ii = 1:numel(topIdx)
    r = topIdx(ii);
    text(xR + 0.02*(xR - min(t)), r, NODE_NAMES{r}, ...
        'Color',[0 0 0], 'FontSize',9, 'HorizontalAlignment','left', ...
        'VerticalAlignment','middle', 'Interpreter','none');
end
% Mark the stim row (short red tick + label)
xl = xlim;
plot([xl(1) xl(1)+0.03*(xl(2)-xl(1))], [stim_node stim_node], 'r-', 'LineWidth', 2.5);
text(xl(1)-0.05*(xl(2)-xl(1)), stim_node, sprintf('→ %s', NODE_NAMES{stim_node}), ...
    'Color',[0.85 0 0], 'FontWeight','bold', ...
    'HorizontalAlignment','right','VerticalAlignment','middle', 'Interpreter','none');
hold off;


subplot(1,2,2); hold on;
% Plot traces (choose: linear or viz)
TRACES_USE_VIZ = false;   % set true to plot log-like traces
if TRACES_USE_VIZ
    plot_data = apply_transform(It(topIdx,:), TRANSFORM);
    ylabel_str = sprintf('%s(|current|)', TRANSFORM.name);
else
    plot_data = It(topIdx,:);
    ylabel_str = '|Current|';
end
for k = 1:size(plot_data,1)
    plot(t, plot_data(k,:), 'LineWidth', 1.5);
end
xlabel('Time (ms)'); ylabel(ylabel_str);
title('Top receiver traces');
legend(NODE_NAMES(topIdx), 'Location','northeastoutside');
box off;
set(gca,'FontName','Helvetica','FontSize',11);
fname = sprintf('fig_strip_and_traces_stim_%03d_%s.png', STIM_SET(s_view),MOVIE_STIM_NAME_SAFE);
saveas(gcf, fname);
%% ------------------ Static Plot 4A: All-node current curves (simple) ------------------
s_view = clamp_index(MOVIE_STIM, STIM_SET);
It_all = double(I_node_t(:,:,s_view));  % [n x T]
figure('Color','w','Position',[60 60 1100 420]);
plot(t, It_all', 'LineWidth', 0.5);     % each row becomes a curve
xlabel('Time (ms)'); ylabel('|Current|');
title(sprintf('All-node |current| curves (stim %d)', STIM_SET(s_view)));
box off; set(gca,'FontName','Helvetica','FontSize',11);


% %% ------------------ Static Plot 4: Circular snapshot of top edges ------------------
% t_idx = max(1, round(num_points*0.25));
% V      = double(V_results(:, t_idx, s_view));
% I_mat  = double(G) .* (V - V.');
% [~, order] = sort(abs(I_mat(:)), 'descend');
% [ee_i, ee_j] = ind2sub([n n], order(1:min(SNAP_EDGES, numel(order))));
%
% theta = linspace(0, 2*pi, n+1)'; theta(end)=[];
% xy = [cos(theta), sin(theta)];
%
% figure('Color','w','Position',[80 80 560 560]); hold on; axis equal off;
% plot(xy(:,1), xy(:,2), '.', 'MarkerSize', 10);
% w = abs(I_mat(sub2ind([n n], ee_i, ee_j)));
% % Highlight the stimulation node
% stim_node = STIM_SET(s_view);
% plot(xy(stim_node,1), xy(stim_node,2), 'o', ...
%     'MarkerSize', 14, 'LineWidth', 2.5, 'Color', [1 0 0]); % red circle
% text(xy(stim_node,1)*1.1, xy(stim_node,2)*1.1, sprintf('Stim %d', stim_node), ...
%     'Color', [0.8 0 0], 'FontWeight', 'bold', 'HorizontalAlignment','center');
% % ---- Dynamic labels at this snapshot time ----
% It_now = double(I_node_t(:, t_idx, s_view));               % |current| at this time
% stim_node = STIM_SET(s_view);
%
% % select top-K at this time (optionally exclude stim)
% It_temp  = It_now;
% if EXCLUDE_STIM_IN_TOPK, It_temp(stim_node) = -inf; end
% Klab = min(TOP_K_LABELS_DYNAMIC, n);
% [~, curIdx] = maxk(It_temp, Klab);
%
% % place names slightly outside the circle
% r_offset = 1.08;
% for u = 1:numel(curIdx)
%     node = curIdx(u);
%     text(r_offset*xy(node,1), r_offset*xy(node,2), NODE_NAMES{node}, ...
%         'HorizontalAlignment','center','VerticalAlignment','middle', ...
%         'FontSize',9,'Color',[0 0 0], 'Interpreter','none');
% end
%
% % (keep) label the stim node in bold red
% text(1.12*xy(stim_node,1), 1.12*xy(stim_node,2), NODE_NAMES{stim_node}, ...
%     'Color',[0.85 0 0], 'FontWeight','bold', 'HorizontalAlignment','center', ...
%     'VerticalAlignment','middle', 'Interpreter','none');
%
% if ~isempty(w), w = w ./ max(w); else, w = []; end
% for e = 1:numel(w)
%     plot(xy([ee_i(e) ee_j(e)],1), xy([ee_i(e) ee_j(e)],2), '-', ...
%         'LineWidth', 0.5 + 3*w(e), 'Color', [0 0 0 0.15 + 0.7*w(e)]);
% end
% title(sprintf('Top-%d edges by |I| at t = %.1f ms (stim %d)', numel(w), t(t_idx), STIM_SET(s_view)));
% saveas(gcf, sprintf('fig_snapshot_stim_%03d_tidx_%04d.png', STIM_SET(s_view), t_idx));
%% ------------------ Static Plot 4: Real-space snapshot of top edges (MNI) ------------------
t_idx = max(1, round(num_points*0.25));
V      = double(V_results(:, t_idx, s_view));
I_mat  = double(G) .* (V - V.');
[~, order] = sort(abs(I_mat(:)), 'descend');
[ee_i, ee_j] = ind2sub([n n], order(1:min(SNAP_EDGES, numel(order))));
w = abs(I_mat(sub2ind([n n], ee_i, ee_j)));
if ~isempty(w), w = w ./ max(w); end

if USE_REAL_COORDS && coords_ok
    % --- Choose projection
    switch lower(COORD_PROJ)
        case 'axial'    % X (L–R), Y (A–P)
            XX = coords_cur(:,1); YY = coords_cur(:,2); ZZ = coords_cur(:,3);
            use3d = false; xlab='X (MNI)'; ylab='Y (MNI)';
        case 'sagittal' % Y (A–P), Z (I–S)
            XX = coords_cur(:,2); YY = coords_cur(:,3); ZZ = coords_cur(:,1);
            use3d = false; xlab='Y (MNI)'; ylab='Z (MNI)';
        case 'coronal'  % X (L–R), Z (I–S)
            XX = coords_cur(:,1); YY = coords_cur(:,3); ZZ = coords_cur(:,2);
            use3d = false; xlab='X (MNI)'; ylab='Z (MNI)';
        otherwise       % '3d'
            XX = coords_cur(:,1); YY = coords_cur(:,2); ZZ = coords_cur(:,3);
            use3d = true;  xlab='X'; ylab='Y'; % zlab set below
    end

    % Only plot nodes with valid coordinates
    good = all(isfinite([XX YY]), 2);
    if use3d, good = good & isfinite(ZZ); end

    figure('Color','w','Position',[80 80 700 600]); hold on;
    if use3d
        plot3(XX(good), YY(good), ZZ(good), '.', 'MarkerSize', 12, 'Color',[0.1 0.1 0.1]);
    else
        plot(XX(good), YY(good), '.', 'MarkerSize', 12, 'Color',[0.1 0.1 0.1]);
    end

    % edges
    E = numel(w);
    for e = 1:E
        i = ee_i(e); j = ee_j(e);
        if ~(good(i) && good(j)), continue; end
        lw = 0.5 + 3*w(e);
        if use3d
            plot3([XX(i) XX(j)], [YY(i) YY(j)], [ZZ(i) ZZ(j)], '-', 'LineWidth', lw, 'Color', [0 0 0 0.15 + 0.7*w(e)]);
        else
            plot([XX(i) XX(j)], [YY(i) YY(j)], '-', 'LineWidth', lw, 'Color', [0 0 0 0.15 + 0.7*w(e)]);
        end
    end

    % highlight stim
    stim_node = STIM_SET(s_view);
    if good(stim_node)
        if use3d
            plot3(XX(stim_node), YY(stim_node), ZZ(stim_node), 'o', 'MarkerSize', 14, 'LineWidth', 2.5, 'Color',[1 0 0]);
            text(XX(stim_node), YY(stim_node), ZZ(stim_node), sprintf('  %s', NODE_NAMES{stim_node}), ...
                'Color',[0.85 0 0], 'FontWeight','bold','Interpreter','none');
        else
            plot(XX(stim_node), YY(stim_node), 'o', 'MarkerSize', 14, 'LineWidth', 2.5, 'Color',[1 0 0]);
            text(XX(stim_node), YY(stim_node), sprintf('  %s', NODE_NAMES{stim_node}), ...
                'Color',[0.85 0 0], 'FontWeight','bold','Interpreter','none');
        end
    end

    xlabel(xlab); ylabel(ylab);
    if use3d, zlabel('Z (MNI)'); view(45,20); grid on; end
    title(sprintf('Top edges by |I| at t = %.1f ms (stim %d) – real MNI %s', t(t_idx), STIM_SET(s_view), COORD_PROJ));
    axis equal; box on;
    fname = sprintf('fig_snapshot_real_stim_%03d_%s_tidx_%04d.png', STIM_SET(s_view), MOVIE_STIM_NAME_SAFE, t_idx);
    saveas(gcf, fname);
else
    % --- Fallback to circular layout if coords missing ---
    theta = linspace(0, 2*pi, n+1)'; theta(end)=[];
    xy = [cos(theta), sin(theta)];
    figure('Color','w','Position',[80 80 560 560]); hold on; axis equal off;
    plot(xy(:,1), xy(:,2), '.', 'MarkerSize', 10);
    for e = 1:numel(w)
        plot(xy([ee_i(e) ee_j(e)],1), xy([ee_i(e) ee_j(e)],2), '-', ...
            'LineWidth', 0.5 + 3*w(e), 'Color', [0 0 0 0.15 + 0.7*w(e)]);
    end
    stim_node = STIM_SET(s_view);
    plot(xy(stim_node,1), xy(stim_node,2), 'o', 'MarkerSize', 14, 'LineWidth', 2.5, 'Color', [1 0 0]);
    text(xy(stim_node,1)*1.1, xy(stim_node,2)*1.1, sprintf('Stim %d', stim_node), ...
        'Color', [0.8 0 0], 'FontWeight', 'bold', 'HorizontalAlignment','center');
    title(sprintf('Top-%d edges by |I| at t = %.1f ms (stim %d) – circular', numel(w), t(t_idx), STIM_SET(s_view)));
     fname = sprintf('fig_snapshot_circle_stim_%03d_%s_tidx_%04d.png', STIM_SET(s_view), MOVIE_STIM_NAME_SAFE, t_idx);
    saveas(gcf, fname );
end

%% ------------------ Movie (3-panel; uses transformed heatmap) ------------------
if MAKE_MOVIE
    s_anim_idx = clamp_index(MOVIE_STIM, STIM_SET);
    ItA   = double(I_node_t(:,:,s_anim_idx));     % linear
    ItA_v = apply_transform(ItA, TRANSFORM);      % viz transform
    ItA_v = maybe_smooth(ItA_v, SMOOTH);

    [peakValsA, ~] = max(ItA, [], 2);
    [~, topIdxA]   = maxk(peakValsA, min(6, n));

    vidName = sprintf('stim_%03d_%s.mp4', MOVIE_STIM, MOVIE_STIM_NAME_SAFE);
    vw = VideoWriter(vidName, 'MPEG-4'); vw.FrameRate = 30; open(vw);
    fig = figure('Color','w','Position',[40 40 1300 650]);

    % Precompute layout for panel 3
    theta = linspace(0, 2*pi, n+1)'; theta(end)=[];
    xy = [cos(theta), sin(theta)];

    cHi = prctile(ItA_v(:), HEAT_PCT); if isempty(cHi) || cHi<=0, cHi = max(ItA_v(:)); end
    if ~isfinite(cHi) || cHi==0, cHi = 1; end
    cLo = min(ItA_v(:));

    hHMTexts = gobjects(0);     % dynamic text handles on heatmap (subplot 1)
    hCircTxt = gobjects(0);     % dynamic text handles on circle (subplot 3)

    for k = 1:MOVIE_STRIDE:num_points
        clf;

        % Panel 1: node×time heatmap (growing window)
        subplot(1,3,1);
        imagesc(t(1:k), 1:n, ItA_v(:,1:k));
        axis tight; set(gca,'YDir','normal');
        xlabel('Time (ms)'); ylabel('Node');
        title(sprintf('Node currents (stim %d) [%s]', MOVIE_STIM, MOVIE_STIM_NAME_SAFE));
        caxis([cLo cHi]); colormap(parula); colorbar;
        set(gca,'FontName','Helvetica','FontSize',11);


        % ---- Mark stim node on heatmap ----
        stim_node = STIM_SET(s_anim_idx);
        hold on;
        yl = ylim; xl = xlim;
        plot([xl(1) xl(1)+0.03*(xl(2)-xl(1))], [stim_node stim_node], 'r-', 'LineWidth', 2.5);
        text(xl(1)-0.05*(xl(2)-xl(1)), stim_node, sprintf('→ Stim %d', stim_node), ...
            'Color',[0.85 0 0],'FontWeight','bold', ...
            'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',9);
        hold off;
        % ---- Dynamic top-name labels on the heatmap’s right margin ----
        % clean previous labels
        if ~isempty(hHMTexts), delete(hHMTexts(ishandle(hHMTexts))); end
        hHMTexts = gobjects(0);

        It_now = ItA(:, k);                                  % |current| at time k
        stim_node = STIM_SET(s_anim_idx);
        It_temp   = It_now;
        if EXCLUDE_STIM_IN_TOPK, It_temp(stim_node) = -inf; end
        Klab = min(TOP_K_LABELS_DYNAMIC, n);
        [~, curIdx] = maxk(It_temp, Klab);

        xl = xlim; xr = xl(2);
        for u = 1:numel(curIdx)
            r = curIdx(u);
            hHMTexts(end+1) = text(xr + 0.02*(xr - xl(1)), r, NODE_NAMES{r}, ...
                'Color',[0 0 0], 'FontSize',8, 'HorizontalAlignment','left', ...
                'VerticalAlignment','middle', 'Interpreter','none');
        end

        % Panel 2: top traces (choose linear or viz)
        subplot(1,3,2); hold on;
        TRACES_USE_VIZ = false;
        if TRACES_USE_VIZ
            pd = apply_transform(ItA(topIdxA,:), TRANSFORM);
            ylab = sprintf('%s(|current|)', TRANSFORM.name);
        else
            pd = ItA(topIdxA,:);
            ylab = '|Current|';
        end
        for kk = 1:size(pd,1)
            plot(t, pd(kk,:), 'LineWidth', 1.2);
        end
        yl = ylim; plot([t(k) t(k)], yl, 'k--');
        xlabel('Time (ms)'); ylabel(ylab);
        title('Top receiver traces');
        legend(NODE_NAMES(topIdxA), 'Location','northeastoutside');

        box off; set(gca,'FontName','Helvetica','FontSize',11);

        % % Panel 3: snapshot edges at current time (linear)
        % subplot(1,3,3); hold on; axis equal off;
        % Vt   = double(V_results(:, k, s_anim_idx));
        % Imat = double(G) .* (Vt - Vt.');
        % [~, ord2] = sort(abs(Imat(:)), 'descend');
        % [ii, jj] = ind2sub([n n], ord2(1:min(SNAP_EDGES, numel(ord2))));
        % plot(xy(:,1), xy(:,2), '.', 'MarkerSize', 8);
        % w2 = abs(Imat(sub2ind([n n], ii, jj)));
        % % Highlight the stimulation node
        % stim_node = STIM_SET(s_anim_idx);
        % plot(xy(stim_node,1), xy(stim_node,2), 'o', ...
        %     'MarkerSize', 14, 'LineWidth', 2.5, 'Color', [1 0 0]);
        % text(xy(stim_node,1)*1.1, xy(stim_node,2)*1.1, sprintf('Stim %d', stim_node), ...
        %     'Color', [0.8 0 0], 'FontWeight', 'bold', 'HorizontalAlignment','center');
        % % ---- Dynamic top-name labels around the circle at time k ----
        % % clean previous labels
        % if ~isempty(hCircTxt), delete(hCircTxt(ishandle(hCircTxt))); end
        % hCircTxt = gobjects(0);
        %
        % % same top-K indices at this frame (reuse It_now / stim_node if in scope)
        % It_now  = double(I_node_t(:, k, s_anim_idx));
        % It_temp = It_now;
        % if EXCLUDE_STIM_IN_TOPK, It_temp(stim_node) = -inf; end
        % Klab = min(TOP_K_LABELS_DYNAMIC, n);
        % [~, curIdx] = maxk(It_temp, Klab);
        %
        % r_offset = 1.08;
        % for u = 1:numel(curIdx)
        %     node = curIdx(u);
        %     hCircTxt(end+1) = text(r_offset*xy(node,1), r_offset*xy(node,2), NODE_NAMES{node}, ...
        %         'HorizontalAlignment','center','VerticalAlignment','middle', ...
        %         'FontSize',8,'Color',[0 0 0], 'Interpreter','none');
        % end
        %
        % % (keep) bold red name for the stim node
        % text(1.12*xy(stim_node,1), 1.12*xy(stim_node,2), NODE_NAMES{stim_node}, ...
        %     'Color',[0.85 0 0], 'FontWeight','bold', 'HorizontalAlignment','center', ...
        %     'VerticalAlignment','middle', 'Interpreter','none');
        %
        %
        % if ~isempty(w2), w2 = w2 ./ max(w2); end
        % for e = 1:numel(w2)
        %     plot(xy([ii(e) jj(e)],1), xy([ii(e) jj(e)],2), '-', ...
        %         'LineWidth', 0.5 + 3*w2(e), 'Color', [0 0 0 0.15 + 0.7*w2(e)]);
        % end
        % title(sprintf('Edges at t = %.1f ms (stim %d)', t(k), STIM_SET(s_anim_idx)));
        % set(gca,'FontName','Helvetica','FontSize',11);
        % Panel 3: real-space edges at current time
        subplot(1,3,3); hold on; box on;
        Vt   = double(V_results(:, k, s_anim_idx));
        Imat = double(G) .* (Vt - Vt.');
        [~, ord2] = sort(abs(Imat(:)), 'descend');
        [ii, jj]  = ind2sub([n n], ord2(1:min(SNAP_EDGES, numel(ord2))));
        w2 = abs(Imat(sub2ind([n n], ii, jj)));
        if ~isempty(w2), w2 = w2 ./ max(w2); end

        if USE_REAL_COORDS && coords_ok
            switch lower(COORD_PROJ)
                case 'axial',    XX = coords_cur(:,1); YY = coords_cur(:,2); use3d=false; xlab='X (MNI)'; ylab='Y (MNI)';
                case 'sagittal', XX = coords_cur(:,2); YY = coords_cur(:,3); use3d=false; xlab='Y (MNI)'; ylab='Z (MNI)';
                case 'coronal',  XX = coords_cur(:,1); YY = coords_cur(:,3); use3d=false; xlab='X (MNI)'; ylab='Z (MNI)';
                otherwise,       XX = coords_cur(:,1); YY = coords_cur(:,2); ZZ = coords_cur(:,3); use3d=true; xlab='X'; ylab='Y';
            end
            good = all(isfinite([XX YY]), 2); if exist('ZZ','var'), good = good & isfinite(ZZ); end

            if use3d
                plot3(XX(good), YY(good), ZZ(good), '.', 'MarkerSize', 8, 'Color',[0.2 0.2 0.2]);
            else
                plot(XX(good), YY(good), '.', 'MarkerSize', 8, 'Color',[0.2 0.2 0.2]);
            end
            for e = 1:numel(w2)
                i = ii(e); j = jj(e);
                if ~(good(i) && good(j)), continue; end
                lw = 0.5 + 3*w2(e);
                if use3d
                    plot3([XX(i) XX(j)], [YY(i) YY(j)], [ZZ(i) ZZ(j)], '-', 'LineWidth', lw, 'Color', [0 0 0 0.15 + 0.7*w2(e)]);
                else
                    plot([XX(i) XX(j)], [YY(i) YY(j)], '-', 'LineWidth', lw, 'Color', [0 0 0 0.15 + 0.7*w2(e)]);
                end
            end

            % stim highlight
            stim_node = STIM_SET(s_anim_idx);
            if good(stim_node)
                if use3d
                    plot3(XX(stim_node), YY(stim_node), ZZ(stim_node), 'o', 'MarkerSize', 14, 'LineWidth', 2.5, 'Color', [1 0 0]);
                    text(XX(stim_node), YY(stim_node), ZZ(stim_node), sprintf('  %s', NODE_NAMES{stim_node}), ...
                        'Color',[0.85 0 0], 'FontWeight','bold','Interpreter','none');
                else
                    plot(XX(stim_node), YY(stim_node), 'o', 'MarkerSize', 14, 'LineWidth', 2.5, 'Color', [1 0 0]);
                    text(XX(stim_node), YY(stim_node), sprintf('  %s', NODE_NAMES{stim_node}), ...
                        'Color',[0.85 0 0], 'FontWeight','bold','Interpreter','none');
                end
            end
            xlabel(xlab); ylabel(ylab);
            if use3d, zlabel('Z (MNI)'); view(45,20); grid on; end
            title(sprintf('Edges at t = %.1f ms (stim %d) – real MNI %s', t(k), STIM_SET(s_anim_idx), COORD_PROJ));
        else
            % Fallback: circle
            theta = linspace(0, 2*pi, n+1)'; theta(end)=[];
            xy = [cos(theta), sin(theta)];
            plot(xy(:,1), xy(:,2), '.', 'MarkerSize', 8);
            for e = 1:numel(w2)
                plot(xy([ii(e) jj(e)],1), xy([ii(e) jj(e)],2), '-', ...
                    'LineWidth', 0.5 + 3*w2(e), 'Color', [0 0 0 0.15 + 0.7*w2(e)]);
            end
            stim_node = STIM_SET(s_anim_idx);
            plot(xy(stim_node,1), xy(stim_node,2), 'o', 'MarkerSize', 14, 'LineWidth', 2.5, 'Color', [1 0 0]);
            text(xy(stim_node,1)*1.1, xy(stim_node,2)*1.1, sprintf('Stim %d', stim_node), ...
                'Color', [0.8 0 0], 'FontWeight','bold','HorizontalAlignment','center');
            title(sprintf('Edges at t = %.1f ms (stim %d) – circular', t(k), STIM_SET(s_anim_idx)));
        end

        drawnow;
        writeVideo(vw, getframe(fig));
    end

    close(vw);
    fprintf('Saved movie: %s\n', vidName);
end
end

% %% ------------------ Feature exports (linear space) ------------------
% % Per-stim features commonly useful for SEEG comparison:
% %   - Peak (I_node_max already)
% %   - AUC, Energy
% AUC    = trapz(t, double(I_node_t), 2);         % [n x 1 x n_stims]
% AUC    = squeeze(AUC);                           % n x n_stims
% Energy = trapz(t, double(I_node_t).^2, 2);
% Energy = squeeze(Energy);
% 
% save('rc_outputs_minimal.mat', ...
%     'I_node_max','TTP','Tau','HalfLife','AUC','Energy','STIM_SET', ...
%     '-v7.3');
% fprintf('Saved rc_outputs_minimal.mat\n');



%% ------------------ Helpers ------------------
function idx = clamp_index(target_stim, STIM_SET)
% Return index in STIM_SET that matches target_stim; if not found, first
s_idx = find(STIM_SET==target_stim, 1);
if isempty(s_idx), s_idx = 1; end
idx = s_idx;
end

function Xs = maybe_smooth(X, SMOOTH)
% Optional temporal smoothing along dim=2
if ~isfield(SMOOTH,'win') || SMOOTH.win<=0
    Xs = X; return;
end
w = SMOOTH.win;
if isfield(SMOOTH,'kind') && strcmpi(SMOOTH.kind,'median')
    Xs = movmedian(X, w, 2);
else
    Xs = movmean(X, w, 2);
end
end

function Y = apply_transform(X, T)
% Apply visualization transform to nonnegative currents X (n×T)
% X is expected >= 0; if negatives appear, symlog handles them.
name = lower(T.name);
switch name
    case 'linear'
        Y = X;

    case 'log1p'
        Y = log1p(max(X,0));

    case 'gamma'
        g = T.gamma; if isempty(g), g = 0.5; end
        Y = max(X,0).^g;

    case 'zscore'
        mu  = mean(X, 2);
        sig = std(X, 0, 2) + eps;
        Y = (X - mu) ./ sig;

    case 'symlog'
        % Allow signs; for currents from I_node_t it's nonnegative, but keep generic.
        c = T.symlog_c;
        if isempty(c)
            c = prctile(abs(X(:)), 90);
            if ~isfinite(c) || c==0, c = 1; end
        end
        Y = sign(X) .* log1p(abs(X)/c);

    otherwise
        warning('Unknown TRANSFORM "%s" -> using linear', T.name);
        Y = X;
end
end
