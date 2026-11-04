
clear all;
% %% -------------------- Scatter Plot: TTP vs Empirical Delay --------------------
%% ------------------ User parameters ------------------
time_scale     = 200;           % ms
num_points     = 2000;          % time steps
use_initial_pulse = true;       % true: V0(stim)=1; false: I_input(stim)=1
C_value        = 1.0;           % uniform capacitance (scalar -> vector internally)
STIM_SET       = [];            % [] = all nodes, else e.g., 1:100 or [5 17 123]
MAKE_MOVIE     = true;
MOVIE_STIM     = 279;
MOVIE_STRIDE   = 10;

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

% load('SC_data_neuron_2023_yihang.mat');
% load('idx_sc.mat');
atlasDir = '/Users/yihjiao/Documents/MATLAB/F-tract-2';
parcellation = 'MNI-HCP-MMP1';
%parcellation = 'Lausanne2008-500';

% 
% 


%% ------------------ Time vector ------------------
t = linspace(0, time_scale, num_points).';  % column, in ms

load('SC_data_neuron_2023_yihang.mat');
load('idx_sc.mat');
atlasDir = '/Users/yihjiao/Documents/MATLAB/F-tract-2';
parcellation = 'MNI-HCP-MMP1';
%parcellation = 'Lausanne2008-500';

% % SC =  X_struct.v8.W;
% % D = X_struct.v8.ED;
% % F = X_struct.v8.SCP;

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

amp_file = fullfile(atlasDir,parcellation,'amplitude__median.txt.gz');
output_dir = fullfile(atlasDir,parcellation);  % You can also use atlasDir if you prefer
gunzip(amp_file, output_dir);

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

N_file = fullfile(atlasDir,parcellation,'peak_delay__median.txt.gz');
output_dir = fullfile(atlasDir,parcellation);  % You can also use atlasDir if you prefer
gunzip(N_file, output_dir);

txt_file = fullfile(output_dir, 'peak_delay__median.txt');
N_data_matrix = readtable(txt_file); 
N_data_matrix_num = N_data_matrix(8:end,1); 

% Convert table (360x1) of strings to cell array of strings
N_row_strings = table2cell(N_data_matrix_num);  % now 360x1 cell

% Parse each row string into a numeric row
row_values = cellfun(@(s) sscanf(s, '%f')', N_row_strings, 'UniformOutput', false);

% Combine into final matrix
delay_peak_matrix = vertcat(row_values{:});  % should now be 360 x 360
% 
% max(max(delay_peak_matrix))
% figure;imagesc(delay_peak_matrix);
%%
idx_sc = 1:374;
sc_order_1 = idx_sc; 
d_order_1 = idx_sc;
f_order_1 = idx_sc;
emp_idx = idx_sc(1:360);

Amp_sorted = amp_final_matrix;
Prob_sorted = prob_final_matrix;

N_prob = prob_final_matrix;
N_sti = N_sti_final_matrix;
N_imp = N_sti.*Prob_sorted;


amplitude_value_cor = Amp_sorted;
prob_value_cor = Prob_sorted;

SC_sorted = SC(sc_order_1, sc_order_1);
D_sorted = D(d_order_1,d_order_1);
F_sorted = F(f_order_1,f_order_1);
SC = SC_sorted; % strength
D = D_sorted;
SC(F_sorted<0.25)=0;

%%

SC(isnan(SC) | isinf(SC)) = 0;
D_sorted(isnan(D_sorted) | isinf(D_sorted)) = 0;



%%
% load('SC_data_neuron_2023_yihang.mat');
% load('idx_sc.mat');
% load("Y_struct_v1-v17_PAL_res.mat");
% load("X_struct_v1-19_FCO.mat");
% amp_final_matrix = Y_struct.v5.A;
% prob_final_matrix = Y_struct.v5.P;
% emp_idx = Y_struct.v5.parcel_list_sidx;
% idx_sc = 1:374;
% sc_order_1 = idx_sc; 
% d_order_1 = idx_sc;
% f_order_1 = idx_sc;
% % SC =  X_struct.v8.W;
% % D = X_struct.v8.ED;
% % F = X_struct.v8.SCP;
% Amp_sorted = amp_final_matrix;
% Prob_sorted = prob_final_matrix;
% 
% N_prob = Y_struct.v5.Ni;
% N_sti = Y_struct.v5.Ns;
% N_imp = Y_struct.v5.Ni;
% N_resp = Y_struct.v5.Nr; 
% 
% 
% amplitude_value_cor = Amp_sorted;
% prob_value_cor = Prob_sorted;
% 
% % SC_sorted = SC(sc_order_1, sc_order_1);
% % D_sorted = D(d_order_1,d_order_1);
% % F_sorted = F(f_order_1,f_order_1);
% SC = SC_sorted; % strength
% D_sorted = D;
% F_sorted = F;
% %SC(F_sorted<0.25)=0;
% 
% %
% 
% SC(isnan(SC) | isinf(SC)) = 0;
% D_sorted(isnan(D_sorted) | isinf(D_sorted)) = 0;

% Avoid division by zero
R =  D_sorted./ SC;
R(isnan(R) | isinf(R)) = 0;


% 
% N_file = fullfile(atlasDir,parcellation,'N_implantations.txt.gz');
% output_dir = fullfile(atlasDir,parcellation);  % You can also use atlasDir if you prefer
% gunzip(N_file, output_dir);
% 
% txt_file = fullfile(output_dir, 'N_implantations.txt');
% N_data_matrix = readtable(txt_file); 
% N_data_matrix_num = N_data_matrix(8:end,1); 
% 
% % Convert table (360x1) of strings to cell array of strings
% N_row_strings = table2cell(N_data_matrix_num);  % now 360x1 cell
% 
% % Parse each row string into a numeric row
% row_values = cellfun(@(s) sscanf(s, '%f')', N_row_strings, 'UniformOutput', false);
% 
% % Combine into final matrix
% N_final_matrix = vertcat(row_values{:});  % should now be 360 x 360
% 
% amp_file = fullfile(atlasDir,parcellation,'amplitude__median.txt.gz');
% output_dir = fullfile(atlasDir,parcellation);  % You can also use atlasDir if you prefer
% gunzip(amp_file, output_dir);
% 
% txt_file = fullfile(output_dir, 'amplitude__median.txt');
% amp_data_matrix = readtable(txt_file); 
% amp_data_matrix_num = amp_data_matrix(8:end,1); 
% 
% % Convert table (360x1) of strings to cell array of strings
% amp_row_strings = table2cell(amp_data_matrix_num);  % now 360x1 cell
% 
% % Parse each row string into a numeric row
% row_values = cellfun(@(s) sscanf(s, '%f')', amp_row_strings, 'UniformOutput', false);
% 
% % Combine into final matrix
% amp_final_matrix = vertcat(row_values{:});  % should now be 360 x 360
% 
% 
% prob_file = fullfile(atlasDir,parcellation,'probability.txt.gz');
% output_dir = fullfile(atlasDir,parcellation);  % You can also use atlasDir if you prefer
% gunzip(prob_file, output_dir);
% 
% txt_file = fullfile(output_dir, 'probability.txt');
% prob_data_matrix = readtable(txt_file); 
% prob_data_matrix_num = prob_data_matrix(8:end,1); 
% 
% % Convert table (360x1) of strings to cell array of strings
% prob_row_strings = table2cell(prob_data_matrix_num);  % now 360x1 cell
% 
% % Parse each row string into a numeric row
% row_values = cellfun(@(s) sscanf(s, '%f')', prob_row_strings, 'UniformOutput', false);
% 
% % Combine into final matrix
% prob_final_matrix = vertcat(row_values{:});  % should now be 360 x 360
% 
% 
% N_file = fullfile(atlasDir,parcellation,'N_stimulations.txt.gz');
% output_dir = fullfile(atlasDir,parcellation);  % You can also use atlasDir if you prefer
% gunzip(N_file, output_dir);
% 
% txt_file = fullfile(output_dir, 'N_stimulations.txt');
% N_data_matrix = readtable(txt_file); 
% N_data_matrix_num = N_data_matrix(8:end,1); 
% 
% % Convert table (360x1) of strings to cell array of strings
% N_row_strings = table2cell(N_data_matrix_num);  % now 360x1 cell
% 
% % Parse each row string into a numeric row
% row_values = cellfun(@(s) sscanf(s, '%f')', N_row_strings, 'UniformOutput', false);
% 
% % Combine into final matrix
% N_sti_final_matrix = vertcat(row_values{:});  % should now be 360 x 360
% 
N_file = fullfile(atlasDir,parcellation,'peak_delay__median.txt.gz');
output_dir = fullfile(atlasDir,parcellation);  % You can also use atlasDir if you prefer
gunzip(N_file, output_dir);

txt_file = fullfile(output_dir, 'peak_delay__median.txt');
N_data_matrix = readtable(txt_file); 
N_data_matrix_num = N_data_matrix(8:end,1); 

% Convert table (360x1) of strings to cell array of strings
N_row_strings = table2cell(N_data_matrix_num);  % now 360x1 cell

% Parse each row string into a numeric row
row_values = cellfun(@(s) sscanf(s, '%f')', N_row_strings, 'UniformOutput', false);

% Combine into final matrix
delay_peak_matrix = vertcat(row_values{:});  % should now be 360 x 360
% 
% % % max(max(delay_peak_matrix))
% % % figure;imagesc(delay_peak_matrix);
% % %%
% % sc_order_1 = idx_sc; 
% % d_order_1 = idx_sc;
% % f_order_1 = idx_sc;
% % 
% % emp_idx = idx_sc(1:360);
% Amp_sorted = amp_final_matrix(emp_idx,emp_idx);
% Prob_sorted = prob_final_matrix(emp_idx,emp_idx);
% 
% N_prob = prob_final_matrix(emp_idx,emp_idx);
% N_sti = N_sti_final_matrix(emp_idx,emp_idx);
% N_imp = N_sti.*Prob_sorted;
% 
% 
% amplitude_value_cor = Amp_sorted;
% prob_value_cor = Prob_sorted;
% 
% SC_sorted = SC(sc_order_1, sc_order_1);
% D_sorted = D(d_order_1,d_order_1);
% F_sorted = F(f_order_1,f_order_1);
% SC = SC_sorted; % strength
% D = D_sorted;
% %D_cor = D_sorted(1:360,1:360);
% F_cor = F_sorted(1:360,1:360);
% %SC(F_sorted<0.25)=0;

%%
F_cor = F(1:360,1:360);
t = linspace(0, time_scale, num_points).';
delay_peak_matrix_re = delay_peak_matrix(emp_idx,emp_idx);
SC(isnan(SC) | isinf(SC)) = 0;
D_sorted(isnan(D_sorted) | isinf(D_sorted)) = 0;



n = size(SC,1);
mask = (SC>0) & (D>0);
R = zeros(n,n);
R(mask) = D(mask) ./ SC(mask);
R(~isfinite(R)) = 0;

G = zeros(n);
G(mask) = 1 ./ R(mask);              % conductance

L = diag(sum(G,2)) - G;              % graph Laplacian (conductance form)
C = C_value * ones(n,1) ;

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
% 4-arg (defaults):
DEL = rc_compute_delays_from_V(V_results, R, t(:), SC);

% or 5-arg with options:
DEL = rc_compute_delays_from_V(V_results, R, t(:), SC, ...
    struct('smooth_ms',3,'baseline_ms',10,'guard_ms',5, ...
    'thr_frac',0.10,'direct_thr',0,'min_area',1e-9));


% Your final RC delay choice (n x n)
T_delay_RC = DEL.chosen;   % rows=targets, cols=stims
function DEL = rc_compute_delays_from_V(V_results, R, t, SC, opts)
% RC delays from RC simulation outputs (version-agnostic: supports 4 or 5 args)
% Usage:
%   DEL = rc_compute_delays_from_V(V_results, R, t, SC);
%   DEL = rc_compute_delays_from_V(V_results, R, t, SC, struct('smooth_ms',3,...));
%
% OUTPUT fields: DEL.mtt_I, DEL.q50_I, DEL.fpt_I, DEL.centroid_V, DEL.q50_V,
%                DEL.chosen, DEL.used_method  (all n x n, rows=targets, cols=stims)

% -------- argument handling --------
if nargin < 5 || isempty(opts), opts = struct(); end
p = inputParser; p.KeepUnmatched = true;
addParameter(p, 'smooth_ms',   3);
addParameter(p, 'baseline_ms', 10);
addParameter(p, 'guard_ms',    5);
addParameter(p, 'thr_frac',    0.10);
addParameter(p, 'direct_thr',  0);
addParameter(p, 'min_area',    1e-9);
parse(p, opts);
opt = p.Results;

[n, T, n2] = size(V_results); %#ok<ASGLU>
assert(isequal(size(R), [n n]), 'R must be n×n');
assert(isequal(size(SC),[n n]), 'SC must be n×n');
t   = t(:);                  % ensure column
dt  = mean(diff(t));
w   = @(ms) max(1, round(ms/dt));
bwin  = 1:min(T, w(opt.baseline_ms));
gidx0 = find(t >= opt.guard_ms, 1, 'first'); if isempty(gidx0), gidx0 = 1; end

% Conductance & Laplacian
G = zeros(size(R)); mask = R > 0 & isfinite(R); G(mask) = 1./R(mask);
L = diag(sum(G,2)) - G;

% Pre-alloc
mtt_I      = nan(n,n);
q50_I      = nan(n,n);
fpt_I      = nan(n,n);
centroid_V = nan(n,n);
q50_V      = nan(n,n);
chosen     = nan(n,n);
used       = strings(n,n);

for s = 1:n
    % --- voltage traces for this stim ---
    V = squeeze(V_results(:,:,s));                         % n x T
    V = movmean(V, w(opt.smooth_ms), 2);
    V = max(V - median(V(:,bwin),2,'omitnan'), 0);         % baseline>=0

    % --- node inflow current: | -L * V |
    I = abs((-L) * V);                                     % n x T

    % --- mean transit time from current ---
    areaI  = trapz(t, I, 2);                               % n x 1
    numI   = trapz(t, I .* t', 2);
    mtt    = numI ./ areaI;
    mtt(areaI < opt.min_area) = NaN;


    % --- time to 50% cumulative charge (robust to plateaus) ---
    Q      = cumtrapz(t, I, 2);         % n x T
    Qend   = Q(:,end);                  % n x 1
    q50    = nan(n,1);
    for i = 1:n
        if Qend(i) > opt.min_area
            q50(i) = safe_qfrac_time(Q(i,:), t, 0.5*Qend(i));
        end
    end

    % --- first-passage time after guard window ---
    fpt    = nan(n,1);
    for i = 1:n
        thr = opt.thr_frac * max(I(i,:));
        if thr > 0
            k = find(I(i, gidx0:end) >= thr, 1, 'first');
            if ~isempty(k), fpt(i) = t(gidx0 + k - 1); end
        end
    end

    % --- voltage-based backups ---
    areaV  = trapz(t, V, 2);
    numV   = trapz(t, V .* t', 2);
    cV     = numV ./ areaV;  cV(areaV < opt.min_area) = NaN;

    % --- time to 50% cumulative voltage area (robust) ---
    AV     = cumtrapz(t, V, 2);
    AVend  = AV(:,end);
    q50v   = nan(n,1);
    for i = 1:n
        if AVend(i) > opt.min_area
            q50v(i) = safe_qfrac_time(AV(i,:), t, 0.5*AVend(i));
        end
    end


    % store
    mtt_I(:,s)      = mtt;
    q50_I(:,s)      = q50;
    fpt_I(:,s)      = fpt;
    centroid_V(:,s) = cV;
    q50_V(:,s)      = q50v;

    % choose final metric per pair
    direct = SC(:,s) > opt.direct_thr; direct(s) = false;
    for i = 1:n
        if direct(i)
            val = mtt(i);   meth = "mtt_I";
            if isnan(val), val = q50(i);  meth = "q50_I"; end
            if isnan(val), val = fpt(i);  meth = "fpt_I"; end
            if isnan(val), val = cV(i);   meth = "centroid_V"; end
        else
            val = q50(i);   meth = "q50_I";
            if isnan(val), val = fpt(i);  meth = "fpt_I"; end
            if isnan(val), val = mtt(i);  meth = "mtt_I"; end
            if isnan(val), val = cV(i);   meth = "centroid_V"; end
        end
        chosen(i,s) = val; used(i,s) = meth;
    end
end

DEL = struct('mtt_I', mtt_I, 'q50_I', q50_I, 'fpt_I', fpt_I, ...
    'centroid_V', centroid_V, 'q50_V', q50_V, ...
    'chosen', chosen, 'used_method', used);
end
function tq = safe_qfrac_time(cumvec, t, target)
% Return the earliest time when cumulative curve 'cumvec' reaches 'target'.
% Handles non-strictly-increasing curves (plateaus) without interp1.
%
% cumvec : 1 x T cumulative (non-decreasing) array
% t      : T x 1 time vector (ms)
% target : scalar target cumulative value (e.g., 0.5*cumvec(end))
%
% Rules:
%  - If target <= cumvec(1): return t(1)
%  - Else find first k with cumvec(k) >= target
%     * If k==1 -> t(1)
%     * If cumvec(k) == cumvec(k-1) (plateau): return t(k)
%     * Else linear interpolate between (k-1,k)

t = t(:);                % ensure column
c = double(cumvec(:));   % column

if ~isfinite(target) || ~any(isfinite(c))
    tq = NaN; return;
end

% Guard: monotonic non-decreasing (enforce numerically)
% (Optional; usually not needed)
% c = max(c, cummax([c(1); c(1:end-1)]));

if target <= c(1)
    tq = t(1); return;
end

k = find(c >= target, 1, 'first');
if isempty(k)
    tq = t(end); return;                 % target not reached -> end
elseif k == 1
    tq = t(1); return;
else
    c1 = c(k-1); c2 = c(k);  t1 = t(k-1); t2 = t(k);
    if c2 <= c1 || ~isfinite(c1) || ~isfinite(c2)
        % plateau or bad segment: pick first crossing time
        tq = t(k);
    else
        % linear interpolate within [k-1,k]
        alpha = (target - c1) / (c2 - c1);   % in (0,1]
        tq = t1 + alpha * (t2 - t1);
    end
end
end
%%
% Replace cutoff values with NaN
TTP(TTP >= 200) = NaN;
sim_TTP = DEL.q50_I(1:360, 1:360);   % simulated time-to-peak (ms)
emp_TTP = delay_peak_matrix_re;         % empirical delay (ms)
Prob    = Prob_sorted;               % empirical probability matrix

% Parameters
prob_thr = 0;                     % probability threshold
time_max = 200;                      % upper bound for valid empirical delay

% Mask: valid pairs only
mask_valid = isfinite(emp_TTP) & emp_TTP > 0 & emp_TTP <= time_max & ...
             N_resp >= 100;
%mask_valid = isfinite(emp_TTP) & emp_TTP > 0 & emp_TTP <= time_max &  N_imp >= 100;
emp_vals = delay_peak_matrix_re(mask_valid);
sim_vals = sim_TTP(mask_valid);
D_cor    = D(1:360,1:360);
d_vals   = D_cor(mask_valid);

% Extract vectors
sim_vec = sim_TTP(mask_valid);
emp_vec = emp_TTP(mask_valid);

[r_thr, p_thr] = corr(sim_vec, emp_vec, 'Type','Pearson');
fprintf('Overall (Prob > %.2f) Spearman ρ = %.3f (p = %.3g, n = %d)\n', ...
        prob_thr, r_thr, p_thr, sum(mask_valid(:)));

[pr_spear, pp_spear] = partialcorr(sim_vec, emp_vec, d_vals, ...
        'Type','Pearson', 'Rows','complete');

fprintf('Partial (Prob > %.2f) Spearman ρ = %.3f (p = %.3g, n = %d)\n', ...
        prob_thr, pr_spear, pp_spear, sum(mask_valid(:)));

% ---------------- Spearman correlation (main scatter) ----------------
valid = isfinite(sim_vec) & isfinite(emp_vec);
[r,p] = corr(sim_vec(valid), emp_vec(valid), 'Type','Spearman');
%p = corr(sim_vec(valid), emp_vec(valid), 'Type','Pearson','Rows','pairwise');
n_eff = sum(valid);

%% -------------------- Figure 1: Raw TTP vs SEEG delay --------------------
bigFS   = 22;      % base font size (≈2x bigger)
smallFS = 20;      % for annotations

figure('Color','w','Position',[100 100 520 450]);

% Scatter: x = empirical, y = simulated
hSc = scatter(emp_vec, sim_vec, 25, Prob(mask_valid), 'filled'); 
colormap(turbo);
c = colorbar; 
c.Label.String = 'Empirical Prob.';

hold on;

% x-axis limits
xlim([0 time_max]);
% y-axis limits: you can uncomment if you want symmetric limits
% ylim([0 time_max]);

%% %%% Linear fit + 95% CI band (no toolboxes needed)
x_data = emp_vec(valid);
y_data = sim_vec(valid);

% Design matrix and OLS fit (y = b0 + b1*x)
X  = [ones(numel(x_data),1), x_data(:)];
b  = X \ y_data(:);

% Grid for smooth line
x_fit  = linspace(min(x_data), max(x_data), 200)';
X_fit  = [ones(numel(x_fit),1), x_fit];
y_fit  = X_fit * b;

% Residual-based variance estimate
resid = y_data(:) - X * b;
dof   = numel(y_data) - 2;
s2    = sum(resid.^2) / dof;

% Covariance of coefficients and SE of predictions
covb = s2 * inv(X.' * X);
se_y = sqrt(diag(X_fit * covb * X_fit.'));

% 95% CI using normal approximation (n is large)
z_crit  = 1.96;
y_lo = y_fit - z_crit * se_y;
y_hi = y_fit + z_crit * se_y;

% CI band (behind points)
hFill = fill([x_fit; flipud(x_fit)], [y_lo; flipud(y_hi)], [0.3 0.3 0.3], ...
             'FaceAlpha',0.15, 'EdgeColor','none');
uistack(hFill,'bottom');

% Regression line
plot(x_fit, y_fit, 'k-', 'LineWidth', 2);

% Bring scatter to front
uistack(hSc,'top');

xlabel('Empirical delay-to-peak (ms)');
ylabel('Simulated time-to-peak (ms)');

ttl = sprintf('TTP vs SEEG delay (\\rho = %.2f, p = %.1e, n = %d)', r, p, n_eff);
title(ttl);

% Corner label
xl = xlim; yl = ylim;
% text(xl(1) + 0.02*range(xl), yl(2) - 0.06*range(yl), ...
%     sprintf('\\rho = %.2f\\n p = %.1e\\n n = %d', r, p, n_eff), ...
%     'FontName','Helvetica','FontSize',smallFS, ...
%     'BackgroundColor','w','Margin',6,'EdgeColor',[0.85 0.85 0.85]);

% ---- Fonts 2x larger ----
set(gca,'FontName','Helvetica','FontSize',bigFS);
set(get(gca,'XLabel'),'FontSize',bigFS);
set(get(gca,'YLabel'),'FontSize',bigFS);
set(get(gca,'Title'), 'FontSize',bigFS+2);
set(c,'FontSize',bigFS-2);
c.Label.FontSize = bigFS;

grid on;

%% -------- Partial Spearman residual–residual scatter (| Distance) --------
Z = d_vals(:);  % confound (distance)

% Keep finite values
ok = isfinite(sim_vec) & isfinite(emp_vec) & isfinite(Z);
x = sim_vec(ok);      
y = emp_vec(ok);      
Z = Z(ok);

% Rank-transform for Spearman partial
xr = tiedrank(x);     
yr = tiedrank(y);     
Zr = tiedrank(Z);

% Regress out Z from both x and y
Xz = [ones(size(Zr,1),1), Zr];
bx = Xz \ xr;   rx = xr - Xz*bx;   % residuals of Sim | Distance
by = Xz \ yr;   ry = yr - Xz*by;   % residuals of Emp | Distance

% Partial Spearman = Pearson between residual ranks
[rp, pp]  = corr(rx, ry, 'type','Pearson');
n_partial = numel(rx);

%% -------------------- Figure 2: Residual–residual scatter --------------------
figure('Color','w','Position',[100 100 520 450]);

% x = empirical residuals, y = simulated residuals
hSc2 = scatter(ry, rx, 18, Z, 'filled'); 
colormap(turbo); 
cb = colorbar; 
cb.Label.String = 'Distance (mm)';

hold on;
grid on;

xlabel('Empirical delay (rank residuals | Distance)');
ylabel('Sim TTP (rank residuals | Distance)');
title(sprintf('Partial Spearman (|Distance): \\rho = %.3f (p = %.2g, n = %d)', ...
      rp, pp, n_partial));

%% %%% Linear fit + 95% CI band on residuals
x_res = ry(:);   % x-axis
y_res = rx(:);   % y-axis

% Valid
valid_res = isfinite(x_res) & isfinite(y_res);
x_res = x_res(valid_res);
y_res = y_res(valid_res);

Xr  = [ones(numel(x_res),1), x_res];
br  = Xr \ y_res;

x_fit_r = linspace(min(x_res), max(x_res), 200)';
X_fit_r = [ones(numel(x_fit_r),1), x_fit_r];
y_fit_r = X_fit_r * br;

resid_r = y_res - Xr * br;
dof_r   = numel(y_res) - 2;
s2_r    = sum(resid_r.^2) / dof_r;

covb_r = s2_r * inv(Xr.' * Xr);
se_y_r = sqrt(diag(X_fit_r * covb_r * X_fit_r.'));
y_lo_r = y_fit_r - z_crit * se_y_r;
y_hi_r = y_fit_r + z_crit * se_y_r;

hFill2 = fill([x_fit_r; flipud(x_fit_r)], ...
              [y_lo_r;  flipud(y_hi_r)], [0.3 0.3 0.3], ...
              'FaceAlpha',0.15, 'EdgeColor','none');
uistack(hFill2,'bottom');

plot(x_fit_r, y_fit_r, 'k-', 'LineWidth', 2);
uistack(hSc2,'top');

% Optional text annotation with stats
xl = xlim; yl = ylim;
% text(xl(1)+0.03*range(xl), yl(2)-0.08*range(yl), ...
%      sprintf('\\rho = %.3f\\n p = %.2g\\n n = %d', rp, pp, n_partial), ...
%      'FontName','Helvetica','FontSize',smallFS, ...
%      'BackgroundColor','w','EdgeColor',[0.8 0.8 0.8],'Margin',6);

% ---- Fonts 2x larger ----
set(gca,'FontName','Helvetica','FontSize',bigFS);
set(get(gca,'XLabel'),'FontSize',bigFS);
set(get(gca,'YLabel'),'FontSize',bigFS);
set(get(gca,'Title'), 'FontSize',bigFS+2);
set(cb,'FontSize',bigFS-2);
cb.Label.FontSize = bigFS;

grid on;
%% -------------------- Figure 2: Partial Pearson (control distance) --------------------
% valid already defined:
% valid = isfinite(sim_vec) & isfinite(emp_vec) & isfinite(d_vals);

x_raw = emp_vec(valid);
y_raw = sim_vec(valid);
d_raw = d_vals(valid);

% Regress out distance from x and y (linear)
Xd = [ones(numel(d_raw),1), d_raw(:)];

bx = Xd \ x_raw(:);
by = Xd \ y_raw(:);

x_res = x_raw(:) - Xd*bx;
y_res = y_raw(:) - Xd*by;

% Partial Pearson equals corr(x_res, y_res)
[r_par, p_par] = corr(x_res, y_res, 'Type','Pearson');
n_par = numel(x_res);

figure('Color','w','Position',[650 100 520 450]);

% Optional: color by Prob (same masking)
hSc = scatter(x_res, y_res, 25, Prob(valid), 'filled');
colormap(turbo);
c = colorbar;
c.Label.String = 'Empirical Prob.';
hold on;

% OLS fit on residuals (y_res = b0 + b1*x_res)
X  = [ones(n_par,1), x_res];
b  = X \ y_res;

x_fit = linspace(min(x_res), max(x_res), 200)';
X_fit = [ones(numel(x_fit),1), x_fit];
y_fit = X_fit * b;

resid = y_res - X*b;
dof   = n_par - 2;
s2    = sum(resid.^2) / dof;

covb = s2 * inv(X.' * X);
se_y = sqrt(diag(X_fit * covb * X_fit.'));

% Prefer t critical (more correct than z)
t_crit = tinv(0.975, dof);
y_lo = y_fit - t_crit * se_y;
y_hi = y_fit + t_crit * se_y;

hFill = fill([x_fit; flipud(x_fit)], [y_lo; flipud(y_hi)], [0.3 0.3 0.3], ...
             'FaceAlpha',0.15, 'EdgeColor','none');
uistack(hFill,'bottom');

plot(x_fit, y_fit, 'k-', 'LineWidth', 2);
uistack(hSc,'top');

xlabel('Empirical delay residual (ms, distance regressed out)');
ylabel('Simulated TTP residual (ms, distance regressed out)');

title(sprintf('Partial Pearson (control distance): r = %.2f, p = %.1e, n = %d', ...
      r_par, p_par, n_par));

grid on;

% ---- Fonts (match your style) ----
bigFS = 22;
set(gca,'FontName','Helvetica','FontSize',bigFS);
set(get(gca,'XLabel'),'FontSize',bigFS);
set(get(gca,'YLabel'),'FontSize',bigFS);
set(get(gca,'Title'), 'FontSize',bigFS+2);
set(c,'FontSize',bigFS-2);
c.Label.FontSize = bigFS;
