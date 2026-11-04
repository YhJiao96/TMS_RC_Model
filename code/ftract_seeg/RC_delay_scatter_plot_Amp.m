%% ============================================================
% Scatter plot: simulated RC delay vs empirical SEEG delay
% Main analysis:
%   - Valid empirical delay-to-peak pairs
%   - Empirical delay <= 200 ms
%   - N_resp >= 100
%   - No probability threshold
%
% Statistics:
%   - Raw Spearman correlation
%   - Partial Spearman correlation controlling Euclidean distance
%
% Figure:
%   - Panel A: raw scatter
%   - Panel B: distance-controlled residual scatter
% ============================================================

%% ---------------- Prepare matrices ----------------

% Optional: remove boundary/cutoff simulated values
TTP(TTP >= 200) = NaN;

% Simulated delay metric
% rows = target/response regions, columns = stimulation regions
sim_TTP = TTP(1:360, 1:360);

% Empirical F-TRACT delay-to-peak
emp_TTP = delay_peak_matrix_re;

% Euclidean distance matrix
D_cor = D(1:360, 1:360);

% Upper bound for valid empirical delay
time_max = 200;

%% ---------------- Valid-pair mask ----------------
% IMPORTANT:
% This keeps your original analysis threshold.
% There is no probability threshold here.
% Probability can be used elsewhere, but not for inclusion.

mask_valid = isfinite(emp_TTP) & ...
             isfinite(sim_TTP) & ...
             isfinite(D_cor) & ...
             emp_TTP > 0 & ...
             emp_TTP <= time_max &...
             N_imp >= 100;

%% ---------------- Extract vectorised values ----------------

emp_vec = emp_TTP(mask_valid);
sim_vec = sim_TTP(mask_valid);
d_vals  = D_cor(mask_valid);

% Final complete-case mask
valid = isfinite(emp_vec) & isfinite(sim_vec) & isfinite(d_vals);

emp_vec = emp_vec(valid);
sim_vec = sim_vec(valid);
d_vals  = d_vals(valid);

n_eff = numel(emp_vec);

%% ---------------- Raw Spearman correlation ----------------

[r_spear, p_spear] = corr(emp_vec, sim_vec, ...
    'Type', 'Pearson', ...
    'Rows', 'complete');

%% ---------------- Partial Spearman controlling Euclidean distance ----------------

[r_part, p_part] = partialcorr(emp_vec, sim_vec, d_vals, ...
    'Type', 'Spearman', ...
    'Rows', 'complete');

fprintf('\n');
fprintf('Overall valid response pairs (N_resp >= 100): Spearman rho = %.3f (p = %.3g, n = %d)\n', ...
    r_spear, p_spear, n_eff);

fprintf('Distance-controlled valid response pairs (N_resp >= 100): partial Spearman rho = %.3f (p = %.3g, n = %d)\n', ...
    r_part, p_part, n_eff);

%% ============================================================
% Make figure
%% ============================================================

figure('Color','w','Position',[100 100 1180 480]);
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

bigFS     = 16;
pointSize = 18;

%% ============================================================
% Panel A: Raw empirical delay vs simulated delay
%% ============================================================

nexttile;

scatter(emp_vec, sim_vec, pointSize, ...
    'filled', ...
    'MarkerFaceAlpha', 0.35, ...
    'MarkerEdgeAlpha', 0.10);
hold on;

% Linear fit for visual guidance only
x = emp_vec(:);
y = sim_vec(:);

X = [ones(numel(x),1), x];
b = X \ y;

x_fit = linspace(min(x), max(x), 200)';
y_fit = [ones(numel(x_fit),1), x_fit] * b;

plot(x_fit, y_fit, 'k-', 'LineWidth', 2);

% Binned median ± IQR
nBins = 12;
edges = linspace(min(emp_vec), max(emp_vec), nBins + 1);

bin_center = nan(nBins,1);
bin_med    = nan(nBins,1);
bin_q25    = nan(nBins,1);
bin_q75    = nan(nBins,1);

for ii = 1:nBins
    if ii < nBins
        idx_bin = emp_vec >= edges(ii) & emp_vec < edges(ii+1);
    else
        idx_bin = emp_vec >= edges(ii) & emp_vec <= edges(ii+1);
    end

    if sum(idx_bin) >= 10
        bin_center(ii) = mean(edges(ii:ii+1));
        bin_med(ii)    = median(sim_vec(idx_bin), 'omitnan');
        bin_q25(ii)    = prctile(sim_vec(idx_bin), 25);
        bin_q75(ii)    = prctile(sim_vec(idx_bin), 75);
    end
end

ok_bin = isfinite(bin_center) & isfinite(bin_med);

errorbar(bin_center(ok_bin), bin_med(ok_bin), ...
    bin_med(ok_bin) - bin_q25(ok_bin), ...
    bin_q75(ok_bin) - bin_med(ok_bin), ...
    'ko-', ...
    'LineWidth', 1.8, ...
    'MarkerFaceColor', 'w', ...
    'CapSize', 5);

% y = x reference line
plot([0 time_max], [0 time_max], '--', ...
    'Color', [0.5 0.5 0.5], ...
    'LineWidth', 1);

xlabel('Empirical SEEG delay-to-peak (ms)');
ylabel('Simulated RC delay, q50 current time (ms)');

title(sprintf('Raw association: \\rho = %.3f, p = %.1e, n = %d', ...
    r_spear, p_spear, n_eff));

xlim([0 time_max]);
ylim([0 time_max]);

grid on;
box off;

set(gca, 'FontName', 'Helvetica', 'FontSize', bigFS);

legend({'Region pairs', 'Linear fit', 'Binned median ± IQR', 'y = x'}, ...
    'Location', 'best', ...
    'Box', 'off');

%% ============================================================
% Panel B: Distance-controlled partial Spearman residual plot
%% ============================================================

nexttile;

% For a partial Spearman residual plot:
% 1. Rank-transform empirical delay, simulated delay, and Euclidean distance.
% 2. Regress ranked distance out of ranked empirical delay.
% 3. Regress ranked distance out of ranked simulated delay.
% 4. Correlate the two residuals.
%
% This visualises the same logic as:
% partialcorr(emp_vec, sim_vec, d_vals, 'Type','Spearman').

emp_rank = (emp_vec);
sim_rank = (sim_vec);
d_rank   = (d_vals);

Xd = [ones(numel(d_rank),1), d_rank(:)];

b_emp = Xd \ emp_rank(:);
b_sim = Xd \ sim_rank(:);

emp_res = emp_rank(:) - Xd * b_emp;
sim_res = sim_rank(:) - Xd * b_sim;

% Pearson correlation between rank residuals.
% This should match the partial Spearman result approximately.
[r_res, p_res] = corr(emp_res, sim_res, ...
    'Type', 'Pearson', ...
    'Rows', 'complete');

scatter(emp_res, sim_res, pointSize, ...
    'filled', ...
    'MarkerFaceAlpha', 0.35, ...
    'MarkerEdgeAlpha', 0.10);
hold on;

% Linear fit in residual space
Xr = [ones(numel(emp_res),1), emp_res(:)];
br = Xr \ sim_res(:);

x_fit_r = linspace(min(emp_res), max(emp_res), 200)';
y_fit_r = [ones(numel(x_fit_r),1), x_fit_r] * br;

plot(x_fit_r, y_fit_r, 'k-', 'LineWidth', 2);

% Zero reference lines
xline(0, '--', ...
    'Color', [0.5 0.5 0.5], ...
    'LineWidth', 1);

yline(0, '--', ...
    'Color', [0.5 0.5 0.5], ...
    'LineWidth', 1);

xlabel('Empirical delay residuals after ED control');
ylabel('Simulated RC delay residuals after ED control');

title(sprintf('Distance-controlled partial Spearman: \\rho = %.3f, p = %.2g, n = %d', ...
    r_part, p_part, n_eff));

grid on;
box off;

set(gca, 'FontName', 'Helvetica', 'FontSize', bigFS);

%% ---------------- Optional: save figure ----------------

% exportgraphics(gcf, 'RC_TTP_vs_empirical_delay_scatter.png', 'Resolution', 300);
% exportgraphics(gcf, 'RC_TTP_vs_empirical_delay_scatter.pdf', 'ContentType', 'vector');