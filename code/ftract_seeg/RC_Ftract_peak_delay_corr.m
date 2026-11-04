%% -------------------- Scatter Plot: TTP vs Empirical Delay --------------------

% Replace cutoff values with NaN
TTP(TTP >= 200) = NaN;
sim_TTP = DEL.q50_I(1:360, 1:360);                  % simulated time-to-peak (ms)
emp_TTP = delay_peak_matrix;       % empirical delay (ms)
Prob    = Prob_sorted;             % empirical probability matrix

% Parameters
prob_thr = 0.25;                   % probability threshold
time_max = 200;                    % upper bound for valid empirical delay

% Mask: valid pairs only
mask_valid = isfinite(emp_TTP) & emp_TTP > 0 & emp_TTP <= time_max &  N_imp >= 100 & F_cor(1:360,1:360) < 0.25;
emp_vals = delay_peak_matrix(mask_valid);
sim_vals = sim_TTP(mask_valid);
D_cor = D(1:360,1:360);
d_vals = D_cor(mask_valid);
% Extract vectors
sim_vec = sim_TTP(mask_valid);
emp_vec = emp_TTP(mask_valid);
[r_thr, p_thr] = corr(sim_vec, emp_vec, 'Type','Spearman');
fprintf('Overall (Prob > %.2f) Spearman ρ = %.3f (p = %.3g, n = %d)\n', ...
        prob_thr, r_thr, p_thr, sum(mask_valid(:)));
[pr_spear, pp_spear] = partialcorr(sim_vec, emp_vec, d_vals, ...
        'Type','Spearman', 'Rows','complete');
fprintf('Partial (Prob > %.2f) Spearman ρ = %.3f (p = %.3g, n = %d)\n', ...
        prob_thr, pr_spear, pp_spear, sum(mask_valid(:)));
% Compute Spearman correlation
[r, p] = corr(sim_vec, emp_vec, 'Type','Spearman');

% Plot
figure('Color','w','Position',[100 100 520 450]);
scatter(emp_vec, sim_vec, 25, Prob(mask_valid), 'filled'); % color = probability
colormap(turbo);
c = colorbar; c.Label.String = 'Empirical Prob.';
xlabel('Empirical delay-to-peak (ms)');
ylabel('Simulated time-to-peak (ms)');
% Recompute valid mask for correlation (guards against NaNs)
valid = isfinite(sim_vec) & isfinite(emp_vec);
r = corr(sim_vec(valid), emp_vec(valid), 'Type','Spearman');
p = corr(sim_vec(valid), emp_vec(valid), 'Type','Spearman','Rows','pairwise'); % same as above; fine to keep single call

n_eff = sum(valid);
ttl = sprintf('TTP vs SEEG delay (\\rho = %.2f, p = %.1e, n = %d)', r, p, n_eff);
title(ttl);

% Optional corner label
xl = xlim; yl = ylim;
text(xl(1)+0.02*range(xl), yl(2)-0.06*range(yl), ...
    sprintf('\\rho = %.2f\\n p = %.1e\\n n = %d', r, p, n_eff), ...
    'FontName','Helvetica','FontSize',10, ...
    'BackgroundColor','w','Margin',6,'EdgeColor',[0.85 0.85 0.85]);

title(sprintf('TTP vs SEEG delay (ρ = %.2f, p = %.1e, n = %d)', r, p, numel(emp_vec)));

% Visual line y = x (perfect match)
hold on;
xlim([0 time_max]); %ylim([0 time_max]);
%plot([0 time_max],[0 time_max],'k--','LineWidth',1.0);

set(gca,'FontName','Helvetica','FontSize',11);
grid on;

%% -------- Partial Spearman residual–residual scatter (| Distance, color = Distance) --------
Z = d_vals(:);  % confound (distance)

% Keep finite values
ok = isfinite(sim_vec) & isfinite(emp_vec) & isfinite(Z);
x = sim_vec(ok);      y = emp_vec(ok);      Z = Z(ok);

% Rank-transform for Spearman partial
xr = tiedrank(x);     yr = tiedrank(y);     Zr = tiedrank(Z);

% Regress out Z from both x and y
Xz = [ones(size(Zr,1),1) Zr];
bx = Xz \ xr;   rx = xr - Xz*bx;   % residuals of Sim | Distance
by = Xz \ yr;   ry = yr - Xz*by;   % residuals of Emp | Distance

% Partial Spearman = Pearson between residual ranks
[rp, pp] = corr(rx, ry, 'type','Pearson');
n_partial = numel(rx);

% --- Plot residual–residual scatter color-coded by distance ---
figure('Color','w','Position',[100 100 520 450]);
scatter(ry, rx, 18, Z, 'filled'); 
colormap(turbo);
cb = colorbar; 
cb.Label.String = 'Distance (mm)';

hold on; lsline; grid on;
xlabel('Empirical delay (rank residuals | Distance)');
ylabel('Sim TTP (rank residuals | Distance)');
title(sprintf('Partial Spearman (|Distance): \\rho = %.3f (p = %.2g, n = %d)', rp, pp, n_partial));
set(gca,'FontName','Helvetica','FontSize',11);

% Optional text annotation with stats
xl = xlim; yl = ylim;
text(xl(1)+0.03*range(xl), yl(2)-0.08*range(yl), ...
     sprintf('\\rho = %.3f\\n p = %.2g\\n n = %d', rp, pp, n_partial), ...
     'FontName','Helvetica','FontSize',10, ...
     'BackgroundColor','w','EdgeColor',[0.8 0.8 0.8],'Margin',6);




% % (Optional) Also add a neat corner label on the axes
% xl = xlim; yl = ylim;
% tx = xl(1) + 0.02*range(xl);
% ty = yl(2) - 0.06*range(yl);
% text(tx, ty, sprintf('\\rho = %.3f\\n p = %.2g\\n n = %d', rp, pp, n), ...
%     'FontName','Helvetica','FontSize',10, ...
%     'BackgroundColor','w','Margin',6,'EdgeColor',[0.8 0.8 0.8]);

% 
% %% -------------------- Histogram: Simulated vs Empirical Delay --------------------
% emp_TTP = delay_peak_matrix;     % empirical delay matrix (ms)
% %sim_TTP = DEL.q50_I(1:360, 1:360);                   % simulated time-to-peak matrix (ms)
% Prob    = Prob_sorted;           % empirical probability matrix
% 
% % Parameters
% prob_thr = 0.25;    % threshold for empirical connection
% time_max = 200;     % cutoff (ms)
% bin_edges = 0:1:time_max;  % 5 ms bins
% 
% % Mask: valid empirical entries and strong connections
% mask_valid = isfinite(emp_TTP) & emp_TTP > 0 & emp_TTP <= time_max & Prob > prob_thr;
% 
% % Extract vectors
% emp_vals = emp_TTP(mask_valid);
% sim_vals = sim_TTP(mask_valid);
% 
% % Plot
% figure('Color','w','Position',[100 100 600 420]); hold on;
% 
% % Empirical histogram
% histogram(emp_vals, bin_edges, 'Normalization','pdf', ...
%     'FaceColor',[0.3 0.7 0.9], 'EdgeColor','none', 'FaceAlpha',0.6);
% 
% % Simulated histogram
% histogram(sim_vals, bin_edges, 'Normalization','pdf', ...
%     'FaceColor',[0.9 0.4 0.4], 'EdgeColor','none', 'FaceAlpha',0.5);
% 
% xlabel('Delay-to-peak (ms)');
% ylabel('Probability density');
% legend({'Empirical SEEG delay','Simulated TTP'}, 'Location','northeast');
% title(sprintf('Distribution of delay times (Prob > %.2f, n=%d)', prob_thr, numel(emp_vals)));
% set(gca,'FontName','Helvetica','FontSize',11);
% grid on;
% 
% % Optionally add log-scale if dynamic range is large:
% % set(gca,'YScale','log');
% figure('Color','w','Position',[100 100 520 450]); 
% scatter(ry, rx, 18, Z, 'filled'); 
% colormap(turbo); cb = colorbar; cb.Label.String = 'Distance (mm)';
% hold on; lsline; grid on;
% xlabel('Empirical delay (rank residuals | Distance)');
% ylabel('Sim TTP (rank residuals | Distance)');
% title('Partial Spearman (colored by original Distance)');
% set(gca,'FontName','Helvetica','FontSize',11);
% % ---- Partial Spearman stats from residual ranks ----
% valid = isfinite(rx) & isfinite(ry);
% [rp, pp] = corr(rx(valid), ry(valid), 'type','Pearson');   % partial Spearman
% n = sum(valid);
% 
% % Put stats in the title (or subtitle)
% title(sprintf('Partial Spearman (|Distance): \\rho = %.3f (p = %.2g, n = %d)', rp, pp, n));
% 
% % % (Optional) Also add a neat corner label on the axes
% % xl = xlim; yl = ylim;
% % tx = xl(1) + 0.02*range(xl);
% % ty = yl(2) - 0.06*range(yl);
% % text(tx, ty, sprintf('\\rho = %.3f\\n p = %.2g\\n n = %d', rp, pp, n), ...
% %     'FontName','Helvetica','FontSize',10, ...
% %     'BackgroundColor','w','Margin',6,'EdgeColor',[0.8 0.8 0.8]);
