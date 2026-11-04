

function [x_clean, y_clean] = remove_outliers_IQR(x, y)
    % Remove non-finite
    valid = isfinite(x) & isfinite(y);
    x = x(valid); 
    y = y(valid);

    % IQR thresholds
    Q1  = quantile(x, 0.25);
    Q3  = quantile(x, 0.75);
    IQRv = Q3 - Q1;

    lower = Q1 - 1.5 * IQRv;
    upper = Q3 + 1.5 * IQRv;

    keep = (x >= lower) & (x <= upper);

    x_clean = x(keep);
    y_clean = y(keep);
end

function [x_clean, y_clean] = remove_outliers_residual(x, y, k)
    % Remove points that are far from the fitted regression line.
    % k = threshold for residuals (robust SD units)

    if nargin < 3
        k = 3; % default to 3 SD
    end

    % Finite points only
    valid = isfinite(x) & isfinite(y);
    x = x(valid);
    y = y(valid);

    if numel(x) < 5
        x_clean = x;
        y_clean = y;
        return;
    end

    % Fit line
    p = polyfit(x, y, 1);
    yhat = polyval(p, x);
    r = y - yhat;   % residuals

    % Robust SD of residuals (MAD → SD)
    sigma_r = 1.4826 * mad(r, 1);

    if sigma_r == 0
        x_clean = x;
        y_clean = y;
        return;
    end

    keep = abs(r) <= k*sigma_r;

    x_clean = x(keep);
    y_clean = y(keep);
end


N_thres = N_resp;
emp_mat = Amp_sorted;
sim_mat = I_node_max_cort;
sim_log_mat = log10(sim_mat + eps);

stim_list = [263, 188];
% 
% figure;
% for k = 1:numel(stim_list)
%     s = stim_list(k);
% 
%     emp_row = emp_mat(s, :)';
%     sim_row = sim_log_mat(s, :)';
% 
%     % ----- IQR OUTLIER REMOVAL -----
%     [sim_row_clean, emp_row_clean] = remove_outliers_IQR(sim_row, emp_row);
% 
%     subplot(1, numel(stim_list), k);
%     scatter(emp_row_clean, sim_row_clean, 25, 'filled'); hold on;
% 
%     % % ---- Linear fit ----
%     % pfit = polyfit(emp_row_clean, sim_row_clean, 1);
%     % xfit = linspace(min(emp_row_clean), max(emp_row_clean), 100);
%     % yfit = polyval(pfit, xfit);
%     % plot(xfit, yfit, 'r-', 'LineWidth', 1.5);
% 
% % ---- Linear fit with CI band ----
% [pfit, S] = polyfit(emp_row_clean, sim_row_clean, 1);   % <-- get S as well
% xfit = linspace(min(emp_row_clean), max(emp_row_clean), 100);
% 
% % 95% confidence band for the regression line
% [yfit, delta] = polyconf(pfit, xfit, S, 'alpha', 0.05); % 'curve' is default
% 
% % Shaded confidence band
% fill([xfit fliplr(xfit)], ...
%      [yfit+delta fliplr(yfit-delta)], ...
%      'r', 'FaceAlpha', 0.15, 'EdgeColor', 'none'); 
% hold on;
% 
% % Regression line
% plot(xfit, yfit, 'r-', 'LineWidth', 1.5);
% 
% 
%     % ---- Spearman correlation ----
%     [rho_row, p_row] = corr(emp_row_clean, sim_row_clean, 'Type','Spearman');
% 
%     xlabel('Empirical Amplitude');
%     ylabel('log_{10}(RC I_{node,max})');
%     title(sprintf('Stim %d: \\rho=%.2f, p=%.1e', s, rho_row, p_row));
% end
% 
% sgtitle('Row-wise scatter with log-scale + IQR-cleaning');
stim_list = [263, 188];

%emp_mat = Prob_sorted;
sim_mat = I_node_max_cort;
sim_log = log10(sim_mat + eps);

for k = 1:numel(stim_list)

    s = stim_list(k);

    % --------------------------------------------
    % Extract row: empirical + simulated
    % --------------------------------------------
    emp_row = emp_mat(s, :)';
    sim_raw = sim_mat(s, :)';
    sim_row = sim_log(s, :)';

    % --------------------------------------------
    % ROW MASK → remove zero + NaN + Inf
    % --------------------------------------------
    row_mask = isfinite(emp_row) & isfinite(sim_row) & ...
               (emp_row ~= 0) & (sim_raw ~= 0);

    emp = emp_row(row_mask);
    sim = sim_row(row_mask);

    % ---- IQR cleaning (1D on x) ----
[sim_iqr, emp_iqr] = remove_outliers_IQR(sim, emp);

% ---- 2D residual cleaning (removes isolated points) ----
[emp_clean, sim_clean] = remove_outliers_residual(emp_iqr, sim_iqr, 3);
% ----- IQR OUTLIER REMOVAL -----
% [sim_clean, emp_clean] = remove_outliers_IQR(sim_vec, emp_vec);
% 

    % --------------------------------------------
    % CREATE NEW FIGURE FOR THIS STIM ONLY
    % --------------------------------------------
    figure;
    set(gcf, 'Position', [200 200 800 700]);  % nice size

    scatter(emp_clean, sim_clean, 40, 'filled'); hold on;

    % --------------------------------------------
    % Linear fit + CI band
    % --------------------------------------------
    [p, S] = polyfit(emp_clean, sim_clean, 1);
    xfit = linspace(min(emp_clean), max(emp_clean), 150);
    [yfit, delta] = polyconf(p, xfit, S, 'alpha', 0.05);

    % CI band
    fill([xfit fliplr(xfit)], ...
         [yfit+delta fliplr(yfit-delta)], ...
         'r', 'FaceAlpha', 0.18, 'EdgeColor', 'none');

    % Fit line
    plot(xfit, yfit, 'r-', 'LineWidth', 2.2);

  

    % --------------------------------------------
    % Spearman correlation
    % --------------------------------------------
    [rho, pval] = corr(emp_clean, sim_clean, 'Type', 'Spearman');

    % --------------------------------------------
    % FORMATTING — 2× fonts everywhere
    % --------------------------------------------
    xlabel('Empirical Amplitude', 'FontSize', 28, 'FontWeight', 'Bold');
    ylabel('log_{10}(RC Current)', 'FontSize', 28, 'FontWeight', 'Bold');

    title(sprintf('Stim %d   \\rho = %.2f,   p = %.1e', s, rho, pval), ...
          'FontSize', 32, 'FontWeight', 'Bold');

    set(gca, 'FontSize', 26, 'LineWidth', 2);
    grid on;

end


%%
% Build your mask first
mask = (N_thres > T);
%mask = mask & blkdiag(ones(180), ones(180));
mask = mask & (D_cor > 20);
%mask = mask & (F_cor < 0.25);
mask = mask & isfinite(Amp_sorted);
mask = mask & isfinite(I_node_max_cort);
% mask = mask & (Prob_sorted ~= 0);
% mask = mask & (I_node_max_cort ~= 0);

emp_vec = Amp_sorted(mask);
sim_vec = log10(I_node_max_cort(mask) + eps);

% ---- IQR cleaning (1D on x) ----
[sim_iqr, emp_iqr] = remove_outliers_IQR(sim_vec, emp_vec);

% ---- 2D residual cleaning (removes isolated points) ----
[emp_clean, sim_clean] = remove_outliers_residual(emp_iqr, sim_iqr, 3);


% ---- Spearman -----
[rho_all, p_all] = corr(sim_clean, emp_clean, 'Type','Spearman');
% 
% figure;
% scatter(emp_clean, sim_clean, 10, 'filled'); hold on;
% 
% % % Linear fit
% % pfit = polyfit(emp_clean, sim_clean, 1);
% % xfit = linspace(min(emp_clean), max(emp_clean), 200);
% % yfit = polyval(pfit, xfit);
% % plot(xfit, yfit, 'r-', 'LineWidth', 1.6);
% figure;
% scatter(emp_clean, sim_clean, 10, 'filled'); hold on;
% 
% % Linear fit with CI band
% [pfit_all, S_all] = polyfit(emp_clean, sim_clean, 1);
% xfit = linspace(min(emp_clean), max(emp_clean), 200);
% [yfit, delta] = polyconf(pfit_all, xfit, S_all, 'alpha', 0.05);
% 
% % CI band
% fill([xfit fliplr(xfit)], ...
%      [yfit+delta fliplr(yfit-delta)], ...
%      'r', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
% plot(xfit, yfit, 'r-', 'LineWidth', 1.6);
% 
% xlabel('Empirical Amplitude');
% ylabel('log_{10}(RC I_{node,max})');
% title(sprintf('Masked Whole-Matrix: \\rho=%.2f, p=%.1e', rho_all, p_all));
% legend({'Data','95% CI of fit','Linear fit'}, 'Location','best');
% grid on;
% 
% 
% % % Lowess
% % yy = smooth(emp_clean, sim_clean, 0.15, 'lowess');
% % [x_sorted, idx_sort] = sort(emp_clean);
% % plot(x_sorted, yy(idx_sort), 'g--', 'LineWidth', 1.6);
% 
% xlabel('Empirical Amplitude');
% ylabel('log_{10}(RC I_{node,max})');
% title(sprintf('Masked Whole-Matrix: \\rho=%.2f, p=%.1e', rho_all, p_all));
% legend({'Data','Linear fit','Lowess'});
% grid on;
figure;
set(gcf,'Position',[200 200 800 700]);

scatter(emp_clean, sim_clean, 25, 'filled'); hold on;

% Linear fit + CI band
[pfit_all, S_all] = polyfit(emp_clean, sim_clean, 1);
xfit = linspace(min(emp_clean), max(emp_clean), 200);
[yfit, delta] = polyconf(pfit_all, xfit, S_all, 'alpha', 0.05);

fill([xfit fliplr(xfit)], ...
     [yfit+delta fliplr(yfit-delta)], ...
     'r', 'FaceAlpha', 0.18, 'EdgeColor','none');

plot(xfit, yfit, 'r-', 'LineWidth', 2.5);

xlabel('Empirical Amplitude', 'FontSize', 28, 'FontWeight','Bold');
ylabel('log_{10}(RC Current)', 'FontSize', 28, 'FontWeight','Bold');

title(sprintf('Whole-Matrix: \\rho = %.2f, p = %.1e', rho_all, p_all), ...
      'FontSize', 32, 'FontWeight','Bold');

legend({'Data','95% CI band','Linear fit'}, ...
       'FontSize', 24, 'Location', 'best');

set(gca, 'FontSize', 26, 'LineWidth', 2);
grid on;
