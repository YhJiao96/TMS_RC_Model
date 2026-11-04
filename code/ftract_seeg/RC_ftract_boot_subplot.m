%% =================== Pretty Full vs Partial Barplots (Prob & Amp) ===================
% Prereqs in workspace (360x360 unless noted):
%   Amp_sorted, Prob_sorted, D_sorted, F_sorted, N_sti, N_resp, N_imp
%   I_node_max, diff_eff, SI, CMY, Nav_eff, SPE, SC
% ------------------------------------------------------------------------------
%
%% Only use when using v2311 F-tract(given idx_amp for reorder);
load('idx_sc.mat');
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

%N_final_matrix = Nimp;amp_final_matrix = Amp;prob_final_matrix = prob;N_sti_final_matrix = Nsti;
[~, idx_amp_full] = sort(idx_sc);     % inverse permutation for all 374
idx_amp = idx_amp_full(1:360);        % only keep first 360 cortical nodes
Amp_sorted = amp_final_matrix(idx_amp,idx_amp) ;
Prob_sorted = prob_final_matrix(idx_amp,idx_amp) ;
N_sti = N_sti_final_matrix(idx_amp,idx_amp);
N_imp = N_final_matrix(idx_amp,idx_amp);
N_resp = N_sti.*Prob_sorted;
% ---------------- Fixed config ----------------
T_pairs   = 100;             % N threshold
D_min     = 20;              % minimum Euclidean distance (mm)
same_hemi = blkdiag(ones(180), ones(180))>0;
B         = 1000;            % bootstrap iters
rng(17);                     % reproducibility

% ---------- Print style (A4-ready) ----------
STYLE = struct( ...
  'figPixels',   [100 100 2400 600], ...   % big canvas (A4 landscape @ ~300dpi)
  'titleFS',     20, ...
  'suptitleFS',  22, ...
  'axisFS',      17, ...
  'labelFS',     18, ...
  'legendFS',    16, ...
  'tickAngle',   20, ...
  'barWidth',    0.9, ...                  % thinner bars than default 0.92
  'errDX',       0.08, ...                   % errorbar horizontal offset
  'catGap',      0.10 ...      % <-- category-to-category gap fraction (0 = bars touch)
);


% Fixed color map (stable across figures)
C = containers.Map;
C('I_node_max') = [0.8392 0.1529 0.1569];  % red  (changed)
C('Diff_eff')   = [0.1725 0.6275 0.1725];  % green
C('SI')         = [0.1216 0.4667 0.7059];  % blue (changed)
C('CMY')        = [0.5804 0.4039 0.7412];  % purple
C('Nav_eff')    = [1.0000 0.4980 0.0549];  % orange
C('SPE')        = [0.0902 0.7451 0.8118];  % teal
C('ED')         = [0.9294 0.6941 0.1255];  % yellow (changed)
C('SC')         = [0.5490 0.3373 0.2941];  % brown



% Metrics list
METRICS = {
    'I_node_max', I_node_max(1:360,1:360)
    'Diff_eff',   diff_eff(1:360,1:360)
    'SI',         SI(1:360,1:360)
    'CMY',        CMY(1:360,1:360)
    'Nav_eff',    Nav_eff(1:360,1:360)
    'SPE',        SPE(1:360,1:360)
    'ED',         D_sorted(1:360,1:360)
    'SC',         SC(1:360,1:360)
};

D_cor = D_sorted(1:360,1:360);
F_cor = F_sorted(1:360,1:360);

% ---------------- Mask builders (anonymous; capture workspace vars) ----------------
build_base_mask = @(N_thres) ...
    (N_thres >= T_pairs) & (N_imp >= 5) & ~eye(360) & same_hemi & (D_cor >= D_min);

mask_all         = @(N_thres) build_base_mask(N_thres);
mask_connected   = @(N_thres) build_base_mask(N_thres) & (F_cor >= 0.75);
mask_unconnected = @(N_thres) build_base_mask(N_thres) & (F_cor <= 0.25);

% ❗ pass anonymous function handles directly (no @ here)
maskSet = {mask_all, mask_connected, mask_unconnected};

%% =================== Render: Prob (3 subplots in one figure) ===================
N_thres_prob = N_sti;
emp_prob     = Prob_sorted;
plot_triptych('Prob: Full vs Partial correlation (95% bootstrap CI)', ...
              emp_prob, N_thres_prob, maskSet, METRICS, D_cor, C, B, STYLE);


%% =================== Render: Amp (3 subplots in one figure) ===================
N_thres_amp = N_resp;
emp_amp     = Amp_sorted;
plot_triptych('Amp: Full vs Partial correlation (95% bootstrap CI)', ...
              emp_amp,  N_thres_amp,  maskSet, METRICS, D_cor, C, B, STYLE);
%% =================== Local functions ===================
% 
% function plot_triptych(fig_title, emp_mat, N_thres, maskFns, METRICS, D_cor, C, B)
%     titles = {'All pairs','Connected pairs','Unconnected pairs'};
%     figure('Color','w','Position',[100 100 1350 360]); 
%     tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
%     for k = 1:3
%         nexttile; ax = gca;
%         M = maskFns{k}(N_thres) & isfinite(emp_mat);
%         [names, counts, rF, rP, ciF, ciP] = compute_stats(emp_mat, METRICS, D_cor, M, B);
%         pretty_dualbars(ax, names, rF, rP, ciF, ciP, counts, C);
%         title(sprintf('%s (n = %d)', titles{k}, nnz(M)), 'FontWeight','bold');
%     end
%     title(tl, fig_title, 'FontSize', 12, 'FontWeight','bold');
% end
function plot_triptych(fig_title, emp_mat, N_thres, maskFns, METRICS, D_cor, C, B, STYLE)
    titles = {'All pairs','Connected pairs','Unconnected pairs'};

    % Big canvas for print
    figure('Color','w','Position',STYLE.figPixels);
    tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

    for k = 1:3
        nexttile; ax = gca;
        M = maskFns{k}(N_thres) & isfinite(emp_mat);

        [names, counts, rF, rP, ciF, ciP] = compute_stats(emp_mat, METRICS, D_cor, M, B);

        pretty_dualbars(ax, names, rF, rP, ciF, ciP, counts, C, STYLE);

        title(sprintf('%s (n = %d)', titles{k}, nnz(M)), ...
            'FontWeight','bold','FontSize',STYLE.titleFS);
        ax.FontSize  = STYLE.axisFS;
        ax.LineWidth = 1.2;
        ylabel('|Spearman r|', 'FontSize', STYLE.labelFS);
    end

    title(tl, fig_title, 'FontSize', STYLE.suptitleFS, 'FontWeight','bold');
end

function [names, counts, r_full, r_part, ci_full, ci_part] = compute_stats(emp_mat, METRICS, D_cor, mask_all, B)
    m = size(METRICS,1);
    names   = strings(m,1);
    r_full  = nan(m,1);  r_part = nan(m,1);
    ci_full = nan(m,2);  ci_part= nan(m,2);
    counts  = nan(m,1);
    for i = 1:m
        names(i) = string(METRICS{i,1});
        sim = METRICS{i,2};
        M   = mask_all & isfinite(sim);
        emp = emp_mat(M); simv = sim(M); Dv = D_cor(M);
        counts(i) = numel(emp);
        if counts(i) < 3, continue; end
        oF = boot_corr(emp, simv, [],  'Partial', false, 'B', B);
        oP = boot_corr(emp, simv, Dv, 'Partial', true,  'B', B);
        r_full(i)  = oF.r_hat;   ci_full(i,:) = oF.ci;
        r_part(i)  = oP.r_hat;   ci_part(i,:) = oP.ci;
    end
    % sort by |r_full|
    [~, idx] = sort(abs(r_full), 'descend');
    names = names(idx); counts = counts(idx);
    r_full = r_full(idx); r_part = r_part(idx);
    ci_full = ci_full(idx,:); ci_part = ci_part(idx,:);
end

% function pretty_dualbars(ax, names, rF, rP, ciF, ciP, counts, Cmap)
%     % --- bar heights (absolute r) ---
%     y1 = abs(rF); y2 = abs(rP);
% 
%     % --- build categories & display labels (fix underscores + special rename) ---
%     raw = cellstr(names);                      % keep keys for color lookup
%     dispLabels = raw;
%     for j = 1:numel(dispLabels)
%         if strcmp(raw{j}, 'I_node_max')
%             dispLabels{j} = 'Inodemax';        % requested label
%         end
%         % otherwise show literal text (no TeX parsing)
%         % leave underscores visible for others (e.g., Diff_eff)
%     end
%     cats = categorical(dispLabels); cats = reordercats(cats, dispLabels);
% 
%     % --- bar colors per metric key ---
%     col = zeros(numel(raw),3);
%     for j = 1:numel(raw)
%         key = raw{j};
%         if isKey(Cmap, key), col(j,:) = Cmap(key); else, col(j,:) = [0.5 0.5 0.5]; end
%     end
% 
%     % --- draw bars ---
%     hold(ax,'on');
%     b1 = bar(cats, y1, 0.92, 'FaceColor','flat','EdgeColor','none','FaceAlpha',0.30); b1.CData = col;
%     b2 = bar(cats, y2, 0.92, 'FaceColor','flat','EdgeColor','none','FaceAlpha',0.95); b2.CData = col; uistack(b2,'top');
% 
%     % --- error bars (absolute scale) ---
%     loF = abs(rF) - abs(ciF(:,1)); hiF = abs(ciF(:,2)) - abs(rF);
%     loP = abs(rP) - abs(ciP(:,1)); hiP = abs(ciP(:,2)) - abs(rP);
%     x1 = b1.XEndPoints; x2 = b2.XEndPoints; dx = 0.12;
%     errorbar(x1-dx, y1, loF, hiF, 'k','LineStyle','none','LineWidth',1);
%     errorbar(x2+dx, y2, loP, hiP, 'k','LineStyle','none','LineWidth',1);
% 
%     % --- cosmetics ---
%     ylabel('|Spearman r|');
%     grid on; ax.GridAlpha = 0.2; ax.LineWidth = 1; ax.FontSize = 10;
%     xtickangle(25);
%     set(ax, 'TickLabelInterpreter','none');     % <-- prevents TeX underscore parsing
%     % --- automatic y-axis scaling per subplot ---
%     ymax = max([y1; y2]);
%     ymin = min([y1; y2]);
%     range = ymax - ymin;
%     if range == 0, range = 0.1; end          % avoid zero range
%     ylim([max(0, ymin - 0.05*range), ymax + 0.15*range]);
% 
%     box on;
% 
%     % --- consistent legend swatches (Full = light red, Partial = dark red) ---
%     FULL_SW  = [1.0 0.6 0.6];
%     PART_SW  = [0.8 0.0 0.0];
%     hFull = patch(nan, nan, FULL_SW, 'EdgeColor','none','FaceAlpha',0.30);
%     hPart = patch(nan, nan, PART_SW, 'EdgeColor','none','FaceAlpha',0.95);
%     lg = legend([hFull, hPart], {'Full (translucent)', 'Partial (opaque)'}, ...
%                 'Location','northeast');           % top-right
%     lg.Box = 'off';
% end
function pretty_dualbars(ax, names, rF, rP, ciF, ciP, counts, Cmap, STYLE)
    % heights
    y1 = abs(rF); y2 = abs(rP);
    n  = numel(y1);

    % labels (keep underscores literal; rename I_node_max)
    raw = cellstr(names);
    dispLabels = raw;
    for j = 1:n
        if strcmp(raw{j}, 'I_node_max'), dispLabels{j} = 'Inodemax'; end
    end

    % colors
    col = zeros(n,3);
    for j = 1:n
        if isKey(Cmap, raw{j}), col(j,:) = Cmap(raw{j}); else, col(j,:) = [0.5 0.5 0.5]; end
    end

    % numeric x so we can control spacing between categories
    x = 1:n;
    width = max(0.05, 1 - STYLE.catGap);    % 0<catGap<1 → smaller gap
    hold(ax,'on');
    b1 = bar(ax, x, y1, width, 'FaceColor','flat','EdgeColor','none','FaceAlpha',0.30); b1.CData = col;
    b2 = bar(ax, x, y2, width, 'FaceColor','flat','EdgeColor','none','FaceAlpha',0.95); b2.CData = col; 
    uistack(b2,'top');

    % error bars (absolute scale)
    loF = abs(rF) - abs(ciF(:,1)); hiF = abs(ciF(:,2)) - abs(rF);
    loP = abs(rP) - abs(ciP(:,1)); hiP = abs(ciP(:,2)) - abs(rP);
    x1 = b1.XEndPoints; x2 = b2.XEndPoints; dx = STYLE.errDX;
    errorbar(x1-dx, y1, loF, hiF, 'k','LineStyle','none','LineWidth',1);
    errorbar(x2+dx, y2, loP, hiP, 'k','LineStyle','none','LineWidth',1);

    % cosmetics
    set(ax, 'TickLabelInterpreter','none');
    xticks(x); xticklabels(dispLabels); xtickangle(STYLE.tickAngle);
    grid on; ax.GridAlpha = 0.2; box on;
    xlim([0.5, n+0.5]);

    ymax = max([y1; y2]); ymin = min([y1; y2]); range = max(0.1, ymax - ymin);
    ylim([max(0, ymin - 0.05*range), ymax + 0.18*range]);

    % legend swatches
    FULL_SW  = [1.0 0.6 0.6];
    PART_SW  = [0.8 0.0 0.0];
    hFull = patch(nan, nan, FULL_SW, 'EdgeColor','none','FaceAlpha',0.30);
    hPart = patch(nan, nan, PART_SW, 'EdgeColor','none','FaceAlpha',0.95);
    lg = legend([hFull, hPart], {'Full (translucent)', 'Partial (opaque)'}, 'Location','northeast');
    lg.Box = 'off'; lg.FontSize = STYLE.legendFS;
end

function out = boot_corr(emp, sim, dist, varargin)
    % out = boot_corr(emp, sim, dist, 'Partial',true/false, 'B',1000)
    p = inputParser;
    addParameter(p, 'Partial', false, @(x)islogical(x)&&isscalar(x));
    addParameter(p, 'B', 1000, @(x)isnumeric(x)&&isscalar(x)&&x>10);
    parse(p, varargin{:});
    isPartial = p.Results.Partial; 
    B = p.Results.B;

    emp = emp(:); sim = sim(:);
    if isPartial, dist = dist(:); end

    if isPartial
        valid = isfinite(emp) & isfinite(sim) & isfinite(dist);
    else
        valid = isfinite(emp) & isfinite(sim);
    end
    emp = emp(valid); sim = sim(valid);
    if isPartial, dist = dist(valid); end

    n = numel(emp);
    if n < 3
        out = struct('r_hat',NaN,'ci',[NaN NaN],'se',NaN,'p_boot',NaN,'r_boot',[]);
        return;
    end

    if isPartial
        r_hat = partialcorr(emp, sim, dist, 'type','Spearman');
    else
        r_hat = corr(emp, sim, 'type','Spearman');
    end

    idxMat = randi(n, [n, B]);
    r_boot = nan(B,1);
    for b = 1:B
        I = idxMat(:,b);
        if isPartial
            r_boot(b) = partialcorr(emp(I), sim(I), dist(I), 'type','Spearman');
        else
            r_boot(b) = corr(emp(I), sim(I), 'type','Spearman');
        end
    end

    out.ci = quantile(r_boot, [0.025 0.975]);
    out.se = std(r_boot, 1);
    out.p_boot = mean(abs(r_boot) >= abs(r_hat));
    out.r_hat = r_hat; out.r_boot = r_boot;
end