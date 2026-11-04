T_prob_results = make_barplot_result_table( ...
    'Probability', ...
    Prob_sorted, ...
    N_sti, ...
    maskSet, ...
    METRICS, ...
    D_cor, ...
    B);

fprintf('\n================ NUMERICAL RESULTS: PROBABILITY ================\n');
disp(T_prob_results);

T_amp_results = make_barplot_result_table( ...
    'Amplitude', ...
    Amp_sorted, ...
    N_resp, ...
    maskSet, ...
    METRICS, ...
    D_cor, ...
    B);

fprintf('\n================ NUMERICAL RESULTS: AMPLITUDE ================\n');
disp(T_amp_results);

% Optional: save tables
writetable(T_prob_results, 'RC_probability_barplot_results.csv');
writetable(T_amp_results,  'RC_amplitude_barplot_results.csv');

function T_all = make_barplot_result_table(analysis_name, emp_mat, N_thres, maskFns, METRICS, D_cor, B)

    pair_titles = {'All pairs', 'Connected pairs', 'Unconnected pairs'};

    T_all = table();

    for k = 1:3

        M = maskFns{k}(N_thres) & isfinite(emp_mat);

        [names, counts, rF, rP, ciF, ciP] = compute_stats_for_table( ...
            emp_mat, METRICS, D_cor, M, B);

        T_k = table();

        T_k.Analysis = repmat(string(analysis_name), numel(names), 1);
        T_k.PairSet  = repmat(string(pair_titles{k}), numel(names), 1);
        T_k.Model    = string(names);
        T_k.NPairs   = counts;

        T_k.Full_rho      = rF;
        T_k.Full_CI_low   = ciF(:,1);
        T_k.Full_CI_high  = ciF(:,2);

        T_k.EDpartial_rho     = rP;
        T_k.EDpartial_CI_low  = ciP(:,1);
        T_k.EDpartial_CI_high = ciP(:,2);

        T_all = [T_all; T_k];
    end
end


function [names, counts, r_full, r_part, ci_full, ci_part] = compute_stats_for_table(emp_mat, METRICS, D_cor, mask_all, B)

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

        oF = boot_corr_for_table(emp, simv, [], ...
            'Partial', false, ...
            'B', B);

        oP = boot_corr_for_table(emp, simv, Dv, ...
            'Partial', true, ...
            'B', B);

        r_full(i) = oF.r_hat;
        r_part(i) = oP.r_hat;

        ci_full(i,:) = oF.ci;
        ci_part(i,:) = oP.ci;
    end

    [~, idx] = sort(abs(r_full), 'descend', 'MissingPlacement','last');

    names   = names(idx);
    counts  = counts(idx);
    r_full  = r_full(idx);
    r_part  = r_part(idx);
    ci_full = ci_full(idx,:);
    ci_part = ci_part(idx,:);
end


function out = boot_corr_for_table(emp, sim, dist, varargin)

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