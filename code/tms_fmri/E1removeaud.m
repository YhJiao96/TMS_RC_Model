c1 = cifti_read('Schaefer2018_1000Parcels_7Networks_order.dlabel.nii');
c2 = cifti_read('FMRIsurface.dscalar.nii');
mnipositiont = readtable('mnipositionc.csv');
mniposition = table2array(mnipositiont);
SinaSchOrder = table2array(readtable('Schreorder.csv'));
ReOrder = SinaSchOrder(:,1);
%mnipos = mniposition+(abs(min(min(mniposition))));
%mni = mnipos./max(max(mnipos)).*max(max(coor));
%indices = knnsearch(mni(1:64984,:), coor);
%a2 = c2.cdata(indices);
%corr1 = zeros(1,1000);
sti_loc = 268;
%noise_loc = [];
%noise_loc = [240:300];
noise_loc  = [96 101 92 108 93 109 90 99 238 234 253 256 106 95 89 97 111 98 105 86 259 263 102 265 103 87 83 607 601 597 606 754 613 612 600 603 611 589 598 596 590 608 787 776 756 614 609 610 1 5 394 498 499 500,303 304 305 306 308 309 310 311 314 315 316 317];
%noise_loc = [1 5 394 498 499 500,303 304 305 306 308 309 310 311 314 315 316 317];
%noise_loc  = [96 101 92 108 93 109 90 99 238 234 253 256 106 95 89 97 111 98 105 86 259 263 102 265 103 87 83];
[~,noise_loc_ro] = ismember(noise_loc, ReOrder);

% sti_loc_ro = find(ReOrder == sti_loc);
%no = cat(1,sti_loc,noise_loc',find(ReOrder>1000));
 n_out = cat(1,find(ReOrder == sti_loc),noise_loc_ro',find(ReOrder>1000));
% %array_off = setdiff(1:500,n_out);
 array_off = setdiff(1:1000,n_out);
%n_out = find(ReOrder>1000);
ReOrder(n_out) = [];
Surf1  = abs(c2.cdata);
a2 = zeros(1012,1);
vertices = zeros(1012,3);
a2num = zeros(1012,1);
for i = 1:1012
    for j = 1:size(c1.cdata,1)
        if c1.cdata(j) == i
            a2(i) = a2(i)+Surf1(j);
            a2num(i) = a2num(i)+1;
            
            %Atlas_normE(i)=Atlas_normE(i)+normE(j);
        else 
            continue
        end
    end
end
a2aver = a2./a2num;


%a2aver(isnan(a2aver))=0;
a2n = a2aver(array_off);

max_half = max(find(ReOrder<=1000));
a2n = a2n(1:max_half);


load SC_ds_hcp_tract_pro_ctx_sch7n1000p_sctx_fs14_n_1000_m_1000_group_mean_slcthr_5_nf_na.mat;
load ED_ctx_sch7n1000p_sctx_na.mat;
sti_loc = 268;
%SC(SC<=0.001*max(max(SC))) = 0;
%SC = log(SC+1);
%SC = SC(1:1000,1:1000).*(D.^5);
%SC = SC(1:1000,1:1000).*(D.^5);
% Range of n values
n_values = -1:10;  % You can adjust the step size as needed
r_values = zeros(length(n_values), 5);  % To store correlation coefficients
p_values = zeros(length(n_values), 5);  % To store p-values
ReOrder1 = SinaSchOrder(:,1);
sti_loc1 = ReOrder1(sti_loc);
% Loop over n and compute correlations
for idx = 1:length(n_values)
    n = n_values(idx);

    % Update SC with different power of D
    SC_n = SC(1:1000,1:1000) .* (D.^n);
    SC_n(isnan(SC_n)) = 0;

    % Recalculate all metrics based on the new SC
    MFPT = mean_first_passage_time(SC_n);
    diff_eff = 1 ./ MFPT;
    diff_eff(isinf(diff_eff)) = 0;

    SI = search_information(SC_n, 1 ./ SC_n);
    SI(isinf(SI)) = 0;
    SI(isnan(SI)) = 0;
    
    s = diag(sum(SC_n, 1));
    CMY = expm(s^(-1/2) * SC_n * s^(-1/2));

    L = 1 ./ SC_n(1:1000, 1:1000);
    L(L == Inf) = 0;
    [~, ~, PL_wei, ~, ~] = navigation_wu(L, D);
    Nav_eff = 1 ./ PL_wei;
    Nav_eff(isinf(Nav_eff)) = 0;

    [SPL, ~, ~] = distance_wei_floyd(1 ./ SC_n);
    SPE = 1 ./ SPL;
    SPE(isinf(SPE)) = 0;

    % Calculate correlations for each measure
    metrics = {diff_eff, SI, CMY, Nav_eff, SPE};
    for j = 1:5
        a1_metric = metrics{j}(:, sti_loc1) ;
        a1_metric = a1_metric(ReOrder);
        a1_metric = a1_metric(1:max_half);

        [r, p] = corr(a1_metric, a2n, 'Type', 'Spearman');
        r_values(idx, j) = r;
        p_values(idx, j) = p;
    end
end

% Plotting the results
figure;
for j = 1:5
    subplot(2, 3, j);
    plot(n_values, r_values(:, j), '-o');
    hold on;
    yyaxis right;
    plot(n_values, p_values(:, j), '--*');
    xlabel('n');
    ylabel('Correlation (r)');
    yyaxis right;
    ylabel('p-value');
    title(['Metric ', num2str(j)]);
    legend('r', 'p', 'Location', 'Best');
end

sgtitle('Correlation and p-values for different values of n');




% Range of threshold values
threshold_values = 0:0.0005:0.006;  % Adjust the step size if needed
r_values_thresh = zeros(length(threshold_values), 5);  % To store correlation coefficients
p_values_thresh = zeros(length(threshold_values), 5);  % To store p-values

% Loop over threshold values and compute correlations
for idx = 1:length(threshold_values)
    threshold = threshold_values(idx);
    
    % Apply threshold to SC (elements below the threshold are set to 0)
    SC_thresh = SC(1:1000, 1:1000);
    SC_thresh(SC_thresh <= threshold * max(max(SC_thresh))) = 0;
    
    % Recalculate all metrics based on the new SC
    MFPT = mean_first_passage_time(SC_thresh);
    diff_eff = 1 ./ MFPT;
    diff_eff(isinf(diff_eff)) = 0;
    
    SI = search_information(SC_thresh, 1 ./ SC_thresh);
    SI(isinf(SI)) = 0;
    SI(isnan(SI)) = 0;

    s = diag(sum(SC_thresh, 1));
    CMY = expm(s^(-1/2) * SC_thresh * s^(-1/2));
    
    L = 1 ./ SC_thresh(1:1000, 1:1000);
    L(L == Inf) = 0;
    [~, ~, PL_wei, ~, ~] = navigation_wu(L, D);
    Nav_eff = 1 ./ PL_wei;
    Nav_eff(isinf(Nav_eff)) = 0;
    
    [SPL, ~, ~] = distance_wei_floyd(1 ./ SC_thresh);
    SPE = 1 ./ SPL;
    SPE(isinf(SPE)) = 0;
    
    % Calculate correlations for each measure
    metrics = {diff_eff, SI, CMY, Nav_eff, SPE};
    for j = 1:5
        a1_metric = metrics{j}(:, sti_loc) ;
        a1_metric = a1_metric(ReOrder);
        a1_metric = a1_metric(1:max_half);
        
        [r, p] = corr(a1_metric, a2n, 'Type', 'Spearman');
        r_values_thresh(idx, j) = r;
        p_values_thresh(idx, j) = p;
    end
end

% Plotting the results for different thresholds
figure;
for j = 1:5
    subplot(2, 3, j);
    plot(threshold_values, r_values_thresh(:, j), '-o');
    hold on;
    yyaxis right;
    plot(threshold_values, p_values_thresh(:, j), '--*');
    xlabel('Threshold (fraction of max SC)');
    ylabel('Correlation (r)');
    yyaxis right;
    ylabel('p-value');
    title(['Metric ', num2str(j)]);
    legend('r', 'p', 'Location', 'Best');
end

sgtitle('Correlation and p-values for different thresholds on SC');







% Range of percentages for selecting top n% nodes
percentage_values = 10:10:100;  % Adjust the step size if needed (e.g., 10%, 20%, ..., 100%)
r_values_percent = zeros(length(percentage_values), 5);  % To store correlation coefficients
p_values_percent = zeros(length(percentage_values), 5);  % To store p-values
overlap_rates = zeros(length(percentage_values), 5);  % To store spatial overlap rates

% Loop over different percentages and compute correlations and spatial overlap
for idx = 1:length(percentage_values)
    percentage = percentage_values(idx);
    
    % Loop over each metric and perform correlation analysis
    metrics = {diff_eff, SI, CMY, Nav_eff, SPE};
    for j = 1:5
        % Extract the metric and fMRI data
        a1_metric = metrics{j}(:, sti_loc) ;
        a1_metric = a1_metric(ReOrder);
        a1_metric = a1_metric(1:max_half);
        
        a2_metric = a2n;

        % Sort the values to find the top n%
        [~, sorted_indices_a1] = sort(a1_metric, 'descend');
        [~, sorted_indices_a2] = sort(a2_metric, 'descend');
        
        num_top_nodes = round((percentage / 100) * length(a1_metric));
        top_indices_a1 = sorted_indices_a1(1:num_top_nodes);
        top_indices_a2 = sorted_indices_a2(1:num_top_nodes);
        
        % Use only the top n% of nodes for correlation
        a1_top = a1_metric(top_indices_a1);
        a2_top = a2_metric(top_indices_a1);
        
        % Compute correlation for the top n% nodes
        [r, p] = corr(a1_top, a2_top, 'Type', 'Spearman');
        r_values_percent(idx, j) = r;
        p_values_percent(idx, j) = p;
        
        % Compute overlap rate
        overlap_rate = length(intersect(top_indices_a1, top_indices_a2)) / num_top_nodes;
        overlap_rates(idx, j) = overlap_rate;
    end
end

% Plotting the results for different top n% nodes
figure;
for j = 1:5
    subplot(3, 2, j);
    % Plot correlation values
    yyaxis left;
    plot(percentage_values, r_values_percent(:, j), '-o');
    ylabel('Correlation (r)');
    
    % Plot spatial overlap rates
    yyaxis right;
    plot(percentage_values, overlap_rates(:, j), '--*');
    ylabel('Overlap Rate');
    
    xlabel('Top n% Nodes');
    title(['Metric ', num2str(j)]);
    legend('r', 'Overlap Rate', 'Location', 'Best');
end

sgtitle('Correlation and Spatial Overlap Rates for Different Top n% Nodes');


% 
% c = cifti_read('Schaefer2018_1000Parcels_7Networks_order.dlabel.nii');
% tstep = 1;
% clabel0 = c.cdata;
% clabel1  = zeros(size(clabel0,1),tstep);   
%     for i = 1:size(clabel0,1)
%         if clabel0(i)==0
%             continue
%         % elseif clabel0(i)>500
%         %     continue
%         elseif ismember(clabel0(i),n_out)
%             continue
%         elseif ismember(clabel0(i),array_off)
%             clabel1(i) = (clabel0(i)==array_off);
%             %clabel1(i)=clabel0(i);
%         end
%     end
% 
% 
%     c.cdata = clabel1;
%     filename_var = strcat('Atlas_SI_Res',num2str(sti_loc),'.dscalar.nii');
%     ciftisavereset(c,filename_var);
