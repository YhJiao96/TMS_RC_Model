% Inputs:
% I_node_max: 360x360
% Amp_sorted: 360x360
% Prob_sorted: 360x360
I_node_max_cort = I_node_max(1:360,1:360);
n = size(I_node_max_cort(1:360,1:360),1);

rho_amp  = zeros(n,1);
p_amp    = zeros(n,1);
rho_prob = zeros(n,1);
p_prob   = zeros(n,1);

for i = 1:n
    [rho_amp(i),  p_amp(i)]  = corr(I_node_max_cort(i,:)', Amp_sorted(i,:)', 'Type','Spearman', 'Rows','pairwise');
    [rho_prob(i), p_prob(i)] = corr(I_node_max_cort(i,:)', Prob_sorted(i,:)', 'Type','Spearman', 'Rows','pairwise');
end

rho_amp_col  = zeros(n,1);
p_amp_col    = zeros(n,1);
rho_prob_col = zeros(n,1);
p_prob_col   = zeros(n,1);

for j = 1:n
    [rho_amp_col(j),  p_amp_col(j)]  = corr(I_node_max_cort(:,j), Amp_sorted(:,j), 'Type','Spearman', 'Rows','pairwise');
    [rho_prob_col(j), p_prob_col(j)] = corr(I_node_max_cort(:,j), Prob_sorted(:,j), 'Type','Spearman', 'Rows','pairwise');
end

figure;
subplot(2,2,1)
bar(rho_amp); ylim([0 max(rho_amp)*1.1])
title('Row-wise Spearman ρ: RC → SEEG Amplitude');
xlabel('Stimulation Site'); ylabel('\rho');
subplot(2,2,2)
bar(p_amp); ylim([0 max(p_amp)*1.1])
title('Row-wise Spearman p: RC → SEEG Amplitude');
xlabel('Stimulation Site'); ylabel('\p');
subplot(2,2,3)
bar(rho_prob); ylim([0 max(rho_prob)*1.1])
title('Row-wise Spearman ρ: RC → SEEG Probability');
xlabel('Stimulation Site'); ylabel('\rho');
subplot(2,2,4)
bar(p_prob); ylim([0 max(p_prob)*1.1])
title('Row-wise Spearman p: RC → SEEG Probability');
xlabel('Stimulation Site'); ylabel('\p');
