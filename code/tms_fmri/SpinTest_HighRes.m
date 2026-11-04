permno = 1000;
wsname = 'rotationFS_test.mat';
load('vertex_indices_L.mat');
load('vertex_indices_R.mat');
vertex_indices_L = vertex_indices_L + 1;
vertex_indices_R = vertex_indices_R + 1;

SL_num = 29696;
SR_num = 29716;
down_SL_num = 3091;
down_SR_num = 3101;

% Initialize file names and labels
data_files = {'Imax.mat', 'Ipmax.mat', 'E1max.mat', 'Cpmax.mat'};
data_labels = {'Imax', 'Ipmax', 'E1max', 'Cpmax'};

% Loop through each data file and plot results
for data_idx = 1:length(data_files)
    % Load data
    data_file = ['./ODESolver/', data_files{data_idx}];
    data_struct = load(data_file);
    data_field = data_labels{data_idx};
    data_values = data_struct.(data_field);
    
    % Split into left and right hemisphere data
    datal = data_values(1:down_SL_num);
    datar = data_values(down_SL_num+1:down_SL_num+down_SR_num);

    load('left_sphere_32k.mat');
    load('right_sphere_32k.mat');
    load('group_indices.mat');
    verticesl = verticesl(vertex_indices_L, :);
    verticesr = verticesr(vertex_indices_R, :);

    % Downsample vertices
    num_groups = size(matching_groups_l, 2);
    verticesld = zeros(num_groups, 3);
    for i = 1:num_groups
        nodes_in_group_i = matching_groups_l{i};
        verticesld(i, :) = mean(verticesl(nodes_in_group_i, :), 1);
    end

    num_groups = size(matching_groups_r, 2);
    verticesrd = zeros(num_groups, 3);
    for i = 1:num_groups
        nodes_in_group_i = matching_groups_r{i} - SL_num;
        verticesrd(i, :) = mean(verticesr(nodes_in_group_i, :), 1);
    end

    load('nodes_fmri_values.mat');
    fMRI_data_full = fMRI_values;
    threshold = 59412;  % Set the threshold value
    filtered_groups = {};  % Initialize a cell array to store filtered groups
    filtered_centroids = [];
    fMRI_downsampled = [];
    %fMRI_downsampled = zeros(num_groups, 1);
    for group_id = 1:numel(group_indices)
        % Get the current group
        current_group = group_indices{group_id};

        % Check if all indices in the group satisfy the condition
        if all(current_group <= threshold)
            % If condition is satisfied, add the group to the filtered list
            filtered_groups{end+1} = current_group;
            filtered_centroids{end+1} = coords(group_id,:);
        end
    end
    num_groups = size(filtered_groups,2);
    for i = 1:num_groups
        nodes_in_group_i = filtered_groups{i};
        fMRI_downsampled(i) = mean(abs(fMRI_data_full(nodes_in_group_i))); % Use 'sum' if required
    end
    target_data = fMRI_downsampled;
    % Calculate correlations
    rng(0);
    bigrotl = [];
    bigrotr = [];
    I1 = eye(3, 3);
    I1(1, 1) = -1;
    bl = verticesld;
    br = verticesrd;

    for j = 1:permno
        A = normrnd(0, 1, 3, 3);
        [TL, temp] = qr(A);
        TL = TL * diag(sign(diag(temp)));
        if (det(TL) < 0)
            TL(:, 1) = -TL(:, 1);
        end
        TR = I1 * TL * I1;
        bl = bl * TL;
        br = br * TR;

        Il = nearestneighbour(verticesld', bl');
        Ir = nearestneighbour(verticesrd', br');

        bigrotl = [bigrotl; datal(Il)'];
        bigrotr = [bigrotr; datar(Ir)'];
    end

    save(wsname, 'bigrotl', 'bigrotr');
    [r(data_idx), p(data_idx)] = corr(cat(1, datal, datar), target_data', 'Type', 'Spearman');
    r_spin = zeros(permno, 4);
    for i = 1:permno
        [r_spin(i,data_idx), ~] = corr(cat(2, bigrotl(i, :), bigrotr(i, :))', target_data', 'Type', 'Spearman');
    end

    observed_corr(data_idx) = r(data_idx);
    null_dist = r_spin(:,data_idx);
    p_values_spin(data_idx) = mean((null_dist) >= (observed_corr(data_idx)));

    % Plot null distribution
    figure;
    histogram(null_dist, 30, 'Normalization', 'probability');
    hold on;

    % Plot observed correlation
    xline(observed_corr(data_idx), 'r', 'LineWidth', 2, 'DisplayName', 'Observed Correlation');
    title(sprintf('%s Null Distribution and Observed Correlation', data_field), 'FontSize', 20);
    xlabel('Correlation Coefficient', 'FontSize', 20);
    ylabel('Probability', 'FontSize', 20);
    subtitle(sprintf('r = %.3f, spin test p-value = %.3f', observed_corr(data_idx), p_values_spin(data_idx)),"FontSize",18);

    grid on;
    hold off;
end
