clear all
%% Load coordinates and atlas labels
file_path = './mnipositionc.csv';
coords = readmatrix(file_path);  % Nx3 matrix of 3D node coordinates

label_data = readtable('Subcortex-Sina.csv', 'ReadVariableNames', false);
num_nodes = size(coords, 1);
labels = zeros(num_nodes, 1);  % Initialize labels array

% Assign labels to nodes based on atlas information
for i = 1:height(label_data)
    labels(label_data.Var2(i):label_data.Var4(i)) = i;
end

%% Adjustable parameters
k = 10;                 % Target group size
Dth = 5;               % Distance threshold (in mm)

%% Precompute k-d Tree for Fast Neighbor Search
kdtree = createns(coords, 'NSMethod', 'kdtree');

%% Initialize variables
visited = false(num_nodes, 1);  % Track visited nodes globally
group_indices = {};             % Store indices of nodes in each group
group_centroids = [];           % Store the centroid of each group

%% Cluster nodes within each region
unique_labels = unique(labels);
for region_idx = unique_labels'
    if region_idx == 0  % Skip unlabelled nodes
        continue;
    end

    % Find nodes in the current region
    region_nodes = find(labels == region_idx);
    
    region_coords = coords(region_nodes, :);  % Subset of coordinates for this region

    % Initialize visitation tracking for this region
    region_visited = false(length(region_nodes), 1);

    % While there are unvisited nodes in the region
    while any(~region_visited)
        %% Step 1: Start Node Selection
        if exist('next_start_idx', 'var') && ~isempty(next_start_idx)
            % Use the next nearest neighbor as the starting point
            start_idx = next_start_idx;
        else
            % Default to the first unvisited node
            start_idx = find(~region_visited, 1);
        end
        start_node_global_idx = region_nodes(start_idx);  % Global index of the node

        %% Step 2: Neighbor Search

        unvisited_indices = find(~region_visited);  % Indices of unvisited nodes within the region
        % Perform knnsearch only among unvisited nodes
        [local_neighbor_indices, distances] = knnsearch(region_coords(unvisited_indices, :), region_coords(start_idx, :), 'K', k);

        % Map local indices back to the global indices within the region
        neighbor_indices = unvisited_indices(local_neighbor_indices);

        % Filter neighbors by distance threshold
        valid_neighbors = neighbor_indices(distances <= Dth);
        group_node_indices = region_nodes(valid_neighbors);  % Global indices of valid neighbors

        %% Step 3: Centroid Computation
        group_coords = coords(group_node_indices, :);
        group_centroid = mean(group_coords, 1);

        %% Step 4: Merging Small Groups
        if numel(valid_neighbors) < k / 2 && ~isempty(group_centroids)
            % Check the distance to the last group's centroid
            last_centroid = group_centroids(end, :);
            dist_to_last_centroid = norm(group_centroid - last_centroid);

            % Merge if within distance threshold and size constraint
            if dist_to_last_centroid <= Dth && ...
               (numel(group_indices{end}) + numel(valid_neighbors)) <= round(3 * k / 2)
                % Merge current group into the last group
                group_indices{end} = union(group_indices{end}, group_node_indices);
                merged_coords = coords(group_indices{end}, :);
                group_centroids(end, :) = mean(merged_coords, 1);  % Update last group's centroid

                % Mark nodes as visited
                region_visited(valid_neighbors) = true;
                visited(group_node_indices) = true;

                % Skip creating a new group since it was merged
                % Set next_start_idx to the closest unvisited neighbor
                unvisited_indices = find(~region_visited);
                if ~isempty(unvisited_indices)
                    local_idx = knnsearch(region_coords(unvisited_indices, :), region_coords(start_idx, :), 'K', 1);
                    next_start_idx = unvisited_indices(local_idx);  % Map local index back to global index
                else
                    next_start_idx = [];
                end
                continue;
            end
        end

        %% Step 5: Update Status
        % Mark nodes in the current group as visited
        region_visited(valid_neighbors) = true;
        visited(group_node_indices) = true;

        % Store the group's indices and centroid
        group_indices{end + 1} = group_node_indices;  % Append group indices
        group_centroids = [group_centroids; group_centroid];  % Append group centroid

        %% Step 6: Update Next Start Node
        % Find the closest unvisited neighbor using k-d tree
        unvisited_indices = find(~region_visited);
        if ~isempty(unvisited_indices)
            % Perform knnsearch for the closest unvisited neighbor
            local_idx = knnsearch(region_coords(unvisited_indices, :), region_coords(start_idx, :), 'K', 1);
            next_start_idx = unvisited_indices(local_idx);  % Map local index back to global index
        else
            next_start_idx = [];  % No unvisited nodes remain
        end

    end
end

%% Save Results
% Save group_indices (list of all groups and their node indices)
save('group_indices.mat', 'group_indices');

% Save group centroids
save('group_centroids.mat', 'group_centroids');

% Optionally, save a node-to-group mapping
node_to_group = zeros(num_nodes, 1);
for group_id = 1:numel(group_indices)
    node_to_group(group_indices{group_id}) = group_id;
end
save('node_to_group_mapping.mat', 'node_to_group');

% Save group information in a human-readable format (CSV)
group_size = cellfun(@numel, group_indices);
group_table = table((1:numel(group_indices))', group_size', group_centroids(:, 1), group_centroids(:, 2), group_centroids(:, 3), ...
                    'VariableNames', {'GroupID', 'Size', 'CentroidX', 'CentroidY', 'CentroidZ'});

writetable(group_table, 'group_summary.csv');

disp('Downsampling completed and results saved.');
