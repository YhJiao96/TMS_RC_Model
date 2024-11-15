% Load coordinates
file_path = './mnipositionc.csv';
your_coords = readmatrix(file_path);
coords = your_coords;   % Nx3 matrix of 3D node coordinates

% Adjustable parameters
k = 10;                 % Target cluster size
max_dist = 10;          % Maximum allowable distance within a cluster (in mm)

% Load atlas labels
label_data = readtable('Subcortex-Sina.csv', 'ReadVariableNames', false);
num_nodes = 91282;
labels = zeros(num_nodes, 1);  % Initialize labels array

% Assign labels to nodes based on start and end indices in the atlas file
for i = 1:height(label_data)
    start_idx = label_data.Var2(i);
    end_idx = label_data.Var4(i);
    labels(start_idx:end_idx) = i;  % Assign unique label for each region
end

% Get unique region labels
unique_labels = unique(labels);

% Initialize variables for tracking
visited = false(size(coords, 1), 1);  % Track visited nodes
downsampled_coords = [];  % Store downsampled coordinates here
group_indices = {};       % Cell array to store each downsampled group's node indices

% Loop over each region in the atlas
for region_idx = unique_labels'
    region_nodes = find(labels == region_idx);  % Nodes within this atlas region
    region_visited = visited(region_nodes);     % Track visited nodes within the region

    while any(~region_visited)
        % Find the next unvisited node within the region
        start_idx = region_nodes(find(~region_visited, 1));

        % Find k nearest neighbors
        [neighbors_idx, distances] = knnsearch(coords, coords(start_idx, :), 'K', k);

        % Filter neighbors by distance and within the same region
        valid_neighbors = neighbors_idx(distances <= max_dist & ismember(neighbors_idx, region_nodes));

        % Additional merging condition if the next group is smaller than k/2
        if numel(valid_neighbors) < k / 2
            % Look for additional nodes that are close enough to merge
            [extra_neighbors, extra_distances] = knnsearch(coords, coords(start_idx, :), 'K', k * 2);
            additional_nodes = extra_neighbors(extra_distances <= max_dist);
            if all(ismember(additional_nodes, region_nodes)) && numel(additional_nodes) > numel(valid_neighbors)
                valid_neighbors = [valid_neighbors; additional_nodes];
            end
        end

        % Mark these nodes as visited
        region_visited(valid_neighbors) = true;
        visited(valid_neighbors) = true;

        % Store indices of nodes in the current group
        group_indices{end + 1} = valid_neighbors;

        % Compute the centroid of this cluster
        cluster_coords = coords(valid_neighbors, :);
        centroid = mean(cluster_coords, 1);

        % Append centroid to downsampled coordinates
        downsampled_coords = [downsampled_coords; centroid];
    end
end


%%
% Downsample the SC matrix


load SCCount.mat; % load the SC matrix 
GroupSC_Count = spconvert(SCCount);
DiagGSC = diag(diag(GroupSC_Count));
SC = GroupSC_Count-DiagGSC;

num_groups = numel(group_indices);  % Number of downsampled groups
SC_downsampled = zeros(num_groups); % Downsampled SC matrix

% Loop over each pair of groups to compute the downsampled SC matrix
for i = 1:num_groups
    for j = i:num_groups
        % Get indices for nodes in each group
        nodes_in_group_i = group_indices{i};
        nodes_in_group_j = group_indices{j};
        
        % Extract submatrix for connections between group i and group j
        if i == j
            % Ignore intra-group connections by setting to zero
            SC_downsampled(i, j) = 0;
        else
            % Average SC values between nodes in group i and nodes in group j
            submatrix = SC(nodes_in_group_i, nodes_in_group_j);
            SC_downsampled(i, j) = mean(submatrix(:));
            SC_downsampled(j, i) = SC_downsampled(i, j); % Symmetric matrix
        end
    end
end


%%
% Downsample the Distance matrix

load SCLength.mat; % load the Distance matrix 
GroupSC_Length = spconvert(SCLength);
DiagGSL = diag(diag(GroupSC_Length));
D = GroupSC_Length - DiagGSL;

% Initialize the downsampled distance matrix
D_downsampled = zeros(num_groups);  % Downsampled Distance matrix

% Loop over each pair of groups to compute the downsampled distance matrix
for i = 1:num_groups
    for j = i:num_groups
        % Get indices for nodes in each group
        nodes_in_group_i = group_indices{i};
        nodes_in_group_j = group_indices{j};
        
        % Extract submatrix for distances between group i and group j
        if i == j
            % Ignore intra-group distances by setting to zero
            D_downsampled(i, j) = 0;
        else
            % Average distance values between nodes in group i and nodes in group j
            submatrix = D(nodes_in_group_i, nodes_in_group_j);
            D_downsampled(i, j) = mean(submatrix(:));
            D_downsampled(j, i) = D_downsampled(i, j); % Symmetric matrix
        end
    end
end


%%
% Downsample the Voltage in SimNIBS
% Load the original 91k voltage data
filename = 'subject_coord_91k_normE.csv';
normE = importdata(filename);  % 91,282-element array of voltages

% Initialize the downsampled voltage array
V_downsampled = zeros(num_groups, 1);  % Downsampled voltage array

% Loop through each group in group_indices to calculate the average voltage
for i = 1:num_groups
    % Get indices for nodes in the current group
    nodes_in_group = group_indices{i};
    
    % Compute the average voltage for this group
    V_downsampled(i) = mean(normE(nodes_in_group));
end


%%
% Save the variables into mat file
% Save Downsampled SC Matrix
save('SC_downsampled.mat', 'SC_downsampled');

% Save Downsampled Distance Matrix
save('D_downsampled.mat', 'D_downsampled');

% Save Downsampled Initial Voltage Vector
save('V_downsampled.mat', 'V_downsampled');