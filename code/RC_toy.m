% Define the number of nodes in the RC network
n = 5;

% Initialize capacitance values for each node
C = ones(n, 1);  % In this case all nodes have capacitance = 1

%%%%%%%%%%%% Randomness C %%%%%%%%%%%%%%
% % Initialize random capacitance values for each node
% % Capacitance values range between 1 and 10 for variability
% C = 1 + (10 - 1) * rand(n, 1);  % Generate random capacitance for each node
%%%%%%%%%%%%

% Initialize resistance matrix with zeros
R = zeros(n, n);
% Define resistance values between connected nodes
% Connections between consecutive nodes with resistance = 1
R(1, 2) = 1;
R(2, 3) = 1;
R(3, 4) = 1;
R(4, 5) = 1;

% Set the diagonal (self-loops) to 0 to ensure no self-resistance at each node
for i = 1:n
    R(i, i) = 0;
end

% Make the resistance matrix symmetric (bidirectional connections)
for i = 1:n
    for j = i:n
        R(j, i) = R(i, j);
    end
end

%%%%%%%%%%%%%% Randomness R %%%%%%%%%%%%%%
% % Initialize resistance matrix with random values for each connection
% % Resistance values range between 0 and 10
% R = 0 + (10 - 0) * rand(n, n);

% % Ensure symmetry in the resistance matrix (for bidirectional connections)
% R = triu(R, 1) + triu(R, 1)';  % Keep only upper triangle and mirror it
% 
% % Set the diagonal (self-loops) to 0 to ensure no self-resistance at each node
% for i = 1:n
%     R(i, i) = 0;
% end

%%%%%%%%%%%%%%



% Set initial voltages for each node,in this case we choose the node 1 as
% input node with 1V source injected 
V0 = zeros(n, 1);
locationV = 1;       % Set initial voltage at Node 1 as TMS input
V0(locationV) = 1;   % Voltage 1V applied at Node 1



% Define simulation time range
T = [0 0.4];

% Solve the differential equations for the RC circuit with passed parameter
% R and C of this RC network
sol = ode45(@(t, v) RCtoy(t, v, R, C), T, V0);

% Define time points for analysis
t = linspace(0, 0.4, 100);
[v, vd] = deval(sol, t);

% Initialize current and positive current matrices
I = zeros(n, length(t));
Ip = zeros(n, length(t));

% Calculate current through each edge in the network over time
for i = 1:n
    for j = 1:n
        if i == j || R(j, i) == 0
            continue
        end
        I(i, :) = I(i, :) + (v(j, :) - v(i, :)) ./ R(j, i); % Total current through each node
        
        for k = 1:size(v, 2)
            if v(j, k) >= v(i, k)
                Ip(i, k) = Ip(i, k) + (v(j, k) - v(i, k)) / R(j, i); % Positive current only
            end
        end
    end
end

% Calculate power (P) and stored energy (Q) at each node over time
P = zeros(n, length(t));
Q = zeros(n, length(t));
Qd = zeros(n, length(t));

for i = 1:n
    P(i, :) = v(i, :) .* I(i, :);  % Power at each node
    Q(i, :) = 0.5 * C(i) .* (v(i, :)).^2;  % Stored energy
    Qd(i, :) = Q(i, :) - 0.5 * C(i) * (V0(i)^2); % Delta of stored energy from initial state
end

% Plot results: Voltage over time for each node
figure;
ledgName = {'Node2', 'Node3', 'Node4', 'Node5'};
for i = 2:n
    plot(t, v(i, :), 'LineWidth', 5);
    grid on;
    hold on;
end
legend(ledgName);
xlabel('Time (s)');
ylabel('Voltage (V)');
ax = gca;
ax.FontSize = 30;

% Plotting node 1(as input node) voltage separately in a small subplot
axes('Position', [.6 .25 .35 .35]);
box on;
plot(t, v(1, :), 'Color', "#77AC30", 'LineWidth', 5);
legend('Node1');
ax = gca;
ax.FontSize = 12;

% Plot current over time for each node
figure;
for i = 2:n
    plot(t, I(i, :), 'LineWidth', 5);
    grid on;
    hold on;
end
legend(ledgName);
xlabel('Time (s)');
ylabel('Current (A)');
ax = gca;
ax.FontSize = 30;

% Plot stored energy over time for each node
figure;
for i = 2:n
    plot(t, Q(i, :), 'LineWidth', 5);
    grid on;
    hold on;
end
legend(ledgName);
xlabel('Time (s)');
ylabel('Stored Energy (J)');
ax = gca;
ax.FontSize = 30;

% Plot node 1 stored energy separately in a small subplot
axes('Position', [.6 .25 .35 .35]);
box on;
plot(t, v(1, :), 'Color', "#77AC30", 'LineWidth', 5);
legend('Node1');
ax = gca;
ax.FontSize = 12;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function defining the differential equations for the RC circuit
function dvdt = RCtoy(t, v, R, C)
    n = length(v);  % Number of nodes
    dvdt = zeros(n, 1);  % Initialize the rate of voltage change for each node
    
    % Define the differential equations using Kirchhoff's Current Law
    for i = 1:n
        for j = 1:n
            if i ~= j && R(j, i) ~= 0
                dvdt(i) = dvdt(i) + (v(j) - v(i)) / (R(j, i) * C(i));
            end
        end
    end
end