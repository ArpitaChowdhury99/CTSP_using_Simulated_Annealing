% Clear the environment
clear; clc; close all;

% Parameters
benefit_per_city = 150; % Benefit for each city visited
cost_factor = 0.2; % Cost factor for the cost calculation

% Read the distance matrix from a CSV file
distance_matrix = csvread('sir_dataset.csv');
num_cities = size(distance_matrix, 1);
dataset_number = 1;

% Initialize benefits and cost matrix
benefits = repmat(benefit_per_city, num_cities, 1);
cost_matrix = distance_matrix * cost_factor;

% Initialize all pairs of starting cities
all_pairs = [];
for i = 1:num_cities
    for j = 1:num_cities
        if i ~= j
            all_pairs = [all_pairs; i, j];
        end
    end
end

% Initialize an empty table for results
results = [];

% Iterate through each starting pair
for pair_index = 1:size(all_pairs, 1)
    starting_cities = all_pairs(pair_index, :);
    
    % Optimize the division of this tour between the two agents using Simulated Annealing
    [agent1_path, agent2_path, agent1_benefit, agent2_benefit, agent1_cost, agent2_cost, agent1_timestamps, agent2_timestamps] = optimize_agents_paths_SA(benefits, cost_matrix, distance_matrix, starting_cities);

    % Calculate payoffs
    agent1_payoff = agent1_benefit - agent1_cost;
    agent2_payoff = agent2_benefit - agent2_cost;
    total_payoff = agent1_payoff + agent2_payoff;
    total_benefit = agent1_benefit + agent2_benefit;
    average_payoff = total_payoff / 2;

    % Count the number of cities visited by each agent
    agent1_city_count = numel(agent1_path);
    agent2_city_count = numel(agent2_path);

    % Convert agent paths to strings
    agent1_path_str = num2str(agent1_path, '%d,');
    agent1_path_str = agent1_path_str(1:end - 1); % Remove last comma
    agent2_path_str = num2str(agent2_path, '%d,');
    agent2_path_str = agent2_path_str(1:end - 1); % Remove last comma

    % Append results
    results = [results; {dataset_number, pair_index, starting_cities(1), starting_cities(2), agent1_path_str, agent2_path_str, agent1_benefit, agent2_benefit, agent1_cost, agent2_cost, agent1_payoff, agent2_payoff, total_benefit, total_payoff, average_payoff, agent1_city_count, agent2_city_count, join(string(agent1_timestamps), ','), join(string(agent2_timestamps), ',')}];
end

% Convert to a table for easy saving
resultsTable = cell2table(results, 'VariableNames', {'Dataset_No', 'Pair_Index', 'Starting_City_Agent1', 'Starting_City_Agent2', 'Agent1_Path', 'Agent2_Path', 'Agent1_Benefit', 'Agent2_Benefit', 'Agent1_Cost', 'Agent2_Cost', 'Agent1_Payoff', 'Agent2_Payoff', 'Total_Benefit', 'Total_Payoff', 'Average_Payoff', 'Agent1_City_Count', 'Agent2_City_Count', 'Agent1_Timestamps', 'Agent2_Timestamps'});

% Define the file path
filePath = 'Distance_OUTPUT_StartingPairs_Count_Graph.csv';

% Check if the file exists and append the results if it does
if exist(filePath, 'file')
    % File exists, read existing data
    existingData = readtable(filePath);
    % Append new results
    updatedResultsTable = [existingData; resultsTable];
else
    % File does not exist, use new results directly
    updatedResultsTable = resultsTable;
end

% Save the updated table to the CSV file
writetable(updatedResultsTable, filePath, 'WriteMode', 'overwrite');

fprintf('Results including city counts, payoffs, and timestamps have been updated in %s\n', filePath);

% Extract costs for each agent from the results
agent1_costs = cell2mat(results(:, 9));
agent2_costs = cell2mat(results(:, 10));

% Generate the indices for the x-axis
x_indices = 1:length(agent1_costs);

% Find the maximum cost and its index for each agent
[max_agent1_cost, max_agent1_index] = max(agent1_costs);
[max_agent2_cost, max_agent2_index] = max(agent2_costs);

% Plotting agent 1's costs
figure; % New figure for Agent 1
scatter(x_indices, agent1_costs, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'c', 'DisplayName', 'Agent 1 Cost');
hold on;
% Highlight the maximum cost
scatter(max_agent1_index, max_agent1_cost, 100, 'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'none', 'LineWidth', 2);
text(max_agent1_index, max_agent1_cost, sprintf('Max: Index %d', max_agent1_index), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
hold off;
title('Cost by Agent 1 for Each Pair Index');
xlabel('Pair Index');
ylabel('Cost');
legend show;

% Plotting agent 2's costs
figure; % New figure for Agent 2
scatter(x_indices, agent2_costs, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.4940 0.1840 0.5560], 'DisplayName', 'Agent 2 Cost');
hold on;
% Highlight the maximum cost
scatter(max_agent2_index, max_agent2_cost, 100, 'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'none', 'LineWidth', 2);
text(max_agent2_index, max_agent2_cost, sprintf('Max: Index %d', max_agent2_index), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
hold off;
title('Cost by Agent 2 for Each Pair Index');
xlabel('Pair Index');
ylabel('Cost');
legend show;

fprintf('Cost plots for each agent have been generated and displayed.\n');



% After all pair calculations are completed

% Extract agent payoffs from the results
agent1_payoffs = cell2mat(results(:, 11));
agent2_payoffs = cell2mat(results(:, 12));

% Find indices of maximum payoffs
[max_agent1_payoff, idx1] = max(agent1_payoffs);
[max_agent2_payoff, idx2] = max(agent2_payoffs);

% Plot for Agent 1
figure;
scatter(1:length(agent1_payoffs), agent1_payoffs, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1, 0.4, 0.6]); % Pink dots
hold on;
scatter(idx1, max_agent1_payoff, 100, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.5, 0, 0.5]); % Highlight maximum in purple
text(idx1, max_agent1_payoff, sprintf('Max: Pair %d', idx1), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
title('Agent 1 Payoffs');
xlabel('Pair Index');
ylabel('Payoff');
hold off;

% Plot for Agent 2
figure;
scatter(1:length(agent2_payoffs), agent2_payoffs, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.4, 1, 0.4]); % Green dots
hold on;
scatter(idx2, max_agent2_payoff, 100, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0, 0.5, 0]); % Highlight maximum in dark green
text(idx2, max_agent2_payoff, sprintf('Max: Pair %d', idx2), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
title('Agent 2 Payoffs');
xlabel('Pair Index');
ylabel('Payoff');
hold off;

% Save figures optionally
% saveas(gcf, 'agent2_payoffs.png');

fprintf('Payoff plots for each agent have been generated and displayed with maximum payoffs highlighted.\n');




% Simulated Annealing optimization function with fixed starting cities
function [agent1_path, agent2_path, agent1_benefit, agent2_benefit, agent1_cost, agent2_cost, agent1_timestamps, agent2_timestamps] = optimize_agents_paths_SA(benefits, cost_matrix, distance_matrix, starting_cities)
    num_cities = size(distance_matrix, 1);
    T_start = 1000; % Starting temperature
    T_end = 1; % Ending temperature
    alpha = 0.99; % Cooling rate
    max_iters = 10; % Maximum iterations at each temperature

    % Initialize the solution, reserving starting cities for agents
    current_solution = false(num_cities, 1);
    current_solution(starting_cities(1)) = true;
    current_solution(starting_cities(2)) = false;
    
    % Initialize a boolean array representing the assignment for remaining cities
    remaining_indices = setdiff(1:num_cities, starting_cities);
    remaining_solution = rand(length(remaining_indices), 1) > 0.5;
    current_solution(remaining_indices) = remaining_solution;

    current_fitness = calculate_agent_fitness(1:num_cities, current_solution, benefits, cost_matrix, distance_matrix, starting_cities);

    T = T_start;
    while T > T_end
        for iter = 1:max_iters
            % Generate a new solution by flipping one bit (within the remaining cities only)
            new_solution = current_solution;
            flip_index = remaining_indices(randi(length(remaining_indices)));
            new_solution(flip_index) = ~new_solution(flip_index);

            % Calculate the fitness of the new solution
            new_fitness = calculate_agent_fitness(1:num_cities, new_solution, benefits, cost_matrix, distance_matrix, starting_cities);

            % Acceptance probability
            if new_fitness > current_fitness
                current_solution = new_solution;
                current_fitness = new_fitness;
            else
                delta = current_fitness - new_fitness;
                if exp(-delta / T) > rand()
                    current_solution = new_solution;
                    current_fitness = new_fitness;
                end
            end
        end
        T = T * alpha; % Cool down
    end

    % Decode the final solution to paths, ensuring starting cities are the first in each path
    agent1_path = [starting_cities(1), remaining_indices(current_solution(remaining_indices))];
    agent2_path = [starting_cities(2), remaining_indices(~current_solution(remaining_indices))];

    % Calculate final benefits and costs using timestamps
    if ~isempty(agent1_path)
        [agent1_benefit, agent1_cost, agent1_timestamps] = calculate_benefit_and_cost_with_timestamps(agent1_path, agent2_path, benefits, cost_matrix, distance_matrix, starting_cities);
    else
        agent1_benefit = 0; agent1_cost = 0; agent1_timestamps = [];
    end
    if ~isempty(agent2_path)
        [agent2_benefit, agent2_cost, agent2_timestamps] = calculate_benefit_and_cost_with_timestamps(agent2_path, agent1_path, benefits, cost_matrix, distance_matrix, starting_cities);
    else
        agent2_benefit = 0; agent2_cost = 0; agent2_timestamps = [];
    end
end

% Calculate benefits, costs, and timestamps with full starting city benefits
function [total_benefit, total_cost, timestamps] = calculate_benefit_and_cost_with_timestamps(agent_path, other_agent_path, benefits, cost_matrix, distances, starting_cities)
    total_benefit = 0;
    total_cost = 0;
    timestamps = zeros(1, length(agent_path));
    if isempty(agent_path)
        return; % Skip processing for an empty path
    end
    timestamps(1) = 0; % Start from time 0
    other_timestamps = generate_timestamps(other_agent_path, distances);
    
    % Add full benefit for starting cities
    if ismember(agent_path(1), starting_cities)
        total_benefit = total_benefit + benefits(agent_path(1));
    end
    
    for i = 1:(length(agent_path) - 1)
        total_cost = total_cost + cost_matrix(agent_path(i), agent_path(i + 1));
        timestamps(i + 1) = timestamps(i) + distances(agent_path(i), agent_path(i + 1));
        
        % Determine benefit based on timestamps of the other agent
        idx_other = find(other_agent_path == agent_path(i + 1), 1);
        if ~isempty(idx_other) && other_timestamps(idx_other) == timestamps(i + 1)
            total_benefit = total_benefit + benefits(agent_path(i + 1)) / 2; % Share benefit
        elseif isempty(idx_other) || other_timestamps(idx_other) > timestamps(i + 1)
            total_benefit = total_benefit + benefits(agent_path(i + 1)); % Full benefit
        end
    end
    
    % Ensure agent_path has more than one city before calculating the return cost
    if length(agent_path) > 1
        total_cost = total_cost + cost_matrix(agent_path(end), agent_path(1));
        timestamps(end + 1) = timestamps(end) + distances(agent_path(end), agent_path(1));
    end
end

% Function to calculate the fitness of each agent using timestamps
function fitness = calculate_agent_fitness(complete_tour, division, benefits, cost_matrix, distances, starting_cities)
    % This function calculates the fitness of a given division of cities between two agents
    agent1_path = complete_tour(division == 1);
    agent2_path = complete_tour(division == 0);
    [agent1_benefit, agent1_cost] = calculate_benefit_and_cost_with_timestamps(agent1_path, agent2_path, benefits, cost_matrix, distances, starting_cities);
    [agent2_benefit, agent2_cost] = calculate_benefit_and_cost_with_timestamps(agent2_path, agent1_path, benefits, cost_matrix, distances, starting_cities);
    fitness = (agent1_benefit + agent2_benefit) - (agent1_cost + agent2_cost); % Objective function
end

% Function to generate timestamps for each agent's path
function timestamps = generate_timestamps(path, distances)
    timestamps = zeros(1, length(path));
    if isempty(path)
        return; % If the path is empty, return early
    end
    timestamps(1) = 0; % Start from time 0
    for i = 2:length(path)
        timestamps(i) = timestamps(i - 1) + distances(path(i - 1), path(i));
    end
end
