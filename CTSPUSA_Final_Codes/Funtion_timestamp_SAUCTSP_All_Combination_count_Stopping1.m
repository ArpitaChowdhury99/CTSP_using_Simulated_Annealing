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

% Initialize variables to track max, min, and average benefits and payoffs
agent1_benefit_max = -Inf;
agent1_benefit_min = Inf;
agent1_benefit_sum = 0;
agent1_payoff_max = -Inf;
agent1_payoff_min = Inf;
agent1_payoff_sum = 0;
agent2_benefit_max = -Inf;
agent2_benefit_min = Inf;
agent2_benefit_sum = 0;
agent2_payoff_max = -Inf;
agent2_payoff_min = Inf;
agent2_payoff_sum = 0;

% Initialize variables for iteration and counter
max_iterations = 500;
counter = 0;
previous_best_cost = Inf;

% Iterate through each starting pair
for pair_index = 1:size(all_pairs, 1)
    starting_cities = all_pairs(pair_index, :);
    
    % Optimize the division of this tour between the two agents using Simulated Annealing
    [agent1_path, agent2_path, agent1_benefit, agent2_benefit, agent1_cost, agent2_cost, agent1_timestamps, agent2_timestamps, best_cost] = optimize_agents_paths_SA(benefits, cost_matrix, distance_matrix, starting_cities, max_iterations);

    % Calculate payoffs
    agent1_payoff = agent1_benefit - agent1_cost;
    agent2_payoff = agent2_benefit - agent2_cost;
    total_payoff = agent1_payoff + agent2_payoff;
    total_benefit = agent1_benefit + agent2_benefit;
    average_payoff = total_payoff / 2;

    % Update max, min, and sum for agent1
    if agent1_benefit > agent1_benefit_max
        agent1_benefit_max = agent1_benefit;
        agent1_benefit_max_city = starting_cities(1);
    end
    if agent1_benefit < agent1_benefit_min
        agent1_benefit_min = agent1_benefit;
        agent1_benefit_min_city = starting_cities(1);
    end
    agent1_benefit_sum = agent1_benefit_sum + agent1_benefit;
    if agent1_payoff > agent1_payoff_max
        agent1_payoff_max = agent1_payoff;
        agent1_payoff_max_city = starting_cities(1);
    end
    if agent1_payoff < agent1_payoff_min
        agent1_payoff_min = agent1_payoff;
        agent1_payoff_min_city = starting_cities(1);
    end
    agent1_payoff_sum = agent1_payoff_sum + agent1_payoff;

    % Update max, min, and sum for agent2
    if agent2_benefit > agent2_benefit_max
        agent2_benefit_max = agent2_benefit;
        agent2_benefit_max_city = starting_cities(2);
    end
    if agent2_benefit < agent2_benefit_min
        agent2_benefit_min = agent2_benefit;
        agent2_benefit_min_city = starting_cities(2);
    end
    agent2_benefit_sum = agent2_benefit_sum + agent2_benefit;
    if agent2_payoff > agent2_payoff_max
        agent2_payoff_max = agent2_payoff;
        agent2_payoff_max_city = starting_cities(2);
    end
    if agent2_payoff < agent2_payoff_min
        agent2_payoff_min = agent2_payoff;
        agent2_payoff_min_city = starting_cities(2);
    end
    agent2_payoff_sum = agent2_payoff_sum + agent2_payoff;

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

    % Check if the best cost improved
    if best_cost < previous_best_cost
        previous_best_cost = best_cost;
        counter = 0;
    else
        counter = counter + 1;
    end

    % Break the loop if no improvement for 50 iterations
    if counter > 50
        break;
    end
end

% Calculate averages
agent1_benefit_avg = agent1_benefit_sum / size(all_pairs, 1);
agent1_payoff_avg = agent1_payoff_sum / size(all_pairs, 1);
agent2_benefit_avg = agent2_benefit_sum / size(all_pairs, 1);
agent2_payoff_avg = agent2_payoff_sum / size(all_pairs, 1);

% Save the detailed results to the first CSV file
resultsTable = cell2table(results, 'VariableNames', {'Dataset_No', 'Pair_Index', 'Starting_City_Agent1', 'Starting_City_Agent2', 'Agent1_Path', 'Agent2_Path', 'Agent1_Benefit', 'Agent2_Benefit', 'Agent1_Cost', 'Agent2_Cost', 'Agent1_Payoff', 'Agent2_Payoff', 'Total_Benefit', 'Total_Payoff', 'Average_Payoff', 'Agent1_City_Count', 'Agent2_City_Count', 'Agent1_Timestamps', 'Agent2_Timestamps'});
writetable(resultsTable, 'Detailed_Results1.csv', 'WriteMode', 'overwrite');

% Save the summary results to the second CSV file
summaryResults = {
    'Agent1', agent1_benefit_max_city, agent1_benefit_max, agent1_benefit_min_city, agent1_benefit_min, agent1_benefit_avg, agent1_payoff_max_city, agent1_payoff_max, agent1_payoff_min_city, agent1_payoff_min, agent1_payoff_avg;
    'Agent2', agent2_benefit_max_city, agent2_benefit_max, agent2_benefit_min_city, agent2_benefit_min, agent2_benefit_avg, agent2_payoff_max_city, agent2_payoff_max, agent2_payoff_min_city, agent2_payoff_min, agent2_payoff_avg;
};
summaryTable = cell2table(summaryResults, 'VariableNames', {'Agent', 'Max_Benefit_City', 'Max_Benefit', 'Min_Benefit_City', 'Min_Benefit', 'Avg_Benefit', 'Max_Payoff_City', 'Max_Payoff', 'Min_Payoff_City', 'Min_Payoff', 'Avg_Payoff'});
writetable(summaryTable, 'Summary_Results1.csv', 'WriteMode', 'overwrite');

fprintf('Detailed results saved to Detailed_Results1.csv\n');
fprintf('Summary results saved to Summary_Results1.csv\n');

% Simulated Annealing optimization function with fixed starting cities
function [agent1_path, agent2_path, agent1_benefit, agent2_benefit, agent1_cost, agent2_cost, agent1_timestamps, agent2_timestamps, best_cost] = optimize_agents_paths_SA(benefits, cost_matrix, distance_matrix, starting_cities, max_iterations)
    num_cities = size(distance_matrix, 1);
    initial_temp = 100; % Initial temperature
    final_temp = 1; % Final temperature
    alpha = 0.95; % Cooling rate

    % Initialize solution
    current_solution = false(1, num_cities);
    remaining_indices = setdiff(1:num_cities, starting_cities);
    current_solution(remaining_indices) = rand(1, length(remaining_indices)) > 0.5;

    % Ensure starting cities are correctly assigned in the solution
    current_solution(starting_cities(1)) = true;
    current_solution(starting_cities(2)) = false;

    % Evaluate initial solution
    current_cost = calculate_agent_fitness(1:num_cities, current_solution, benefits, cost_matrix, distance_matrix, starting_cities);
    best_cost = current_cost;

    temp = initial_temp;

    % Simulated Annealing algorithm
    while temp > final_temp
        for iter = 1:max_iterations
            % Generate a new solution by flipping a bit
            new_solution = current_solution;
            idx = randi(num_cities);
            if ~ismember(idx, starting_cities)
                new_solution(idx) = ~new_solution(idx);
            end

            % Ensure starting cities are correctly assigned in the new solution
            new_solution(starting_cities(1)) = true;
            new_solution(starting_cities(2)) = false;

            % Evaluate new solution
            new_cost = calculate_agent_fitness(1:num_cities, new_solution, benefits, cost_matrix, distance_matrix, starting_cities);

            % Acceptance probability
            if new_cost > current_cost || exp((new_cost - current_cost) / temp) > rand()
                current_solution = new_solution;
                current_cost = new_cost;
                if new_cost > best_cost
                    best_cost = new_cost;
                end
            end
        end
        % Cool down the temperature
        temp = temp * alpha;
    end

    % Decode the final solution to paths, ensuring starting cities are the first in each path
    remaining_indices = setdiff(1:num_cities, starting_cities);
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
