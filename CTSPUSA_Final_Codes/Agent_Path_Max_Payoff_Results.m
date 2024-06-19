% Define fixed starting city combinations for the agents
start_city_pairs = [
    1, 5;
    1, 10;
    1, 15;
    1, 20;
    5, 1;
    5, 10;
    5, 15;
    5, 20;
    10, 1;
    10, 5;
    10, 15;
    10, 20;
    15, 1;
    15, 5;
    15, 10;
    15, 20;
    20, 1;
    20, 5;
    20, 10;
    20, 15
];

num_pairs = size(start_city_pairs, 1);

% Initialize results table
path_results = cell(num_pairs, 3);

% Iterate over each pair and process the 20 CSV files
for pair_index = 1:num_pairs
    agent1_start_city = start_city_pairs(pair_index, 1);
    agent2_start_city = start_city_pairs(pair_index, 2);
    
    total_payoffs = [];
    agent1_paths = {};
    agent2_paths = {};
    
    for run = 1:20
        % Read the CSV file for the current run
        filePath = sprintf('All_Distance_OUTPUT_Run%d.csv', run);
        if exist(filePath, 'file')
            resultsTable = readtable(filePath);

            % Filter results for the current pair
            pair_results = resultsTable(resultsTable.Pair_Index == pair_index, :);
            
            % Collect total payoffs and paths
            if ~isempty(pair_results)
                total_payoffs = [total_payoffs; pair_results.Total_Payoff];
                agent1_paths = [agent1_paths; pair_results.Agent1_Path];
                agent2_paths = [agent2_paths; pair_results.Agent2_Path];
            end
        end
    end

    % Find the index of the maximum total payoff
    if ~isempty(total_payoffs)
        [max_total_payoff, max_index] = max(total_payoffs);
        max_agent1_path = agent1_paths{max_index};
        max_agent2_path = agent2_paths{max_index};
    else
        max_total_payoff = NaN;
        max_agent1_path = '';
        max_agent2_path = '';
    end

    % Store results
    path_results(pair_index, :) = {
        max_agent1_path, max_agent2_path, max_total_payoff
    };
end

% Convert to table and save to CSV
path_results_table = cell2table(path_results, 'VariableNames', {
    'Agent1_Path', 'Agent2_Path', 'Max_Total_Payoff'
});
writetable(path_results_table, 'Agent_Path_Max_Payoff_Results.csv', 'WriteMode', 'overwrite');
fprintf('Agent path and max payoff results saved in Agent_Path_Max_Payoff_Results.csv\n');
