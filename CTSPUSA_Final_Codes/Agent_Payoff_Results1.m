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
final_results = cell(num_pairs, 4);

% Iterate over each pair and process the 20 CSV files
for pair_index = 1:num_pairs
    agent1_start_city = start_city_pairs(pair_index, 1);
    agent2_start_city = start_city_pairs(pair_index, 2);
    
    total_payoffs = [];
    
    for run = 1:20
        % Read the CSV file for the current run
        filePath = sprintf('All_Distance_OUTPUT_Run%d.csv', run);
        if exist(filePath, 'file')
            resultsTable = readtable(filePath);

            % Filter results for the current pair
            pair_results = resultsTable(resultsTable.Pair_Index == pair_index, :);
            
            % Collect total payoffs
            if ~isempty(pair_results)
                total_payoffs = [total_payoffs; pair_results.Total_Payoff];
            end
        end
    end

    % Calculate the maximum total payoff
    if ~isempty(total_payoffs)
        max_total_payoff = max(total_payoffs);
    else
        max_total_payoff = NaN;
    end

    % Store results
    final_results(pair_index, :) = {
        agent1_start_city, agent2_start_city, ...
        total_payoffs, max_total_payoff
    };
end

% Convert to table and save to CSV
final_results_table = cell2table(final_results, 'VariableNames', {
    'Agent1_Starting_City', 'Agent2_Starting_City', ...
    'Total_Payoffs', 'Max_Total_Payoff'
});
writetable(final_results_table, 'Agent_Payoff_Results1.csv', 'WriteMode', 'overwrite');
fprintf('Agent results saved in Agent_Payoff_Results1.csv\n');
