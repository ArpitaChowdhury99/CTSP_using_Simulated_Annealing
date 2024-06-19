% Clear the environment
clear; clc; close all;

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
final_results = cell(num_pairs, 18);

% Iterate over each pair and process the 20 CSV files
for pair_index = 1:num_pairs
    agent1_start_city = start_city_pairs(pair_index, 1);
    agent2_start_city = start_city_pairs(pair_index, 2);
    
    all_benefits_agent1 = [];
    all_benefits_agent2 = [];
    all_payoffs_agent1 = [];
    all_payoffs_agent2 = [];
    all_timestamps = [];

    for run = 1:20
        % Read the CSV file for the current run
        filePath = sprintf('All_Distance_OUTPUT_Run%d.csv', run);
        if exist(filePath, 'file')
            resultsTable = readtable(filePath);

            % Filter results for the current pair
            pair_results = resultsTable(resultsTable.Pair_Index == pair_index, :);
            
            % Collect benefits and payoffs
            if ~isempty(pair_results)
                all_benefits_agent1 = [all_benefits_agent1; pair_results.Agent1_Benefit];
                all_benefits_agent2 = [all_benefits_agent2; pair_results.Agent2_Benefit];
                all_payoffs_agent1 = [all_payoffs_agent1; pair_results.Agent1_Payoff];
                all_payoffs_agent2 = [all_payoffs_agent2; pair_results.Agent2_Payoff];

                % Collect timestamps
                if ~isempty(pair_results.Agent1_Timestamps{1})
                    agent1_timestamps = str2num(pair_results.Agent1_Timestamps{1})'; %#ok<ST2NM>
                else
                    agent1_timestamps = [];
                end
                if ~isempty(pair_results.Agent2_Timestamps{1})
                    agent2_timestamps = str2num(pair_results.Agent2_Timestamps{1})'; %#ok<ST2NM>
                else
                    agent2_timestamps = [];
                end
                % Ensure both are column vectors and handle empty cases
                agent1_timestamps = agent1_timestamps(:);
                agent2_timestamps = agent2_timestamps(:);
                combined_timestamps = [agent1_timestamps; agent2_timestamps];
                % Concatenate combined_timestamps to all_timestamps if all_timestamps is empty or has matching columns
                if isempty(all_timestamps)
                    all_timestamps = combined_timestamps';
                else
                    if size(all_timestamps, 2) == size(combined_timestamps, 1)
                        all_timestamps = [all_timestamps; combined_timestamps'];
                    end
                end
            end
        end
    end

    % Calculate statistics for agent1
    if ~isempty(all_benefits_agent1)
        max_benefit_agent1 = max(all_benefits_agent1);
        min_benefit_agent1 = min(all_benefits_agent1);
        avg_benefit_agent1 = mean(all_benefits_agent1);
        max_payoff_agent1 = max(all_payoffs_agent1);
        min_payoff_agent1 = min(all_payoffs_agent1);
        avg_payoff_agent1 = mean(all_payoffs_agent1);
        max_payoff_path_agent1_idx = find(all_payoffs_agent1 == max_payoff_agent1, 1);
        if ~isempty(max_payoff_path_agent1_idx) && max_payoff_path_agent1_idx <= height(pair_results)
            max_payoff_path_agent1 = pair_results.Agent1_Path{max_payoff_path_agent1_idx};
        else
            max_payoff_path_agent1 = '';
        end
    else
        max_benefit_agent1 = NaN;
        min_benefit_agent1 = NaN;
        avg_benefit_agent1 = NaN;
        max_payoff_agent1 = NaN;
        min_payoff_agent1 = NaN;
        avg_payoff_agent1 = NaN;
        max_payoff_path_agent1 = '';
    end
    
    % Calculate statistics for agent2
    if ~isempty(all_benefits_agent2)
        max_benefit_agent2 = max(all_benefits_agent2);
        min_benefit_agent2 = min(all_benefits_agent2);
        avg_benefit_agent2 = mean(all_benefits_agent2);
        max_payoff_agent2 = max(all_payoffs_agent2);
        min_payoff_agent2 = min(all_payoffs_agent2);
        avg_payoff_agent2 = mean(all_payoffs_agent2);
        max_payoff_path_agent2_idx = find(all_payoffs_agent2 == max_payoff_agent2, 1);
        if ~isempty(max_payoff_path_agent2_idx) && max_payoff_path_agent2_idx <= height(pair_results)
            max_payoff_path_agent2 = pair_results.Agent2_Path{max_payoff_path_agent2_idx};
        else
            max_payoff_path_agent2 = '';
        end
    else
        max_benefit_agent2 = NaN;
        min_benefit_agent2 = NaN;
        avg_benefit_agent2 = NaN;
        max_payoff_agent2 = NaN;
        min_payoff_agent2 = NaN;
        avg_payoff_agent2 = NaN;
        max_payoff_path_agent2 = '';
    end
    
    % Calculate timestamp statistics
    if ~isempty(all_timestamps)
        avg_timestamp = mean(all_timestamps, 'all');
        std_timestamp = std(all_timestamps, 0, 'all');
    else
        avg_timestamp = NaN;
        std_timestamp = NaN;
    end

    % Store results
    final_results(pair_index, :) = {
        agent1_start_city, max_benefit_agent1, min_benefit_agent1, avg_benefit_agent1, ...
        max_payoff_agent1, min_payoff_agent1, avg_payoff_agent1, max_payoff_path_agent1, ...
        agent2_start_city, max_benefit_agent2, min_benefit_agent2, avg_benefit_agent2, ...
        max_payoff_agent2, min_payoff_agent2, avg_payoff_agent2, max_payoff_path_agent2, ...
        avg_timestamp, std_timestamp
    };
end

% Convert to table and save to CSV
final_results_table = cell2table(final_results, 'VariableNames', {
    'Agent1_Starting_City', 'Agent1_Max_Benefit', 'Agent1_Min_Benefit', 'Agent1_Avg_Benefit', ...
    'Agent1_Max_Payoff', 'Agent1_Min_Payoff', 'Agent1_Avg_Payoff', 'Agent1_Max_Payoff_Path', ...
    'Agent2_Starting_City', 'Agent2_Max_Benefit', 'Agent2_Min_Benefit', 'Agent2_Avg_Benefit', ...
    'Agent2_Max_Payoff', 'Agent2_Min_Payoff', 'Agent2_Avg_Payoff', 'Agent2_Max_Payoff_Path', ...
    'Avg_Timestamp', 'Std_Timestamp'
});
writetable(final_results_table, 'Final_Results.csv', 'WriteMode', 'overwrite');
fprintf('Final results saved in Final_Results.csv\n');
