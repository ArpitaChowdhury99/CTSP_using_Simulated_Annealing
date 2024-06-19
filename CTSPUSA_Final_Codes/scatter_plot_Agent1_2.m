% Clear the environment
clear; clc; close all;

% Load the final results from the CSV file
final_results = readtable('Final_Results.csv');

% Extract relevant data for plotting
agent1_start_city = final_results.Agent1_Starting_City;
agent2_start_city = final_results.Agent2_Starting_City;

agent1_max_benefit = final_results.Agent1_Max_Benefit;
agent1_min_benefit = final_results.Agent1_Min_Benefit;
agent1_max_payoff = final_results.Agent1_Max_Payoff;
agent1_min_payoff = final_results.Agent1_Min_Payoff;

agent2_max_benefit = final_results.Agent2_Max_Benefit;
agent2_min_benefit = final_results.Agent2_Min_Benefit;
agent2_max_payoff = final_results.Agent2_Max_Payoff;
agent2_min_payoff = final_results.Agent2_Min_Payoff;

% Scatter plot for Agent1 benefits
figure;
scatter(agent1_start_city, agent1_max_benefit, 50, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'DisplayName', 'Agent1 Max Benefit');
hold on;
scatter(agent1_start_city, agent1_min_benefit, 50, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b', 'DisplayName', 'Agent1 Min Benefit');

% Mark the maximum and minimum benefit values
[~, idx_max_benefit] = max(agent1_max_benefit);
[~, idx_min_benefit] = min(agent1_min_benefit);
text(agent1_start_city(idx_max_benefit), agent1_max_benefit(idx_max_benefit), 'Max Benefit', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontWeight', 'bold', 'Color', 'r');
text(agent1_start_city(idx_min_benefit), agent1_min_benefit(idx_min_benefit), 'Min Benefit', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'FontWeight', 'bold', 'Color', 'b');

xlabel('Starting City');
ylabel('Benefit');
title('Agent1 Maximum and Minimum Benefits');
legend('show');
grid on;

% Scatter plot for Agent1 payoffs
figure;
scatter(agent1_start_city, agent1_max_payoff, 50, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'DisplayName', 'Agent1 Max Payoff');
hold on;
scatter(agent1_start_city, agent1_min_payoff, 50, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b', 'DisplayName', 'Agent1 Min Payoff');

% Mark the maximum and minimum payoff values
[~, idx_max_payoff] = max(agent1_max_payoff);
[~, idx_min_payoff] = min(agent1_min_payoff);
text(agent1_start_city(idx_max_payoff), agent1_max_payoff(idx_max_payoff), 'Max Payoff', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontWeight', 'bold', 'Color', 'r');
text(agent1_start_city(idx_min_payoff), agent1_min_payoff(idx_min_payoff), 'Min Payoff', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'FontWeight', 'bold', 'Color', 'b');

xlabel('Starting City');
ylabel('Payoff');
title('Agent1 Maximum and Minimum Payoffs');
legend('show');
grid on;

% Scatter plot for Agent2 benefits
figure;
scatter(agent2_start_city, agent2_max_benefit, 50, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'DisplayName', 'Agent2 Max Benefit');
hold on;
scatter(agent2_start_city, agent2_min_benefit, 50, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b', 'DisplayName', 'Agent2 Min Benefit');

% Mark the maximum and minimum benefit values
[~, idx_max_benefit] = max(agent2_max_benefit);
[~, idx_min_benefit] = min(agent2_min_benefit);
text(agent2_start_city(idx_max_benefit), agent2_max_benefit(idx_max_benefit), 'Max Benefit', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontWeight', 'bold', 'Color', 'r');
text(agent2_start_city(idx_min_benefit), agent2_min_benefit(idx_min_benefit), 'Min Benefit', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'FontWeight', 'bold', 'Color', 'b');

xlabel('Starting City');
ylabel('Benefit');
title('Agent2 Maximum and Minimum Benefits');
legend('show');
grid on;

% Scatter plot for Agent2 payoffs
figure;
scatter(agent2_start_city, agent2_max_payoff, 50, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'DisplayName', 'Agent2 Max Payoff');
hold on;
scatter(agent2_start_city, agent2_min_payoff, 50, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b', 'DisplayName', 'Agent2 Min Payoff');

% Mark the maximum and minimum payoff values
[~, idx_max_payoff] = max(agent2_max_payoff);
[~, idx_min_payoff] = min(agent2_min_payoff);
text(agent2_start_city(idx_max_payoff), agent2_max_payoff(idx_max_payoff), 'Max Payoff', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontWeight', 'bold', 'Color', 'r');
text(agent2_start_city(idx_min_payoff), agent2_min_payoff(idx_min_payoff), 'Min Payoff', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'FontWeight', 'bold', 'Color', 'b');

xlabel('Starting City');
ylabel('Payoff');
title('Agent2 Maximum and Minimum Payoffs');
legend('show');
grid on;
