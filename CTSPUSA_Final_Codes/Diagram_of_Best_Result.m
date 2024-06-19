% Clear the environment
clear; clc; close all;

% Load the best results from the CSV file
best_results = readtable('Best_Distance_OUTPUT1.csv');

% Extract paths for both agents
agent1_path = str2num(best_results.Agent1_Path{1});
agent2_path = str2num(best_results.Agent2_Path{1});

% Load the distance matrix to get the coordinates of the cities
distance_matrix = csvread('sir_dataset.csv');
num_cities = size(distance_matrix, 1);

% Generate random coordinates for cities (assuming coordinates are not provided in the distance matrix)
coordinates = rand(num_cities, 2) * 100;

% Plot the paths
figure;
hold on;

% Plot the path for Agent 1
for i = 1:length(agent1_path)-1
    plot([coordinates(agent1_path(i), 1), coordinates(agent1_path(i+1), 1)], ...
         [coordinates(agent1_path(i), 2), coordinates(agent1_path(i+1), 2)], ...
         'r-', 'LineWidth', 2, 'Marker', '>');
end
plot([coordinates(agent1_path(end), 1), coordinates(agent1_path(1), 1)], ...
     [coordinates(agent1_path(end), 2), coordinates(agent1_path(1), 2)], ...
     'r-', 'LineWidth', 2, 'Marker', '>');

% Plot the path for Agent 2
for i = 1:length(agent2_path)-1
    plot([coordinates(agent2_path(i), 1), coordinates(agent2_path(i+1), 1)], ...
         [coordinates(agent2_path(i), 2), coordinates(agent2_path(i+1), 2)], ...
         'g-', 'LineWidth', 2, 'Marker', '>');
end
plot([coordinates(agent2_path(end), 1), coordinates(agent2_path(1), 1)], ...
     [coordinates(agent2_path(end), 2), coordinates(agent2_path(1), 2)], ...
     'g-', 'LineWidth', 2, 'Marker', '>');

% Plot the cities
scatter(coordinates(:, 1), coordinates(:, 2), 50, 'filled');
text(coordinates(:, 1), coordinates(:, 2), num2str((1:num_cities)'), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

% Highlight the starting cities
scatter(coordinates(agent1_path(1), 1), coordinates(agent1_path(1), 2), 100, 'k', 'filled');
scatter(coordinates(agent2_path(1), 1), coordinates(agent2_path(1), 2), 100, 'k', 'filled');

% Add legend
legend('Agent 1 Path', 'Agent 2 Path', 'Cities', 'Starting Cities');

% Add title and labels
title('Optimized Paths for Two Agents');
xlabel('X Coordinate');
ylabel('Y Coordinate');
hold off;
