% Set the number of cities
numCities = input('Enter the number of cities: ');
%numCities = 20;

% Set the size of the square area
areaSize = 200;

% Generate random coordinates for the cities within the square area
cityCoordinates = rand(numCities, 2) * areaSize;

% Calculate Euclidean distance between each pair of cities (vectorized)
xDiff = cityCoordinates(:, 1) - cityCoordinates(:, 1)';
yDiff = cityCoordinates(:, 2) - cityCoordinates(:, 2)';
distances = sqrt(xDiff.^2 + yDiff.^2);

% Round distances to the nearest integer
distances = round(distances);

% Save the distances matrix to a CSV file (using a relative path)
csvwrite('C:\Users\Pritoma Saha\Desktop\dataset\distance_matrix.csv', distances);

% Display the generated coordinates and distances
disp('City Coordinates:');
disp(cityCoordinates);
disp('Euclidean Distances (as integers):');
disp(distances);

% Plot the cities on a scatter plot
figure;
scatter(cityCoordinates(:, 1), cityCoordinates(:, 2), 'filled');
title('City Locations');
xlabel('X-coordinate');
ylabel('Y-coordinate');
axis([0 areaSize 0 areaSize]);
grid on;
