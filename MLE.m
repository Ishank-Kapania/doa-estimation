function est_DOAs = MLE(y, m, n, grid_resolution)
    % Implements the MLE algorithm using brute force grid search.
    %
    % Inputs:
    %   y               : Observed data matrix (m x N)
    %   m               : Number of sensors in the array
    %   n               : Number of sources
    %   grid_resolution : Grid resolution in degrees (e.g., 0.001)
    %
    % Outputs:
    %   est_DOAs        : Estimated Directions of Arrival (DOAs)

    % Generate finely spaced search grid
    search_grid = -90:grid_resolution:90;

    % Sample covariance matrix
    R = (y * y') / size(y, 2);

    % Initialize minimum likelihood function value
    min_FML = Inf;
    est_DOAs = zeros(1, n); % Placeholder for estimated DOAs

    % Iterate over all combinations of DOAs in the search grid
    combinations = nchoosek(search_grid, n); % All combinations of n DOAs
    num_combinations = size(combinations, 1);
    fprintf('Evaluating %d combinations...\n', num_combinations);

    for i = 1:num_combinations
        % Get current combination of DOAs
        angles = combinations(i, :);

        % Construct the steering matrix for the current angles
        A_test = exp(1j * (0:m-1)' * pi * sin(deg2rad(angles)));

        % Compute the likelihood function
        FML = trace(eye(m) - A_test * ((A_test' * A_test) \ A_test') * R);

        % Check if this is the minimum
        if FML < min_FML
            min_FML = FML;
            est_DOAs = angles; % Update estimated DOAs
        end
    end

    fprintf('MLE estimation completed.\n');
end
