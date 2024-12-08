function [estimated_DOAs, P_MODE1_dB] = mode1_estimation(X_noisy, M, K, d, lambda, theta_search, snr, snapshot)
    % Step 1: Sample Covariance Matrix
    R = (1 / snapshot) * (X_noisy * X_noisy');
    
    % Eigen-decomposition of R
    [E, D_eig] = eig(R);
    [D_eig_sorted, idx] = sort(diag(D_eig), 'descend');
    E = E(:, idx);

    % Signal and Noise Subspaces
    Es = E(:, 1:K); % Signal subspace
    En = E(:, K+1:end); % Noise subspace

    % Step 2: Obtain initial estimates using MUSIC
    P_MUSIC = zeros(size(theta_search));

    for i = 1:length(theta_search)
        theta = theta_search(i);
        a_theta = exp(-1i * 2 * pi * d * sin(theta * pi / 180) / lambda * (0:M-1)');
        P_MUSIC(i) = 1 / (a_theta' * (En * En') * a_theta);
    end

    % Find peaks in MUSIC spectrum
    [peaks, locs] = findpeaks(abs(P_MUSIC), 'NPeaks', K, 'SortStr', 'descend');
    initial_estimates = theta_search(locs);

    % Step 3: Compute the initial steering matrix using initial estimates
    A_init = zeros(M, K);
    for k = 1:K
        A_init(:, k) = exp(-1i * 2 * pi * d * sin(initial_estimates(k) * pi / 180) / lambda * (0:M-1)');
    end

    % Compute \hat{P} using the corrected computation
    P_hat_inv = inv(Es' * A_init);  % Corrected line

    % Step 4: Compute the MODE-1 Spectrum
    P_MODE1 = zeros(size(theta_search));

    for i = 1:length(theta_search)
        theta = theta_search(i);
        a_theta = exp(-1i * 2 * pi * d * sin(theta * pi / 180) / lambda * (0:M-1)');
        A_theta = a_theta;
        temp = (A_theta' * En) * (En' * A_theta) * P_hat_inv;
        F_MODE1 = trace(temp);
        P_MODE1(i) = abs(1 / F_MODE1); % Take reciprocal to get peaks
    end

    % Normalize and convert to dB scale
    P_MODE1_dB = 10 * log10(P_MODE1 / max(P_MODE1));

    % Find peaks in MODE-1 spectrum
    [peaks_MODE1, locs_MODE1] = findpeaks(P_MODE1_dB, theta_search, 'NPeaks', K, 'SortStr', 'descend');
    estimated_DOAs = locs_MODE1'; % Return the estimated DOAs
end
