clc;
clear all;

% Define parameters
angles_true = [10, 25]; % True DOAs in degrees
K = length(angles_true); % Number of sources
angles_rad = angles_true * pi / 180; % Convert to radians

M = 5; % Number of array elements
d = 150 / 2; % Element spacing (half-wavelength)
lambda = 150; % Wavelength
snapshot = 256; % Number of snapshots
snr = 20; % Signal-to-Noise Ratio in dB

% Frequencies of the sources
w = [pi/4, pi/3]'; % Frequencies (arbitrary)

% Generate Steering Matrix for True DOAs
A = zeros(M, K);
for k = 1:K
    A(:, k) = exp(-1i * 2 * pi * d * sin(angles_rad(k)) / lambda * (0:M-1)');
end

% Generate signals
F = 2 * exp(1j * (w * (1:snapshot)));
X = A * F;

% Add noise
X_noisy = awgn(X, snr, 'measured');

% Define the search grid for MUSIC and MODE-1
theta_search = -90:0.1:90;

% Call the MODE-1 estimation function
[estimated_DOAs, P_MODE1_dB] = mode1_estimation(X_noisy, M, K, d, lambda, theta_search, snr, snapshot);

% Display the results
disp('Estimated DOAs from MODE-1:');
disp(estimated_DOAs);

