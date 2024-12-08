clc;
clear all;

% Parameters
angles = [10, 25] / 180 * pi; % DOA (true angles in radians)
snapshot = 256; % Number of time samples
w = [pi/4, pi/3]'; % Frequencies of incident wavefronts
M = 5; % Number of array elements (sensors)
D = length(w); % Number of sources
lambda = 150; % Wavelength of signals
spacing = lambda / 2; % Spacing between elements in the ULA
snr = 0; % SNR in dB

% Generate the steering matrix
A = zeros(D, M);
for k = 1:D
    A(k, :) = exp(-1i * 2 * pi * spacing * sin(angles(k)) / lambda * (0:M-1));
end
A = A'; % Transpose to match dimensionality

% Generate the signals
F = 2 * exp(1j * (w * (1:snapshot))); % Incident signals
X = A * F; % Received signal at the array
X = X + awgn(X, snr); % Add Gaussian noise to simulate real-world scenario

% Estimate DOAs using MLE
grid_resolution = 0.5; % Grid resolution in degrees
est_DOAs = MLE(X, M, D, grid_resolution);

% Display the estimated DOAs
disp(['Estimated DOAs using MLE: ', num2str(est_DOAs)]);


