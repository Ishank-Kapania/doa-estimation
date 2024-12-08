clc
clear all

angles = [10, 25]/180*pi; % DOA
snapshot = 256; % Number of instants of time
w = [pi/4, pi/3]'; % Frequency of incident wavefronts
M = 5; % Number of array elements
D = length(w); % Number of incident wavefronts
lambda = 150; % Wavelength of incident wavefronts
spacing = lambda/2; % Spacing between array elements of ULA
snr = 0; % SNR
A = zeros(D,M);

for k=1:D
    A(k,:) = exp(-1i*2*pi*spacing*sin(angles(k))/lambda*(0:M-1)); % Steering/Mode vectors
end

A = A';
F = 2*exp(1j*(w*(1:snapshot))); % Incident signals
X = A*F;
X = X+awgn(X,snr); % Observed signal
S = X*X'; % Covariance
theta = -90:0.5:90; % Range of thetas to simulate for

Pmusic = MUSIC(S,M,D,lambda,spacing,theta); % The function P_{MU}(theta)

[Peaks,indices] = sort(Pmusic,'descend'); % Finding peaks
message = zeros(1,D);
for k = 1:D
   message(1,k) = theta(indices(k)); 
end
message = ['The Estimated angles are: ', int2str(message)]; % Displaying the obtained DOA
disp(message);