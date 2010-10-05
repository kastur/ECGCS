% Prepare the workspace
clear; close all; home;

% Dimensions of the problem
N = 2000;   % Signal length
n = 400;    % Number of measurements
k = 80;     % Number of non-zero elements

% Number of iterations
T = 1000;

% Tolerance
tol = 0.001;

% Generate the problem instance
A = (1/sqrt(n)).*randn(n, N);

% Sparse signal
x = sign(rand(k,1) - 0.5);
x = [x; zeros(N-k,1)];
x = x(randperm(N));

% Generate the measurement vector
y = A*x;

% Estimate using Iterative Soft Thresholding
%xist = RecommendedIST(A,y, T, tol);

% Estimate using Approximate Message Passing
xamp = reconstructAmp(A, y, T, tol, x, 0);

% Compute MSE
erramp = mean((xamp - x).*(xamp - x));

% Print the result
fprintf(1, 'Mean-Squared Error AMP: %.4f\n', erramp);