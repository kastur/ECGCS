% This script executes the function reconstructAmp M times to determine the
% sparsity-undersampling tradeoff as described in the paper "Message
% Passing Algorithms for Compressed Sensing" by Donoho et al.
% This code is open to anyone interested in using it.
% Ulugbek Kamilov, LTHC, EPFL, 2010.

% Prepare the workspace
clear; close all; home;

%--------------------------------------------------------------------------
% Set initial parameters
%--------------------------------------------------------------------------

% Dimension of the original signal
N = 1000;

% Undersampling grid
nDeltas = 50;
minDelta = 0.05;
maxDelta = 0.95;
ddelta = (maxDelta-minDelta)/(nDeltas-1);
deltaGrid = minDelta:ddelta:maxDelta;

% Sparsity grid
nRhos = 100;
minRho = 0.01;
maxRho = 0.99;
drho = (maxRho-minRho)/(nRhos-1);
rhoGrid = minRho:drho:maxRho;

% Number of iterations of the algorithm
T = 1000;

% Stopping criterion (tolerance for successfull decoding)
tol = 1e-4;

% Number of repetitions
M = 20;

% Array to store the phase transition values for each delta
rhos = zeros(nDeltas, 1);

%--------------------------------------------------------------------------
% Print to command line
%--------------------------------------------------------------------------
fprintf(1, '-----------------------------------\n');
fprintf(1, 'PhaseTransitionPlot: \n');
fprintf(1, '-----------------------------------\n');
fprintf(1, 'Parameters: \n');
fprintf(1, 'N = %d \n', N);
fprintf(1, 'tol = %f \n', tol);
fprintf(1, 'M = %d \n', M);
fprintf(1, 'nDeltas = %d \n', nDeltas);
fprintf(1, 'nRhos = %d \n', nRhos);

%--------------------------------------------------------------------------
% Perform simulations
%--------------------------------------------------------------------------
for i = 1:nDeltas

    % Set the delta
    d = deltaGrid(i);

    % Number of measurements corresponding to this delta
    n = ceil(d * N);

    % Print to command line
    fprintf(1, '-----------------------------------\n');
    fprintf(1, 'delta = %.4f\n', d);
    fprintf(1, '-----------------------------------\n');

    % Stores the rho giving the best tradeoff for this delta
    rhoopt = 0;

    for j = 1:nRhos

        % Set the rho
        r = rhoGrid(j);

        % Number of non-zero coefficients of the signal
        k = ceil(r * n);

        % Counter for the number successfull recoveries
        S = 0;

        % Counter for the number of failures
        F = 0;

        % Print to command line
        fprintf(1, 'rho = %.4f (maximum %d tries; o - pass, x - fail): ', r, M);

        for m = 1:M

            % Generate the measurement matrix
            A = (1/sqrt(n)) .* randn(n, N);

            % Sparse signal (with uniform distribution of non-zeros)
            x = [sign(rand(k,1) - 0.5); zeros(N-k,1)];
            x = x(randperm(N));

            % Generate Measurements
            y = A*x;

            % Recover x using AMP
            xhat = reconstructAmp(A, y, T, tol);

            % Verify if correct
            if(norm(y - A*xhat)/norm(y) < tol)
                S = S + 1;
                fprintf(1, 'o');
            else
                F = F + 1;
                fprintf(1, 'x');
            end

            % If too many errors or enough successes stop
            if(S/M >= 0.5 || F/M > 0.5)
                break;
            end
        end
        
        % Jump to the next line
        fprintf(1, '\n');

        % Save the phase transition value for this delta (or move on to the
        % next delta as we are in the failure zone)
        if(S/M >= 0.5 && r > rhoopt)
            rhoopt = r;
        else
            break;
        end
    end
    
    % Add optimal point for this delta to the curve
    rhos(i) = rhoopt;
end

%--------------------------------------------------------------------------
% Plot the curve
%--------------------------------------------------------------------------
figure;
plot(deltaGrid, rhos, 'LineWidth', 2);
title('Phase Transition curve');
xlabel('\delta');
ylabel('\rho');
xlim([0.05, 0.95]);
ylim([0, 1]);
grid on;