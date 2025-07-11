MaxN = 29;  % Number of time points

% DATA
timeData = 1:MaxN;
larvaeData = [0,346.25,331,127.75,79.5,104.25,123.25,65.5,59.5,131.5,...
    105.75,57,61.25,204.25,155.5,72,75,228.5, 132.75, 96.5, 72.5, 358.75, 119, 37, 59.75, 255, 176, 76.50, 40];
pupaeData = [0,0,267.5,72.5,36.25,9.75,4.25,0,0,3,2,1.25,.75,9.25,4,...
    1,1.25,14.75, 6, 2.75, 0.25, 10.25, 21.75, 0.5, 0, 9.5, 20.25, 2.75, 3.75];
adultsData =[19.25,19.25,158,536.5,474.25,401.5,301,294.5,235.75,230.5,...
    209.5,208.5,150.5,145.25,112.5,109.25,98.5,90, 86, 74.25, 72.5, 62.75, 65.6, 56, 36.75, 44.5, 38, 32.75, 30.75];

% Combine data into a single vector for fitting
data = [larvaeData; pupaeData; adultsData];

% Reshape data to match model output format
data = reshape(data', [], 1);

% Fixed parameters
%mu_a = 0.006;
%mu_l = 0.0432;

% Function to fit
fitFunction = @(params, time) reshape(get_population_data(params, MaxN), [], 1);

% Initial guesses for the parameters
initialGuess = [0.005, 0.005, 0.005, 20, 0.005, 0.005];

% Bounds for the parameters
lb = [0, 0, 0, 11, 0, 0];
ub = [2, 2, 2, 100, 2, 2];

% Perform the curve fitting
paramsFit = lsqcurvefit(fitFunction, initialGuess, timeData, data, lb, ub);

% Extract the fitted parameters
c_ea = paramsFit(1);
c_el = paramsFit(2);
c_pa = paramsFit(3);
b = paramsFit(4);
mu_l = paramsFit(5);
mu_a = paramsFit(6);

% Display the fitted parameters
disp('Fitted parameters:');
disp(['c_ea = ', num2str(c_ea)]);
disp(['c_el = ', num2str(c_el)]);
disp(['c_pa = ', num2str(c_pa)]);
disp(['b = ', num2str(b)]);
disp(['mu_l = ', num2str(mu_l)]);
disp(['mu_a = ', num2str(mu_a)]);


% Generate the fitted model
LPAMatFit = simulate_population(paramsFit, MaxN);

% Plot the results
figure;
plot(LPAMatFit(:, 1), LPAMatFit(:, 2), 'LineWidth', 2, 'Color','b');
hold on;
plot(LPAMatFit(:, 1), LPAMatFit(:, 3), 'LineWidth', 2,'Color', 'r');
plot(LPAMatFit(:, 1), LPAMatFit(:, 4), 'LineWidth', 2,'Color', 'k');
plot(timeData, larvaeData,'o','MarkerFaceColor', 'b');
plot(timeData, pupaeData,'o','MarkerFaceColor', 'r');
plot(timeData, adultsData,'o','MarkerFaceColor', 'k');
title('Stage Structured Model Fit');
xlabel('Time (Weekly)');
ylabel('Population');
legend('Fitted Larvae', 'Fitted Pupae', 'Fitted Adults', 'Larvae Data', 'Pupae Data', 'Adults Data');
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

% Extract fitted data
fittedData = reshape(LPAMatFit(:, 2:4)', [], 1);

% Calculate residuals
residuals = data - fittedData;

% Calculate total sum of squares (SST)
SST = sum((data - mean(data)).^2);

% Calculate residual sum of squares (SSR)
SSR = sum(residuals.^2);

% Calculate R-squared
R_squared = 1 - (SSR / SST);

% Calculate RMSE
RMSE = sqrt(mean(residuals.^2));

% Display the goodness of fit metrics
disp(['R-squared: ', num2str(R_squared)]);
disp(['RMSE: ', num2str(RMSE)]);


% Define the model simulation function
function LPAMat = simulate_population(params, MaxN)
    c_ea = params(1);
    c_el = params(2);
    c_pa = params(3);
    b = params(4);
    mu_l = params(5);
    mu_a = params(6);
    
    MaxT = 29;
    ProjMatrix = zeros(29);
    ProjMatrix(1,29) = b;
    P = zeros(29, MaxN);
    P(29,1) = 50;
    
    for n = 1:MaxN
        ProjMatrix(2,1) = exp(-c_ea * P(29, n) - c_el * sum(P(4:11, n)));
        ProjMatrix(29,29) = (1 - mu_a);
        
        for i = 3:29
            if i < 5
                ProjMatrix(i, i-1) = exp(-c_ea * P(29, n) - c_el * sum(P(4:11, n)));
            elseif i < 13
                ProjMatrix(i, i-1) = (1 - mu_l);
            elseif i < 29
                ProjMatrix(i, i-1) = exp(-c_pa * P(29, n));
            elseif i == 29
                ProjMatrix(i, i-1) = exp(-c_pa * P(29, n));
            end
        end
        
        P(:, n+1) = ProjMatrix * P(:, n);
    end
    
    LPAMat = zeros(MaxN, 4);
    LPAMat(:, 1) = 1:MaxN;
    LPAMat(1, 4) = 50;
    
    for j = 2:MaxN
        LPAMat(j, 2) = sum(P(1:14, j));
        LPAMat(j, 3) = sum(P(15:28, j));
        LPAMat(j, 4) = sum(P(29, j));
    end
end
function populationData = get_population_data(params,  MaxN)
    LPAMat = simulate_population(params, MaxN);
    populationData = LPAMat(:, 2:4);  % Extract larvae, pupae, and adults columns
end