% Load the data
load('pResc_3000Rep_mu50_deme.mat'); 

% Extract data
m12 = data_table.m12;
m21 = data_table.m21;
beta = data_table.beta;
Tmig = data_table.Tmig;
mu = data_table.mu;
rescue = data_table.rescue;

% Get unique beta values
unique_beta = unique(beta);

% Plot the data
figure;
hold on;

for i = 1:length(unique_beta)
    % Filter data for the current beta value
    idx = beta == unique_beta(i);
    Tmig_beta = Tmig(idx);
    rescue_beta = rescue(idx);
    
    % Plot rescue values for the current beta value with a unique color
    plot(Tmig_beta, rescue_beta, 'LineWidth', 0.7, 'DisplayName', ['\beta = ' num2str(unique_beta(i))]);
end

% Set logarithmic scale for x-axis
set(gca, 'XScale', 'log');

% Set y-axis range from 0 to 1
ylim([0 1]);

% Add grid
grid on;

% Labels and title
xlabel('Tmig');
ylabel('Rescue');
title('Asymmetric Migration');

% Create legend
legend('Location', 'Best');

% Hold off
hold off;
