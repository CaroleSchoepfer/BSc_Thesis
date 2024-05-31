% Load the data
load('pResc_3000Rep_mu50_geno.mat'); 

% Extract data
mWT = data_table.mWT;
mM = data_table.mM;
alpha = data_table.alpha;
Tmig = data_table.Tmig;
mu = data_table.mu;
rescue = data_table.rescue;

% Get unique alpha values
unique_alpha = unique(alpha);

% Plot the data
figure;
hold on;

for i = 1:length(unique_alpha)
    % Filter data for the current alpha value
    idx = alpha == unique_alpha(i);
    Tmig_alpha = Tmig(idx);
    rescue_alpha = rescue(idx);
    
    % Plot rescue values for the current alpha value with a unique color
    plot(Tmig_alpha, rescue_alpha, 'LineWidth', 0.7, 'DisplayName', ['\alpha = ' num2str(unique_alpha(i))]);
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
title('Genotype-dependent Migration');

% Create legend
legend('Location', 'Best');

% Hold off
hold off;
