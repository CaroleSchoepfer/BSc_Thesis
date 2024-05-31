tic % applied with "toc" at the very end to display how long the code took to run
% Here we perform an automated analysis of the Gillespie_fct_geno (function) with indirectly varying migration rates (T_mig and alpha change)

% set some parameters to fixed values
fA = 0.9;     % wild-type birth rate
gA = 1;       % wild-type death rate
NA10 = 90;    % initial number of wild-type individuals in deme 1
NA20 = 90;    % initial number of wild-type individuals in deme 2
fB = 0.9;     % mutant birth rate
gB = 1;       % mutant death rate
NB10 = 0;     % initial number of mutants in deme 1
NB20 = 0;     % initial number of mutants in deme 2
K1 = 100;     % carrying capacity of deme 1
K2 = 100;     % carrying capacity of deme 2
mu = 1/(50*(K1+K2));       % mutation probability upon reproduction
theta1 = 500;  % time at which deme 1 deteriorates
theta2 = 1000; % time at which deme 2 deteriorates


% set the desired number of replicates:
Nit = 3000;      % number of stochastic replicates (implemented within the function)

% set the parameters that should change value (this translates to varying migration rates in the function
alpha_values = [1/5, 1/2, 1, 2, 5];
Tmig_values = [1e-3 2e-3 5e-3 1e-2 2e-2 5e-2 1e-1 2e-1 5e-1 1e0 2e0 5e0 1e1 2e1 5e1 1e2 2e2 5e2 1e3 2e3 5e3 1e4 2e4];

% Initialize and define the data table - where and how to collect the results
num_values = numel(alpha_values) * numel(Tmig_values);
data_table = zeros(num_values, 6); % Preallocate data_table with zeros
row_counter = 0; % Initialize row counter

% Run the simulation across all chosen parameters
for alpha = alpha_values

    for Tmig = Tmig_values
        % Calculate migration rates based on current values of alpha and Tmig
        mWT = Tmig/(2*K1);   
        mM = alpha * mWT;          

        % Run the function defined above once (it will be automatically repeated "Nit"-times)
        analysis = Gillespie_fct_geno(Nit, fA, gA, NA10, NA20, fB, gB, NB10, NB20, K1, K2, mu, mWT, mM, theta1, theta2);

        % Increment row counter
        row_counter = row_counter + 1;
        
        % Enter the data into the table
        data_table(row_counter, :) = [mWT, mM, alpha, Tmig, mu, analysis];

        % display how many iterations of all planned iterations are done
        Say = sprintf('%d/%d', row_counter, length(Tmig_values)*length(alpha_values));
        disp(Say)

    end
end

% Convert data_table to a table
data_table = array2table(data_table, 'VariableNames', {'mWT', 'mM', 'alpha', 'Tmig', 'mu', 'rescue'});

% Display the last rows of the data table
disp(data_table(end-4:end, :));

% saving the data if desired
save('pResc_3000Rep_mu50_geno.mat', 'data_table');

toc
