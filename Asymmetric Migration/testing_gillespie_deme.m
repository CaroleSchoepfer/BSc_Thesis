
% selecting parameters

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
m12 = 1e-3;    % migration rate from deme 1 to deme 2
m21 = 1e-3;    % migration rate from deme 2 to deme 1
theta1 = 500;  % time at which deme 1 deteriorates
theta2 = 1000; % time at which deme 2 deteriorates
Nit = 1;       % number of stochastic replicates

% test the function

result = Gillespie_fct_deme(Nit, fA, gA, NA10, NA20, fB, gB, NB10, NB20, K1, K2, mu, m12, m21, theta1, theta2);



