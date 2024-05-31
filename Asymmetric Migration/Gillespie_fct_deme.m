function rescue_status = Gillespie_fct_deme(Nit, fA, gA, NA10, NA20, fB, gB, NB10, NB20, K1, K2, mu, m12, m21, theta1, theta2)

  rescue_status = zeros(1, Nit);  % Vector to store rescue status for each replicate

  for i = 1:Nit
    % Initialization:
    NA1 = NA10;  % Initialization of the number of A individuals in deme 1
    NB1 = NB10;  % Initialization of the number of B individuals in deme 1
    NA2 = NA20;  % Initialization of the number of A individuals in deme 2
    NB2 = NB20;  % Initialization of the number of B individuals in deme 2
    t = 0;

    while (NA1 + NA2 ~= 0)
        % Compute the transition rates
        
         % Deme 1
         TrepA1 = (t < theta1) * fA * NA1 * (1 - mu); % as long deme 1 is intact we have reproduction of wt there
         TmutA1 = (t < theta1) * fA * NA1 * mu;       % as long deme 1 is intact we have mutations form wt to mutant
         TrepB1 = fB * NB1;                           % new number of mutant individuals in deme 1
         TdeathA1 = gA * NA1 * (NA1 + NB1) / K1;      % number of dead wt, according to total pop. size and carrying capacity of deme 1
         TdeathB1 = gB * NB1 * (NA1 + NB1) / K1;      % number of dead mutant, according to total pop. size and carrying capacity of deme 1
        
         % Deme 2
         TrepA2 = (t < theta2) * fA * NA2 * (1 - mu); % as long deme 2 is intact we have reproduction of wt there
         TmutA2 = (t < theta2) * fA * NA2 * mu;       % as long deme 2 is intact we have mutations form wt to mutant
         TrepB2 = fB * NB2;                           % new number of mutant individuals in deme 2
         TdeathA2 = gA * NA2 * (NA2 + NB2) / K2;      % number of dead wt, according to total pop. size and carrying capacity of deme 2
         TdeathB2 = gB * NB2 * (NA2 + NB2) / K2;      % number of dead mutant, according to total pop. size and carrying capacity of deme 2
        
         % Migration
         TmigA12 = NA1 * m12;            % number of A migrants from deme 1 to deme 2
         TmigB12 = NB1 * m12;            % number of B migrants from deme 1 to deme 2
         TmigA21 = NA2 * m21;            % number of A migrants from deme 2 to deme 1
         TmigB21 = NB2 * m21;            % number of B migrants from deme 2 to deme 1
       
         % Total Transitions, used later for probabilities for different actions
         T = TrepA1 + TmutA1 + TrepB1 + TdeathA1 + TdeathB1 + TrepA2 + TmutA2 + TrepB2 + TdeathA2 + TdeathB2 + TmigA12 + TmigB12 + TmigA21 + TmigB21;

      % Update time  
      r1 = rand(1);                        % generate one uniform random number; 0 < r1 < 1
      tau = (1 / T) * log(1 / r1);         % finding interval length (when will the next action happen) 
      t = t + tau;                         % Update timer: passed time plus new time interval

      % Build a sampling tower 
      ir2 = 1;			     % reset ir2 to the first position for each new screening of the sampling tower
      r2 = rand(1);                        % generate one uniform random number; 0 < r2 < 1
        cumul = zeros(1, 14);                % initialize cumul array
        cumul(1)  = TrepA1;
        cumul(2)  = TrepA1 + TmutA1;
        cumul(3)  = TrepA1 + TmutA1 + TrepB1;
        cumul(4)  = TrepA1 + TmutA1 + TrepB1 + TdeathA1;
        cumul(5)  = TrepA1 + TmutA1 + TrepB1 + TdeathA1 + TdeathB1;
        cumul(6)  = TrepA1 + TmutA1 + TrepB1 + TdeathA1 + TdeathB1 + TrepA2;
        cumul(7)  = TrepA1 + TmutA1 + TrepB1 + TdeathA1 + TdeathB1 + TrepA2 + TmutA2;
        cumul(8)  = TrepA1 + TmutA1 + TrepB1 + TdeathA1 + TdeathB1 + TrepA2 + TmutA2 + TrepB2;
        cumul(9)  = TrepA1 + TmutA1 + TrepB1 + TdeathA1 + TdeathB1 + TrepA2 + TmutA2 + TrepB2 + TdeathA2;
        cumul(10) = TrepA1 + TmutA1 + TrepB1 + TdeathA1 + TdeathB1 + TrepA2 + TmutA2 + TrepB2 + TdeathA2 + TdeathB2;
        cumul(11) = TrepA1 + TmutA1 + TrepB1 + TdeathA1 + TdeathB1 + TrepA2 + TmutA2 + TrepB2 + TdeathA2 + TdeathB2 + TmigA12;
        cumul(12) = TrepA1 + TmutA1 + TrepB1 + TdeathA1 + TdeathB1 + TrepA2 + TmutA2 + TrepB2 + TdeathA2 + TdeathB2 + TmigA12 + TmigB12;
        cumul(13) = TrepA1 + TmutA1 + TrepB1 + TdeathA1 + TdeathB1 + TrepA2 + TmutA2 + TrepB2 + TdeathA2 + TdeathB2 + TmigA12 + TmigB12 + TmigA21;
        cumul(14) = TrepA1 + TmutA1 + TrepB1 + TdeathA1 + TdeathB1 + TrepA2 + TmutA2 + TrepB2 + TdeathA2 + TdeathB2 + TmigA12 + TmigB12 + TmigA21 + TmigB21;  % corresponds to T

    % Determine which reaction occurs (recovery of the random number in the sampling tower vector) 
    while cumul(ir2) < r2 * T
        ir2 = ir2 + 1; 
    end
        
        % Update individuals according to the reaction 
        % --> screening the sampling tower
        if ir2 == 1
    		NA1 = NA1 + 1;         % Add one A individual to deme 1 (birth of NA1 individual)
	     elseif ir2 == 2
   		    NB1 = NB1 + 1;         % Add one B individual to deme 1 (mutation from A to B)
	     elseif ir2 == 3
    		NB1 = NB1 + 1;         % Add one A individual to deme 1 (birth of NB1 individual)
	     elseif ir2 == 4
    		NA1 = NA1 - 1;         % Subtract one A individual from deme 1 (death of NA1 individual)
	     elseif ir2 == 5
    		NB1 = NB1 - 1;         % Subtract one B individual from deme 1 (death of NB1 individual)
	     elseif ir2 == 6
   		    NA2 = NA2 + 1;         % Add one A individual to deme 2 (birth of NA2 individual)
	     elseif ir2 == 7
    		NB2 = NB2 + 1;         % Add one B individual to deme 2 (mutation from A to B)
	     elseif ir2 == 8
    		NB2 = NB2 + 1;         % Add one B individual to deme 2 (birth of NB1 individual)
	     elseif ir2 == 9
    		NA2 = NA2 - 1;         % Subtract one A individual from deme 2 (death of NA2 individual)
	     elseif ir2 == 10
    		NB2 = NB2 - 1;         % Subtract one A individual from deme 2 (death of NB2 individual)
	     elseif ir2 == 11
    		NA1 = NA1 - 1;         % Transfer one A individual from deme 1 to deme 2 (migration 1->2)
    		NA2 = NA2 + 1;
	     elseif ir2 == 12
    		NB1 = NB1 - 1;         % Transfer one B individual from deme 1 to deme 2 (migration 1->2)
    		NB2 = NB2 + 1;
	     elseif ir2 == 13
    		NA2 = NA2 - 1;         % Transfer one A individual from deme 2 to deme 1 (migration 2->1)
    		NA1 = NA1 + 1;
	     elseif ir2 == 14
    		NB2 = NB2 - 1;         % Transfer one B individual from deme 2 to deme 1 (migration 2->1)
    		NB1 = NB1 + 1;
        end
    end  
        
        final_population = [NA1, NB1, NA2, NB2];   % store the population size of the last interation
        
        % Apply rescue condition
        if final_population(2) + final_population(4) > 0  % if mutants are present at the last timepoint of the simulation 
            rescue_status(i) = 1;                           % we consider the population as rescued
        else                                              % if no mutants left
            rescue_status(i) = 0;                           % we consider the population as not rescued
        end
  end
    
  if Nit > 1  % this part is handy for the automated analysis. In any case we will have only number as an output
        rescue_status = mean(rescue_status);  % If we have more than 1 replicate we directly output the overall rescue probability
  end

end
