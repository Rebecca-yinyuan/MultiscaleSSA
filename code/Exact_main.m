function Exact_main(k)
% This function implements the novel exact algorithm based on
% the kth initial condition.

% To be submitted to the Spartan HPC in the University of Melbourne:
pc = parcluster('local');
pc.JobStorageLocation = getenv('SCRATCH');
pc.NumWorkers = 30; % 30 parallel workers


% BH IC: -----------------------------------------------
S0_values = [100, 200, 400, 800];
I0_values = [10, 20, 40, 80];
S0 = S0_values(k); I0 = I0_values(k);
SI0 = [S0; I0]; 

% BH parameters: -----------------------------------------------
BH_parms = struct;
BH_parms.mu = 10 ^ (-5);           % natural death rate
BH_parms.lambda = 10 ^ (-5);       % S introduction rate
% Coefficient linking the WH viral load to the individual disease
% transmission rate:
l_values = [10^(-9.8), 10^(-9.8), 10^(-9.8), 10^(-9.8)]; 
BH_parms.l = l_values(k);

% WH IC: -----------------------------------------------
T0_infect = 10 ^ (8);          % healthy target cells
Tstar0_infect = 10 ^ (0);      % infectious cells
V0_infect = 10 ^ (6);          % viral particles
TIV0_infect = [T0_infect; Tstar0_infect; V0_infect]; % this are initial
                               % conditions for the INFECTIOUS agents

% WH parameters: -----------------------------------------------
WH_parms = struct;
WH_parms.p = 10 ^ (1);           % viral production rate (p per day)
WH_parms.c = 10 ^ (0);           % viral clearance rate (c per day)
WH_parms.k = 10 ^ (-7);          % infection coefficient of target cells 
                                 % (k per virus per day)
WH_parms.mu_c = 10 ^ (-1);       % mortality rate of cells (mu_c per day)
WH_parms.delta_c = 10 ^ (1);     % extra mortality rate of infectious cells 
                                 % (delta_c per day)
WH_parms.R_0_WH = 1.1;           % within-host basic reproduction number
WH_parms.Lambda_c = WH_parms.R_0_WH * ...
    (WH_parms.mu_c * (WH_parms.mu_c + WH_parms.delta_c) * WH_parms.c)...
    / (WH_parms.k * WH_parms.p);

% Simulation setup: -----------------------------------------------
t0 = 0; 
t_end_values = [20, 20, 20, 20];
t_endBH = t_end_values(k);
t_endWH = t_endBH; 

% Resolution to extract the poplation-level SI dynamics over 800
% repetitions:
rep = 800; 
delta_t_uniform = 0.1;
t_uniform = (t0 : delta_t_uniform : t_endBH);

% The resolution to harvest the within-host information: ------------------
pt_N = 130;
delta_t_values = [0.0000937 * 1.084.^(0:pt_N);
    0.0000937 * 1.084.^(0:pt_N);
    0.0000937 * 1.084.^(0:pt_N);
    0.0000937 * 1.084.^(0:pt_N)];
% To construct a table to store 'Psi' info: -------------------------------
delta_t = delta_t_values(k, :);     


% Store the setup parameters and the IC: --------------------------------
setup_filename = ['exact_main_Setup_S', num2str(S0), 'I', ...
    num2str(I0), '_tend', num2str(t_endBH), '_rep', num2str(rep), '_l', ...
    num2str(BH_parms.l), '_k', num2str(k), '.mat'];
save(setup_filename, "BH_parms", "WH_parms", "SI0", "TIV0_infect", "t_endBH", ...
    "t0", "t_uniform", "t_endWH", "delta_t", "rep", 'k');

% RUN the novel exact algorithm: ------------------------------
exact_Run_deltat(BH_parms, SI0, t_endBH, t0, t_uniform, ...
    WH_parms, TIV0_infect, t_endWH, delta_t, rep, k)


end


function exact_Run_deltat(BH_parms, SI0, t_endBH, t0, t_uniform, ...
    WH_parms, TIV0_infect, t_endWH, delta_t, rep, k)
% This function runs the novel exact multi-scale algorithm for 'rep' times for 
% different delta t grid size (which is used to constuct the Psi table).

% Create an array to store the run times:
times_tabWH_exact = zeros(length(delta_t), 1); 
runs_elapsed_time_effiMulti = zeros(length(delta_t), 1); 

% To store the overall info at 't_uniform':
stocha_info_S_effiMulti_deltat = zeros(rep, length(t_uniform), length(delta_t)); 
stocha_info_I_effiMulti_deltat = zeros(rep, length(t_uniform), length(delta_t));

for i = 1 : length(delta_t)
    
    % Tabulate 'Psi-table':
    time_tabWH = tic; 
    [t_exact_WH, ~, Psi_table] = ...
        WH_exact(TIV0_infect, WH_parms, t_endWH, delta_t(i), BH_parms);
    times_tabWH_exact(i) = toc(time_tabWH);
    
    % Run rep times:
    runs_time_effiMulti = tic; 
    [stocha_info_S_exact, stocha_info_I_exact] = ...
        Exact_Run(BH_parms, SI0, t_endBH, t0, t_exact_WH, Psi_table, t_uniform, rep);
    runs_elapsed_time_effiMulti(i) = toc(runs_time_effiMulti);
    stocha_info_S_effiMulti_deltat(:, :, i) = stocha_info_S_exact;
    stocha_info_I_effiMulti_deltat(:, :, i) = stocha_info_I_exact;
end

% Store info into relevant .mat files:
EffiMultiscale_filename = ['ExactMulti_RunsDeltats_S', num2str(SI0(1)), 'I', ...
    num2str(SI0(2)), '_tendBH', num2str(t_endBH), '_rep', num2str(rep), '_l', ...
    num2str(BH_parms.l), '_k', num2str(k), '.mat'];
save(EffiMultiscale_filename, 't_uniform', 'stocha_info_S_effiMulti_deltat', ...
    'stocha_info_I_effiMulti_deltat', 'times_tabWH_exact', ...
    'runs_elapsed_time_effiMulti', 'k', '-v7.3');

end


function [stocha_info_S_exact, stocha_info_I_exact] = ...
    Exact_Run(BH_parms, SI0, t_endBH, t0, t_exact_WH, ...
    Psi_table, t_uniform, rep)
% This function runs the novel exact algorithm for 'rep' times.

stocha_info_S_exact = zeros(rep, length(t_uniform)); 
stocha_info_I_exact = zeros(rep, length(t_uniform));

parfor run = 1 : rep
    [S_uniform_exact, I_uniform_exact] = ...
    exact_oneRun(BH_parms, SI0, t_endBH, t0, t_exact_WH, Psi_table, t_uniform);
    % Store relevant info at 't_uniform':
    stocha_info_S_exact(run, :) = S_uniform_exact;
    stocha_info_I_exact(run, :) = I_uniform_exact;
end

end


function [S_uniform_exact, I_uniform_exact] = ...
    exact_oneRun(BH_parms, SI0, t_endBH, t0, t_exact_WH, ...
    Psi_table, t_uniform)
% This function runs the novel exact algorithm once.

[S_uniform_exact, I_uniform_exact] = ...
    Exact_Multiscale(BH_parms, SI0, t_endBH, t0, t_exact_WH, ...
    Psi_table, t_uniform);

end


function [BH_stocha_dynamics_t_array, BH_stocha_dynamics_t_index, ...
    BH_stocha_dynamics_S, BH_stocha_dynamics_I, BH_SItimelengths_Ilen, ...
    BH_SItimelengths_Slen, BH_event_infect, BH_event_intro_S, BH_event_Sdie, ...
    BH_event_Idie, BH_event_Q, BH_event_Psi] = ...
    init_multisystem(t0, SI0, BH_parms, Psi_table, t_WH)
% This function initialise the whole system (e.g. WH, BH systems, 
% inter-event time queues, and etc.).

% To store the BH t-S-I data: 
BH_stocha_dynamics_t_array = [];
BH_stocha_dynamics_t_index = 1; % the index of the event time array
BH_stocha_dynamics_t_array(BH_stocha_dynamics_t_index) = t0;
BH_stocha_dynamics_S = [];  BH_stocha_dynamics_I = []; 
BH_stocha_dynamics_S(BH_stocha_dynamics_t_index) = SI0(1); 
BH_stocha_dynamics_I(BH_stocha_dynamics_t_index) = SI0(2); 
% To store age of S or I:
BH_SItimelengths_Ilen = []; BH_SItimelengths_Slen = []; 
% To store the inter-event time queue corresponding to the following 4
% events in order:
% (1): a possible infection -- "infect"
% (2)-(3): one susceptible or infectious agent being removed from the 
% system, respectively -- "Sdie" or "Idie"
% (4): one susceptible individual entering the population -- "intro_S"
% We also use an array, "Q = [infect, Sdie, Idie, intro_S]", to store
% everything. 
BH_event_infect = []; BH_event_intro_S = [];
BH_event_Sdie = []; BH_event_Idie = []; 
BH_event_Q = [];
% To store the nu_value for each I:
BH_event_Psi = [];

% Initialise the inital event-waiting times:
for i = 1 : BH_stocha_dynamics_S(end)
    BH_SItimelengths_Slen(i) = BH_stocha_dynamics_t_array(end) - t0;
    BH_event_Sdie(i) = -log(rand) / BH_parms.mu;
end

for i = 1 : BH_stocha_dynamics_I(end)
    BH_SItimelengths_Ilen(i) = BH_stocha_dynamics_t_array(end) - t0;
    BH_event_Idie(i) = -log(rand) / BH_parms.mu;
    BH_event_Psi(i) = -log(rand) / BH_stocha_dynamics_S(end);
    % find the time length till the next infection based on 'Psi' value:
    [~, index] = min(abs(Psi_table - BH_event_Psi(i))); 
    BH_event_infect(i) = t_WH(index);
end
BH_event_intro_S(1) = -log(rand) / BH_parms.lambda;
BH_event_Q = [BH_event_infect, BH_event_Sdie, BH_event_Idie, BH_event_intro_S];

end


function [BH_stocha_dynamics_t_array, BH_stocha_dynamics_t_index, ...
    BH_stocha_dynamics_S, BH_stocha_dynamics_I, BH_SItimelengths_Ilen, ...
    BH_SItimelengths_Slen, BH_event_infect, BH_event_intro_S, BH_event_Sdie, ...
    BH_event_Idie, BH_event_Q, BH_event_Psi] = ...
    Event_Idie(DeltaT, removedI_index, Psi_table, t_WH, ...
    BH_stocha_dynamics_t_array, BH_stocha_dynamics_t_index, ...
    BH_stocha_dynamics_S, BH_stocha_dynamics_I, BH_SItimelengths_Ilen, ...
    BH_SItimelengths_Slen, BH_event_infect, BH_event_intro_S, BH_event_Sdie, ...
    BH_event_Idie, BH_event_Q, BH_event_Psi)
% This function lets the event of I removal happens and updates
% everything accordingly.

% Update BH population-level dynamics:
BH_stocha_dynamics_S(end + 1) = BH_stocha_dynamics_S(end); 
BH_stocha_dynamics_I(end + 1) = BH_stocha_dynamics_I(end) - 1;

% Update the time lengths of being S and I
BH_SItimelengths_Ilen = Delete_ele(removedI_index, BH_SItimelengths_Ilen')';
BH_SItimelengths_Ilen = BH_SItimelengths_Ilen + DeltaT;
BH_SItimelengths_Slen = BH_SItimelengths_Slen + DeltaT;

% Remove such I as it dies:
BH_event_infect = Delete_ele(removedI_index, BH_event_infect')';
BH_event_Psi = Delete_ele(removedI_index, BH_event_Psi')';
% Update the other infecious events' waiting times:
for i = 1 : BH_stocha_dynamics_I(end)
    % find nu_i_dash based on age of infection after this event of S intro:
    [~, index] = min(abs(t_WH - BH_SItimelengths_Ilen(i))); 
    nu_i_dash = Psi_table(index);
    BH_event_Psi(i) = BH_stocha_dynamics_S(end-1)/BH_stocha_dynamics_S(end)...
        *(BH_event_Psi(i)-nu_i_dash)+nu_i_dash;
    
    % update the time till next infection based on this new 'Psi' value:
    [~, index] = min(abs(Psi_table - BH_event_Psi(i))); 
    BH_event_infect(i) = t_WH(index) - BH_SItimelengths_Ilen(i);
end

% this event has already happened:
BH_event_Idie = Delete_ele(removedI_index, BH_event_Idie')';
% update event time queue:
BH_event_Q = [BH_event_infect, BH_event_Sdie, BH_event_Idie, BH_event_intro_S];

end


function [BH_stocha_dynamics_t_array, BH_stocha_dynamics_t_index, ...
    BH_stocha_dynamics_S, BH_stocha_dynamics_I, BH_SItimelengths_Ilen, ...
    BH_SItimelengths_Slen, BH_event_infect, BH_event_intro_S, BH_event_Sdie, ...
    BH_event_Idie, BH_event_Q, BH_event_Psi] = ...
    Event_Sdie(DeltaT, removedS_index, Psi_table, t_WH, ...
    BH_stocha_dynamics_t_array, BH_stocha_dynamics_t_index, ...
    BH_stocha_dynamics_S, BH_stocha_dynamics_I, BH_SItimelengths_Ilen, ...
    BH_SItimelengths_Slen, BH_event_infect, BH_event_intro_S, BH_event_Sdie, ...
    BH_event_Idie, BH_event_Q, BH_event_Psi)
% This function lets the event of S removal happens and updates
% everything accordingly.

% Update BH population-level dynamics:
BH_stocha_dynamics_S(end + 1) = BH_stocha_dynamics_S(end) - 1; 
BH_stocha_dynamics_I(end + 1) = BH_stocha_dynamics_I(end);

% Update the time lengths of being S and I
BH_SItimelengths_Slen = Delete_ele(removedS_index, BH_SItimelengths_Slen')'; 
BH_SItimelengths_Ilen = BH_SItimelengths_Ilen + DeltaT;
BH_SItimelengths_Slen = BH_SItimelengths_Slen + DeltaT;

% Update the other infecious events' waiting times:
for i = 1 : BH_stocha_dynamics_I(end)
    % Find nu_i_dash based on age of infection after this event of S intro:
    [~, index] = min(abs(t_WH - BH_SItimelengths_Ilen(i))); 
    nu_i_dash = Psi_table(index);
    BH_event_Psi(i) = BH_stocha_dynamics_S(end-1)/BH_stocha_dynamics_S(end)...
        *(BH_event_Psi(i)-nu_i_dash)+nu_i_dash;
    
    % Update the time till next infection based on this new 'Psi' value:
    [~, index] = min(abs(Psi_table - BH_event_Psi(i))); 
    BH_event_infect(i) = t_WH(index) - BH_SItimelengths_Ilen(i);
end

% This event has already happened:
BH_event_Sdie = Delete_ele(removedS_index, BH_event_Sdie')';
% Update event time queue:
BH_event_Q = [BH_event_infect, BH_event_Sdie, BH_event_Idie, BH_event_intro_S];

end


function [BH_stocha_dynamics_t_array, BH_stocha_dynamics_t_index, ...
    BH_stocha_dynamics_S, BH_stocha_dynamics_I, BH_SItimelengths_Ilen, ...
    BH_SItimelengths_Slen, BH_event_infect, BH_event_intro_S, BH_event_Sdie, ...
    BH_event_Idie, BH_event_Q, BH_event_Psi] = ...
    Event_Sintro(DeltaT, BH_parms, Psi_table, t_WH, ...
    BH_stocha_dynamics_t_array, BH_stocha_dynamics_t_index, ...
    BH_stocha_dynamics_S, BH_stocha_dynamics_I, BH_SItimelengths_Ilen, ...
    BH_SItimelengths_Slen, BH_event_infect, BH_event_intro_S, BH_event_Sdie, ...
    BH_event_Idie, BH_event_Q, BH_event_Psi)
% This function lets the event of S introduction happens and updates
% everything accordingly.

% Update BH population-level dynamics:
BH_stocha_dynamics_S(end + 1) = BH_stocha_dynamics_S(end) + 1; 
BH_stocha_dynamics_I(end + 1) = BH_stocha_dynamics_I(end);

% Update the time lengths of being S and I
BH_SItimelengths_Ilen = BH_SItimelengths_Ilen + DeltaT;
BH_SItimelengths_Slen = BH_SItimelengths_Slen + DeltaT;
BH_SItimelengths_Slen(end + 1) = 0; % this is for the newly introduced S

% Include this newly added S's estimated waiting time to be removed:
u = rand;
BH_event_Sdie(end + 1) = -log(u) / BH_parms.mu;

% Update the other infecious events' waiting times:
for i = 1 : BH_stocha_dynamics_I(end)
    % Find nu_i_dash based on age of infection after this event of S intro:
    [~, index] = min(abs(t_WH - BH_SItimelengths_Ilen(i))); 
    nu_i_dash = Psi_table(index);
    BH_event_Psi(i) = BH_stocha_dynamics_S(end-1)/BH_stocha_dynamics_S(end)...
        *(BH_event_Psi(i)-nu_i_dash)+nu_i_dash;
    
    % Update the time till next infection based on this new 'Psi' value:
    [~, index] = min(abs(Psi_table - BH_event_Psi(i))); 
    BH_event_infect(i) = t_WH(index) - BH_SItimelengths_Ilen(i);
end

% This event has already happened:
BH_event_intro_S = Delete_ele(1, BH_event_intro_S')';
% Update event time queue:
BH_event_Q = [BH_event_infect, BH_event_Sdie, BH_event_Idie, BH_event_intro_S];

end


function [BH_stocha_dynamics_t_array, BH_stocha_dynamics_t_index, ...
    BH_stocha_dynamics_S, BH_stocha_dynamics_I, BH_SItimelengths_Ilen, ...
    BH_SItimelengths_Slen, BH_event_infect, BH_event_intro_S, BH_event_Sdie, ...
    BH_event_Idie, BH_event_Q, BH_event_Psi] = ...
    Event_infection(DeltaT, activeI_index, BH_parms, t_WH, Psi_table, ...
    BH_stocha_dynamics_t_array, BH_stocha_dynamics_t_index, ...
    BH_stocha_dynamics_S, BH_stocha_dynamics_I, BH_SItimelengths_Ilen, ...
    BH_SItimelengths_Slen, BH_event_infect, BH_event_intro_S, BH_event_Sdie, ...
    BH_event_Idie, BH_event_Q, BH_event_Psi)
% This function lets the event of infection happens (i.e. infectious agent, 
% "activeI_index", randomly chooses a S to infect with some rejection 
% probability). It then updates everything accordingly.

% Update the time lengths of being S and I:
BH_SItimelengths_Ilen = BH_SItimelengths_Ilen + DeltaT;
BH_SItimelengths_Slen = BH_SItimelengths_Slen + DeltaT;

% The infected susceptible is randomly selected and deleted from the group
% of susceptible people: 
Sj = randi(BH_stocha_dynamics_S(end));
BH_event_Sdie = Delete_ele(Sj, BH_event_Sdie')';
BH_SItimelengths_Slen = Delete_ele(Sj, BH_SItimelengths_Slen')'; 

% Update BH population-level dynamics:
BH_stocha_dynamics_S(end + 1) = BH_stocha_dynamics_S(end) - 1; 
BH_stocha_dynamics_I(end + 1) = BH_stocha_dynamics_I(end) + 1;

% Add the new infectious individual to the system:
BH_event_Idie(end + 1) = -log(rand) / BH_parms.mu;
BH_SItimelengths_Ilen(end + 1) = 0; 
BH_event_Psi(end + 1) = -log(rand) / BH_stocha_dynamics_S(end);
% Find the time length till the next infection based on 'Psi' value:
[~, index] = min(abs(Psi_table - BH_event_Psi(end))); 
BH_event_infect(end + 1) = t_WH(index);

% Re-initialise this 'activeI_index' agent's inter-event time to next
% infection:
BH_event_Psi(activeI_index) = BH_event_Psi(activeI_index)-...
    log(rand)/BH_stocha_dynamics_S(end);
% Update the time till next infection based on this new 'Psi' value:
[~, index] = min(abs(Psi_table - BH_event_Psi(activeI_index))); 
BH_event_infect(activeI_index) = t_WH(index) - BH_SItimelengths_Ilen(activeI_index);

% Update the other infecious events' waiting times while excluding the 
% active infectious agent and the newly infectious person:
for i = 1 : BH_stocha_dynamics_I(end) - 1
    if i ~= activeI_index % exclude the person who just infects one S
        % Find Psi_i_dash based on age of infection after this event of S intro:
        [~, index] = min(abs(t_WH - BH_SItimelengths_Ilen(i))); 
        Psi_i_dash = Psi_table(index);
        BH_event_Psi(i) = BH_stocha_dynamics_S(end-1)/BH_stocha_dynamics_S(end)...
            *(BH_event_Psi(i)-Psi_i_dash)+Psi_i_dash;
        
        % Update the time till next infection based on this new 'Psi' value:
        [~, index] = min(abs(Psi_table - BH_event_Psi(i))); 
        BH_event_infect(i) = t_WH(index) - BH_SItimelengths_Ilen(i);
    end
end

% Update event time queue:
BH_event_Q = [BH_event_infect, BH_event_Sdie, BH_event_Idie, BH_event_intro_S];

end


function [BH_stocha_dynamics_t_array, BH_stocha_dynamics_t_index, ...
    BH_stocha_dynamics_S, BH_stocha_dynamics_I, BH_SItimelengths_Ilen, ...
    BH_SItimelengths_Slen, BH_event_infect, BH_event_intro_S, BH_event_Sdie, ...
    BH_event_Idie, BH_event_Q, BH_event_Psi] = ...
    choose_event(DeltaT, index_event, BH_parms, Psi_table, t_WH, ...
    BH_stocha_dynamics_t_array, BH_stocha_dynamics_t_index, ...
    BH_stocha_dynamics_S, BH_stocha_dynamics_I, BH_SItimelengths_Ilen, ...
    BH_SItimelengths_Slen, BH_event_infect, BH_event_intro_S, BH_event_Sdie, ...
    BH_event_Idie, BH_event_Q, BH_event_Psi)
% This function chooses event to happen depending on "DeltaT" and 
% "index_event"; lets event happen; and update all the information
% accordingly. 

% If there is no S left:
if BH_stocha_dynamics_S(end) <= 0 
    % S introduction happens:
    if index_event > BH_stocha_dynamics_I(end) 
        [BH_stocha_dynamics_t_array, BH_stocha_dynamics_t_index, ...
            BH_stocha_dynamics_S, BH_stocha_dynamics_I, BH_SItimelengths_Ilen, ...
            BH_SItimelengths_Slen, BH_event_infect, BH_event_intro_S, BH_event_Sdie, ...
            BH_event_Idie, BH_event_Q, BH_event_Psi] = ...
            Event_Sintro(DeltaT, BH_parms, Psi_table, t_WH, ...
            BH_stocha_dynamics_t_array, BH_stocha_dynamics_t_index, ...
            BH_stocha_dynamics_S, BH_stocha_dynamics_I, BH_SItimelengths_Ilen, ...
            BH_SItimelengths_Slen, BH_event_infect, BH_event_intro_S, BH_event_Sdie, ...
            BH_event_Idie, BH_event_Q, BH_event_Psi);
    % I removal happens:
    else 
        removedI_index = index_event;
        [BH_stocha_dynamics_t_array, BH_stocha_dynamics_t_index, ...
            BH_stocha_dynamics_S, BH_stocha_dynamics_I, BH_SItimelengths_Ilen, ...
            BH_SItimelengths_Slen, BH_event_infect, BH_event_intro_S, BH_event_Sdie, ...
            BH_event_Idie, BH_event_Q, BH_event_Psi] = ...
            Event_Idie(DeltaT, removedI_index, Psi_table, t_WH, ...
            BH_stocha_dynamics_t_array, BH_stocha_dynamics_t_index, ...
            BH_stocha_dynamics_S, BH_stocha_dynamics_I, BH_SItimelengths_Ilen, ...
            BH_SItimelengths_Slen, BH_event_infect, BH_event_intro_S, BH_event_Sdie, ...
            BH_event_Idie, BH_event_Q, BH_event_Psi);
    end
% There is S left:
else 
    % S introduction happens:
    if index_event > 2*BH_stocha_dynamics_I(end) + BH_stocha_dynamics_S(end) 
        [BH_stocha_dynamics_t_array, BH_stocha_dynamics_t_index, ...
            BH_stocha_dynamics_S, BH_stocha_dynamics_I, BH_SItimelengths_Ilen, ...
            BH_SItimelengths_Slen, BH_event_infect, BH_event_intro_S, BH_event_Sdie, ...
            BH_event_Idie, BH_event_Q, BH_event_Psi] = ...
            Event_Sintro(DeltaT, BH_parms, Psi_table, t_WH, ...
            BH_stocha_dynamics_t_array, BH_stocha_dynamics_t_index, ...
            BH_stocha_dynamics_S, BH_stocha_dynamics_I, BH_SItimelengths_Ilen, ...
            BH_SItimelengths_Slen, BH_event_infect, BH_event_intro_S, BH_event_Sdie, ...
            BH_event_Idie, BH_event_Q, BH_event_Psi);
    % S removal happens:
    elseif (BH_stocha_dynamics_I(end) + 1 <= index_event) && ...
            (index_event <= BH_stocha_dynamics_I(end) + BH_stocha_dynamics_S(end)) 
        removedS_index = index_event - BH_stocha_dynamics_I(end);
        [BH_stocha_dynamics_t_array, BH_stocha_dynamics_t_index, ...
            BH_stocha_dynamics_S, BH_stocha_dynamics_I, BH_SItimelengths_Ilen, ...
            BH_SItimelengths_Slen, BH_event_infect, BH_event_intro_S, BH_event_Sdie, ...
            BH_event_Idie, BH_event_Q, BH_event_Psi] = ...
            Event_Sdie(DeltaT, removedS_index, Psi_table, t_WH, ...
            BH_stocha_dynamics_t_array, BH_stocha_dynamics_t_index, ...
            BH_stocha_dynamics_S, BH_stocha_dynamics_I, BH_SItimelengths_Ilen, ...
            BH_SItimelengths_Slen, BH_event_infect, BH_event_intro_S, BH_event_Sdie, ...
            BH_event_Idie, BH_event_Q, BH_event_Psi);
     % I removal happens:
    elseif (BH_stocha_dynamics_I(end) + BH_stocha_dynamics_S(end) + 1 <= index_event) && ...
            (index_event <= 2*BH_stocha_dynamics_I(end) + BH_stocha_dynamics_S(end)) 
        removedI_index = index_event - BH_stocha_dynamics_I(end) - BH_stocha_dynamics_S(end);
        [BH_stocha_dynamics_t_array, BH_stocha_dynamics_t_index, ...
            BH_stocha_dynamics_S, BH_stocha_dynamics_I, BH_SItimelengths_Ilen, ...
            BH_SItimelengths_Slen, BH_event_infect, BH_event_intro_S, BH_event_Sdie, ...
            BH_event_Idie, BH_event_Q, BH_event_Psi] = ...
            Event_Idie(DeltaT, removedI_index, Psi_table, t_WH, ...
            BH_stocha_dynamics_t_array, BH_stocha_dynamics_t_index, ...
            BH_stocha_dynamics_S, BH_stocha_dynamics_I, BH_SItimelengths_Ilen, ...
            BH_SItimelengths_Slen, BH_event_infect, BH_event_intro_S, BH_event_Sdie, ...
            BH_event_Idie, BH_event_Q, BH_event_Psi);
    % Infection happens:
    else 
        activeI_index = index_event;
        [BH_stocha_dynamics_t_array, BH_stocha_dynamics_t_index, ...
            BH_stocha_dynamics_S, BH_stocha_dynamics_I, BH_SItimelengths_Ilen, ...
            BH_SItimelengths_Slen, BH_event_infect, BH_event_intro_S, BH_event_Sdie, ...
            BH_event_Idie, BH_event_Q, BH_event_Psi] = ...
            Event_infection(DeltaT, activeI_index, BH_parms, t_WH, Psi_table, ...
            BH_stocha_dynamics_t_array, BH_stocha_dynamics_t_index, ...
            BH_stocha_dynamics_S, BH_stocha_dynamics_I, BH_SItimelengths_Ilen, ...
            BH_SItimelengths_Slen, BH_event_infect, BH_event_intro_S, BH_event_Sdie, ...
            BH_event_Idie, BH_event_Q, BH_event_Psi);
    end
end

% Althouth all information has been updated in those "Event_" subfunctions,
% we still perform one overall update depending on whether we've run out S.
% Update the waiting time Q:
if BH_stocha_dynamics_S(end) == 0 % no S left
    BH_event_infect = [];
    BH_event_Sdie = [];
    % note that by setting the event time array "BH_event.infect = []", we
    % no longer to update "BH_event.ProbInfect" properly as the infection
    % is not allowed to happen in the first place
end
BH_event_Q = [BH_event_infect, BH_event_Sdie, BH_event_Idie, BH_event_intro_S];

end


function [S_uniform, I_uniform, SI_lengths_uniform] = ...
    Exact_Multiscale(BH_parms, SI0, t_endBH, t0, t_WH, ...
    Psi_table, t_uniform)
% This function runs the novel exact algorithm once.

% Initilise the system:
[BH_stocha_dynamics_t_array, BH_stocha_dynamics_t_index, ...
    BH_stocha_dynamics_S, BH_stocha_dynamics_I, BH_SItimelengths_Ilen, ...
    BH_SItimelengths_Slen, BH_event_infect, BH_event_intro_S, BH_event_Sdie, ...
    BH_event_Idie, BH_event_Q, BH_event_Psi] = init_multisystem(t0, ...
    SI0, BH_parms, Psi_table, t_WH); % checked

while (BH_stocha_dynamics_t_array(end) < t_endBH)
    % We still have events to happen:
    if ~ isempty(BH_event_Q) 
        % Find the min value in the waiting time queue Q:
        [DeltaT, index_event] = min(BH_event_Q);
        %if BH_stocha_dynamics.t_array(end) + DeltaT > t_endBH
        %    break;
        %else
            % Let the corresponding event happen and update the event queues:
            [BH_stocha_dynamics_t_array, BH_stocha_dynamics_t_index, ...
                BH_stocha_dynamics_S, BH_stocha_dynamics_I, BH_SItimelengths_Ilen, ...
                BH_SItimelengths_Slen, BH_event_infect, BH_event_intro_S, BH_event_Sdie, ...
                BH_event_Idie, BH_event_Q, BH_event_Psi] = ...
                choose_event(DeltaT, index_event, BH_parms, Psi_table, t_WH, ...
                BH_stocha_dynamics_t_array, BH_stocha_dynamics_t_index, ...
                BH_stocha_dynamics_S, BH_stocha_dynamics_I, BH_SItimelengths_Ilen, ...
                BH_SItimelengths_Slen, BH_event_infect, BH_event_intro_S, BH_event_Sdie, ...
                BH_event_Idie, BH_event_Q, BH_event_Psi);

            % Step through time:
            BH_stocha_dynamics_t_index = BH_stocha_dynamics_t_index + 1;
            BH_stocha_dynamics_t_array(BH_stocha_dynamics_t_index) = ...
                BH_stocha_dynamics_t_array(BH_stocha_dynamics_t_index - 1) + DeltaT;

        %end
    end
end

% Post-simulation SI infomation harvesting at 't_uniform':
S_uniform = zeros(size(t_uniform))';
I_uniform = zeros(size(t_uniform))';
SI_lengths_uniform = cell(size(t_uniform))';
for m = 1 : length(t_uniform)
    t_m = t_uniform(m);
    t_stocha_index = max(find(BH_stocha_dynamics_t_array <= t_m));
    % The SI dyanmics at t = t_i are the same as the ones at
    % t_stocha_index:
    S_uniform(m) = BH_stocha_dynamics_S(t_stocha_index);
    I_uniform(m) = BH_stocha_dynamics_I(t_stocha_index);
end

end


function matrix = Delete_ele(index, matrix)
% Remove the 'index_th' row in the matrix (note that if "matrix" is an array, 
% we need have a column vector!!)
if size(matrix, 1) == 1
    matrix = [];
else
    if index > 1 && index < size(matrix, 1)
        matrix = [matrix(1:index-1, :); matrix(index+1:size(matrix, 1), :)];
    else
        if index == 1
            matrix = matrix(index+1:size(matrix, 1), :);
        elseif index == size(matrix, 1)
            matrix = matrix(1:size(matrix, 1)-1, :);
        end
    end
end

end


function dydt = WH_infect_odesystem_Log(t, y, WH_parms)
% The ODE equations describing the dynamics of WH subsystem among I populations. 
% Note that:
% y(1) := T
% y(2) := Tstar
% y(3) := V
dydt = [
    WH_parms.Lambda_c*exp(-y(1)) - WH_parms.k*exp(y(3)) - WH_parms.mu_c;
    WH_parms.k*exp(y(1)+y(3)-y(2)) - WH_parms.mu_c - WH_parms.delta_c;
    -WH_parms.c + WH_parms.p*exp(y(2)-y(3))
];
end


function [t_WH, V_infect, Psi_table] = ...
    WH_exact(TIV0_infect, WH_parms, t_endWH, delta_t, BH_parms)
% This function solves the TIV model for the infectious population and
% outputs V_infect depending on the age of infection from 0 to 't_endWH'.

Initial_Conds_infect = log(TIV0_infect); % log trandsform the ICs to solve ODE
t_WH = (0 : delta_t: t_endWH);
opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
[~, y1] = ode15s(@(t,y) WH_infect_odesystem_Log(t, y, WH_parms), t_WH, ...
    Initial_Conds_infect, opts);
V_infect_Log = y1(: , 3); 
V_infect = exp(V_infect_Log);

% Record the 'Psi' table:
Psi_table = cumtrapz(t_WH, V_infect) * BH_parms.l;

end

