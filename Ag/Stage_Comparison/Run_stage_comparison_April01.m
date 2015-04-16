%% Run_stage_comparison_April01
% Description:
%% INPUTS
Num_sim = 9;                                       % Number of simulations that comparison is run
q_freq_arr = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001];  % Initial frequency of resistant alleles
K_arr = [42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000];       % Carrying capacity
Pref_arr = [0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3];              % Proportion of area that is refuge
WErr_ref_arr = [1, 0.5, 0.4, 0.3, 0.2, 0.5, 0.4, 0.3, 0.2];            % Fitness of RR in refuge with natural enemies
WErs_ref_arr = [1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5];            % Fitness of RS in refuge with natural enemies
WEss_ref_arr= [1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.4, 0.3, 0.2];             % Fitness of SS in refuge with natural enemies
WErr_toxic_arr = [1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5];                  % Fitness of RR in toxic with natural enemies
WErs_toxic_arr = [1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5];                  % Fitness of RS in toxic with natural enemies
WEss_toxic_arr= [1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5];                   % Fitness of SS in toxic with natural enemies
gen_num_arr = [5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000];      % Number of generations models are run];    

%% INITIALIZE
median_generation2thresh_array = [];
runs_below_thresh_array = [];
median_times_crashed_array = [];
median_initial_end_array = [];
median_spread_end_array = [];
difference_gen2thresh_initial_spread = 0;
difference_array = [];

%% CALCULATIONS
% Call stochastic model
for mm = 1:Num_sim
    [median_generations2thresh, median_times_crashed, median_initial_end, runs_below_thresh, median_spread_end]= ...
        Run_stochastic_April01_stage(q_freq_arr(mm), Pref_arr(mm), K_arr(mm), ...
        WErr_ref_arr(mm), WErs_ref_arr(mm), WEss_ref_arr(mm), ...
        WErr_toxic_arr(mm), WErs_toxic_arr(mm), WEss_toxic_arr(mm),...
        gen_num_arr(mm));
    median_generation2thresh_array = [median_generation2thresh_array, median_generations2thresh];
    median_times_crashed_array = [median_times_crashed_array, median_times_crashed];
    median_initial_end_array = [median_initial_end_array, median_initial_end];
    median_spread_end_array = [median_spread_end_array, median_spread_end];
    runs_below_thresh_array = [runs_below_thresh_array, runs_below_thresh];
    difference_gen2thresh_initial_spread = [difference_gen2thresh_initial_spread, median_initial_end_array - median_spread_end];
end

