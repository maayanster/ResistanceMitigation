%% DESCRIPTION
% Comparison: runs both stochastic and deterministic models of
%                        insecticide resistance with differing rates of one 
%                        of the inputs (independent variable). Also
%                        produces a graph that compares the number of
%                        generations to resistance threshold and
%                        independent variable between the two models

%% INPUTS
Num_sim = 5;                                       % Number of simulations that comparison is run
q_freq_arr = [0.005, 0.005, 0.005, 0.005, 0.005];  % Initial frequency of resistant alleles
K_arr = [42000, 42000, 42000, 42000, 42000];       % Carrying capacity
Pref_arr = [0.3, 0.3, 0.3, 0.3, 0.3];              % Proportion of area that is refuge
WErr_ref_arr = [1, 0.5, 0.4, 0.3, 0.2];            % Fitness of RR in refuge with natural enemies
WErs_ref_arr = [1, 0.5, 0.5, 0.5, 0.5];            % Fitness of RS in refuge with natural enemies
WEss_ref_arr= [1, 0.5, 0.5, 0.5, 0.5];             % Fitness of SS in refuge with natural enemies
WErr_toxic_arr = [1, 1, 1, 1, 1];                  % Fitness of RR in toxic with natural enemies
WErs_toxic_arr = [1, 1, 1, 1, 1];                  % Fitness of RS in toxic with natural enemies
WEss_toxic_arr= [1, 1, 1, 1, 1];                   % Fitness of SS in toxic with natural enemies
gen_num_arr = [5000, 5000, 5000, 5000, 5000];      % Number of generations models are run];    

%% INITIALIZE
median_stochastic = zeros(1,length(q_freq_arr));
gen2thresh_mean_stochastic = zeros(1,length(q_freq_arr));
runs_below_thresh_array = zeros(1, length(q_freq_arr));

median_det50 = zeros(1,length(q_freq_arr));
gen2thresh_det_mean = zeros(1, length(q_freq_arr));

%% CALCULATIONS
% Call stochastic model
for mm = 1:Num_sim
    [gen2thresh_25percentile, gen2thresh_50percentile, ...
        gen2thresh_75percentile, gen2thresh_mean, runs_below_thresh]= ...
        Stochastic_Runs(q_freq_arr(mm), Pref_arr(mm), K_arr(mm), ...
        WErr_ref_arr(mm), WErs_ref_arr(mm), WEss_ref_arr(mm), ...
        WErr_toxic_arr(mm), WErs_toxic_arr(mm), WEss_toxic_arr(mm),...
        gen_num_arr(mm));
    median_stochastic(mm) = gen2thresh_50percentile;
    gen2thresh_mean_stochastic(mm) = gen2thresh_mean;
    runs_below_thresh_array(mm) = runs_below_thresh;
end

% Call deterministic model to match number of runs that reach threshold of
% resistance for the 50th percentile stochastic model
for mm = 1:length(median_stochastic)
    [q_freq, gen2thresh] = Deterministic_13March(q_freq_arr(mm),...
        Pref_arr(mm), WErr_ref_arr(mm), WErs_ref_arr(mm), ...
        WEss_ref_arr(mm), WErr_toxic_arr(mm), WErs_toxic_arr(mm),...
        WEss_toxic_arr(mm), gen_num_arr(mm));
    median_det50(mm) = gen2thresh;
end

% Call deterministic model to match number of runs that reach threshold of
% resistance for the mean stochastic model
for mm = 1:length(gen2thresh_mean_stochastic)
    [q_freq, gen2thresh] = Deterministic_13March(q_freq_arr(mm),...
        Pref_arr(mm), WErr_ref_arr(mm), WErs_ref_arr(mm), ...
        WEss_ref_arr(mm), WErr_toxic_arr(mm), WErs_toxic_arr(mm),...
        WEss_toxic_arr(mm), gen_num_arr(mm));
    gen2thresh_det_mean(mm) = gen2thresh;
end
%% COMPARE RESULTS
% Find ratio of stochastic over deterministic for the number of generations
% until resistant allele reaches threshold (frequency of 0.1)
Gen_thresh_50s_d = median_stochastic./median_det50;
Gen_thresh_means_d = gen2thresh_mean_stochastic./gen2thresh_det_mean;

% Assign parameter that is varied on x-axis
x = [0, 0.5, 0.6, 0.7, 0.8];

% Plot calculated stochastic over deterministic ratio of generations to 
% threshold (y_axis) and the parameter that is varied - x (x_axis)
figure
plot(x, Gen_thresh_50s_d)
xlabel('Mortality', 'FontSize', 12);
ylabel('Generations to threshold (stochastic over deterministic ratio)',...
    'FontSize', 12);
