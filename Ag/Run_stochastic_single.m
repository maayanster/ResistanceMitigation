%% DESCRIPTION of Run_stochastic_single 
% Runs a stochastic model of insecticide resistance with a single set of 
% parameters multiple times. The script also plots a histogram of the
% number of generations to resistance for however many times the stochastic
% model is run
%% Structure of stochastic function:
%[generations2thresh, q_array] = Stochastic(q_freq, Pref, K, ...
% WErr_ref, WErs_ref, WEss_ref, WErr_toxic, WErs_toxic, WEss_toxic,...
% gen_num)
% Where
% q_freq -----> Frequency of resistant alleles
% Pref -------> Proportion of area that is refuge (no Bt)
% K ----------> Carrying capacity
% WErr_ref ---> Fitness of RR in refuge with natural enemies
% WErs_ref ---> Fitness of RS in refuge with natural enemies
% WEss_ref ---> Fitness of SS in refuge with natural enemies
% WErr_toxic -> Fitness of RR in field due to natural enemies
% WErs_toxic -> Fitness of RS in field due to natural enemies
% WEss_toxic -> Fitness of SS in field due to natural enemies
% gen_num ----> Number of generations 

%% INPUTS
num_runs = 5000;

%% INITIALIZE
q_freq_arr = zeros(1,num_runs);
gen2thresh_arr = zeros(1,num_runs);

%% Call stochastic function multiple times
for nn=1:num_runs                                        % number of trials
    [generation2thresh, q_freq] = Stochastic(0.005,0.3,42000,...
        0.5, 0.5, 0.5, 0.5, ...
        0.5, 0.5, 2500);
    gen2thresh_arr(nn) = generation2thresh;
end

%% Generate histogram of generations to threshold of multiple runs
hist(gen2thresh_arr)

% Calculate mean, median, and standard deviation
gen2thresh_mean = mean(gen2thresh_arr);
gen2thresh_median = median(gen2thresh_arr);
gen2thresh_std = std(gen2thresh_arr);




