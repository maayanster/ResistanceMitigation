function[gen2thresh_arr, runs_below_thresh, times_crashed_arr] = Stochastic_Runs...
    (q_freq, Pref, K, WErr_ref, WErs_ref, WEss_ref,WErr_toxic, WErs_toxic, ...
    WEss_toxic, gen_num, num_sim)
% Stochastic_Runs runs stochastic insecticide resistance model 
% -------------------------------------------------------------------
% [gen2thresh_arr, runs_below_thresh, times_crashed_arr] = Stochastic_Runs...
%   (q_freq, Pref, K, WErr_ref, WErs_ref, WEss_ref,WErr_toxic, WErs_toxic, ...
%   WEss_toxic, gen_num, num_sim)
% -------------------------------------------------------------------
% Description:   Runs 100 simulations of stochastic insecticide resistance
%                model generating 100 values of generations to threshold,
%                which are sorted. Calculates and outputs the 25th 
%                percentile, median, 75th percentile and the mean
% Input:         {q_freq} Initial frequency of resistant alleles
%                {Pref} Proportion of area that is refuge (no Bt)
%                {WErr_ref} Fitness of RR in refuge with natural enemies
%                {WErs_ref} Fitness of RS in refuge with natural enemies
%                {WEss_ref} Fitness of SS in refuge with natural enemies
%                {gen_num} Number of generation that model runs
% Output:        {gen2thresh_25percentile} 25th percentile number of 
%                generations to resistance threshold
%                {gen2thresh_50percentile} median number of generations to 
%                resistance threshold
%                {gen2thresh_75percentile} 75th percentile number of 
%                generations to resistance threshold
%                {gen2thresh_mean} mean number of generations to 
%                resistance threshold
% Adrian Semmelink
% Classification: Honours project
% Last revision date: 15-April-2015

%% INITIALIZE
gen2thresh_arr = zeros(1, num_sim);
q_cell = {};
final_qs=[];
times_crashed_arr = [];

%% CALCULATIONS
% Runs 100 simulations of stochastic model
for nn = 1:num_sim
        [generations2thresh, q_arr] = Stochastic(q_freq,Pref, K, WErr_ref,...
            WErs_ref, WEss_ref, WErr_toxic, WErs_toxic, WEss_toxic, ...
            gen_num);
        
        % Saves the number of generations to threshold for each simulation
        gen2thresh_arr(nn) = generations2thresh;
        
       % Creates a matrix with all R allele frequencies across the 
       % maximum number of generations for all 100 simulations
        q_cell{nn} = q_arr;
        final_qs(nn) = q_arr(end);
        
        % Calculate times the R pop'n crashed to nothing per simulation
        repeats = diff(q_arr);
        q_log = logical([1, repeats]);
        q_noreps = q_arr(q_log);
        sum_zeros = length(q_noreps(q_noreps==0));
        times_crashed_arr = [times_crashed_arr, sum_zeros];       
end

fprintf('\nNumber of simulation in batch: %d\n',num_sim);
% Calculate number of runs (out of 100) that do not reach resistance
% threshold
runs_below_thresh = sum(isnan(gen2thresh_arr)); 
fprintf('Runs below threshold: %d\n', runs_below_thresh);

% Remove NaNs from gen2thresh
gen2thresh_NoNaNs = gen2thresh_arr(isfinite(gen2thresh_arr)); 

% Calculates 25th, median, 75th, percentile, and mean for number of generations to resistance threshold
gen2thresh_25percentile = prctile(gen2thresh_NoNaNs, 25); 
gen2thresh_50percentile = prctile(gen2thresh_NoNaNs, 50); 
gen2thresh_75percentile = prctile(gen2thresh_NoNaNs, 75); 
gen2thresh_st_dev = std(gen2thresh_NoNaNs);
gen2thresh_mean = mean(gen2thresh_NoNaNs); 

%% Display outputs - not used in comparison runs
fprintf('Median generations to threshold: %d\n25th and 75th percentile generations to threshold: %d, %d\nStandard dev generations to threshold: %d\nMean generations to threshold %d\n\n', ...
    gen2thresh_50percentile,gen2thresh_25percentile,gen2thresh_75percentile,gen2thresh_st_dev,gen2thresh_mean);

% Make histograms for final R frequency, generations to threshold and
% number of times q crashed to 0 per run

figure
hist(final_qs)
xlabel(sprintf('Proportion of resistant alleles at %d generations', gen_num));
ylabel('number of simulations');

figure
hist(times_crashed_arr)
xlabel('times R alleles crash to 0 per sim')
ylabel('number of simulations');

figure
hist(gen2thresh_arr)
xlabel ('generations to threshold')
ylabel ('number of simulations')
                