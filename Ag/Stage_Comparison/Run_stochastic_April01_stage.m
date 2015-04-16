function[median_gen2thresh, median_times_crashed, median_initial_end, runs_below_thresh, median_spread_end] = Run_stochastic_April01_stage...
    (q_freq, Pref, K, WErr_ref, WErs_ref, WEss_ref,WErr_toxic, WErs_toxic, ...
    WEss_toxic, gen_num)
% Run_stochastic_12March runs stochastic insecticide resistance model 
% -------------------------------------------------------------------
% [q_fix_percentileA, q_fix_percentileB, q_fix_percentileC, q_fix_mean] = 
% Run_stochastic_12March(q_freq, Pref, WErr_ref, WErs_ref, WEss_ref, 
% gen_num)
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
% Last revision date: 13-March-2015

%% INITIALIZE
Initial_end_array = [];
Spread_end_array = [];
gen2thresh_array = [];
times_crashed_array = [];
gen2thresh_NoNaNs = 0;

%% CALCULATIONS
% Plot graph of 100 simulations - not used for comparison model
% figure 
% xlabel('generations', 'FontSize', 12);
% ylabel('frequency', 'FontSize', 12);
% legend('q', 'p');
% hold on

% Runs 100 simulations of stochastic model
for nn = 1:100
        [gen2thresh, q_array] = stochastic_15March(q_freq,Pref, K, WErr_ref,...
            WErs_ref, WEss_ref, WErr_toxic, WErs_toxic, WEss_toxic, ...
            gen_num);
        
        % Saves the number of generations to threshold for each simulation
        gen2thresh_array = [gen2thresh_array, gen2thresh];
               
        % Calculate times crashed per simulation
        repeats = diff(q_array);
        q_log = logical([1, repeats]);
        q_noreps = q_array(q_log);
        sum_zeros = length(q_noreps(q_noreps==0));
        times_crashed_array = [times_crashed_array, sum_zeros];
        
        % Calculate at what generation simulation stop crashing (end of
        % initial acquisition stage)
        if isnan(gen2thresh) == 0
        Initial_end = find(q_array==0,1,'last');
        Initial_end_array = [Initial_end_array,Initial_end];
        Spread_end = gen2thresh - Initial_end;
        Spread_end_array = [Spread_end_array, Spread_end];
        end
end

% Remove NaNs from gen2thresh
gen2thresh_NoNaNs = gen2thresh_array(isfinite(gen2thresh_array));

% Calculate median number of crashes per simulation
median_times_crashed = median(times_crashed_array);

% Calculate median generation at which q frequency reaches end of initial
% stage and the threshold
median_initial_end = median(Initial_end_array);
median_spread_end = median(Spread_end_array);
median_gen2thresh = median(gen2thresh_NoNaNs);

% Calculate number of runs (out of 100) that do not reach resistance
% threshold
runs_below_thresh = 100 - length(gen2thresh_NoNaNs); 

               