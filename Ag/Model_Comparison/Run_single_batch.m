function[gen2thresh_sto_median, gen2thresh_det,gen2thresh_s_d,std_stochastic] = ...
    Run_single_batch(q_freq, Pref, K, WErr_ref, WErs_ref, WEss_ref,WErr_toxic, ...
    WErs_toxic,WEss_toxic, gen_num)
%% Description: runs the stochastic model 100 times and pairs this with a
% deterministic run with the same parameters

%% INITIALIZE
gen2thresh_sto = zeros(1, 100);
gen2thresh_det = 0;

%% CALCULATIONS
% Runs 100 simulations of stochastic model
tic
for nn = 1:100
        [generations2thresh] = Stochastic(q_freq,Pref, K, WErr_ref,...
            WErs_ref, WEss_ref, WErr_toxic, WErs_toxic, WEss_toxic, ...
            gen_num);
        
        % Saves the number of generations to threshold for each simulation
        gen2thresh_sto(nn) = generations2thresh;
        
end 
toc
% Plot histogram of 100 runs of stochastic model to evaluate its
% distribution
%hist(gen2thresh_sto)

% Find the median number of generations to threshold of 100 runs of
% stochastic model
gen2thresh_sto_median = prctile(gen2thresh_sto, 50);

% Calculate standard deviation of 100 stochastic runs
std_stochastic = std(gen2thresh_sto);

% Call deterministic model for given parameter
gen2thresh_det = Deterministic(q_freq,...
        Pref, WErr_ref, WErs_ref, ...
        WEss_ref, WErr_toxic, WErs_toxic,...
        WEss_toxic, gen_num);
    
% Find difference between stochastic and deterministic model
diff_s_d = gen2thresh_sto - gen2thresh_det;

gen2thresh_s_d = prctile(diff_s_d, 50);
    






