Num_sim = 11;                                                                         % Number of simulations that comparison is run
q_freq_arr = [0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005,0.005];  % Initial frequency of resistant alleles
K_arr = [42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000,42000];       % Carrying capacity
Pref_arr = [0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3,0.3];                        % Proportion of area that is refuge
WErr_ref_arr = [0.5, 0.45, 0.40, 0.35, 0.30,0.25,0.20,0.15,0.10,0.05,0];                   % Fitness of RR in refuge with natural enemies
WErs_ref_arr = [0.5, 0.45, 0.40, 0.35, 0.30,0.25,0.20,0.15,0.10,0.05,0];                         % Fitness of RS in refuge with natural enemies
WEss_ref_arr = [0.5, 0.45, 0.40, 0.35, 0.30,0.25,0.20,0.15,0.10,0.05,0];                         % Fitness of SS in refuge with natural enemies
WErr_toxic_arr = [1,1,1,1,1,1,1,1,1,1,1];                                      % Fitness of RR in toxic due to natural enemies 
WErs_toxic_arr = [1,1,1,1,1,1,1,1,1,1,1];                                      % Fitness of RS in toxic due to natural enemies
WEss_toxic_arr = [1,1,1,1,1,1,1,1,1,1,1];                                      % Fitness of SS in toxic due to natural enemies
gen_num_arr = [2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500,2500];           % Number of generations models are run];      
num_sim_arr = [100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100];


medians = [];
std_devs = [];
runs_bel_thresh = [];
med_times_pop_crashed = [];
stddevs_pop = [];

for ii=1:length(num_sim_arr)
    [a, b, c] = single_batch_stochastic(q_freq_arr(ii), Pref_arr(ii), K_arr(ii), WErr_ref_arr(ii), WErs_ref_arr(ii), WEss_ref_arr(ii), WErr_toxic_arr(ii), WErs_toxic_arr(ii), ...
    WEss_toxic_arr(ii), gen_num_arr(ii), num_sim_arr(ii));
    
    median_gen2thresh = median(a);
    stddev_gen2thresh = std(a);
    
    medians = [medians, median_gen2thresh];
    std_devs = [std_devs, stddev_gen2thresh];
    runs_bel_thresh = [runs_bel_thresh, b];
    
    median_pop_crashed = median(c);
    stddev_pop_crashed = std(c);
    
    med_times_pop_crashed = [med_times_pop_crashed, median_pop_crashed];
    stddevs_pop = [stddevs_pop, stddev_pop_crashed];
end

xrange = [0.5, 0.45, 0.40, 0.35, 0.30,0.25,0.20,0.15,0.10,0.05,0];
errorbar(xrange, medians, std_devs)

figure
errorbar(xrange, med_times_pop_crashed, stddevs_pop)

figure
plot(xrange, runs_bel_thresh)
