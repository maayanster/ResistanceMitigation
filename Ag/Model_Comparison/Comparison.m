%% DESCRIPTION
% Compares deterministic and stochastic model to find difference which is
% assumed to be the abundance effect for a range of values.
%% INPUTS - NE present only in refuge
Num_sim = 51;                                                                         % Number of simulations that comparison is run
q_freq_arr = [0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005];  % Initial frequency of resistant alleles
K_arr = [42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000];       % Carrying capacity
Pref_arr = [0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3];                        % Proportion of area that is refuge
WErr_ref_arr = [0.5,0.49,0.48,0.47,0.46,0.45,0.44,0.43,0.42,0.41,0.40,0.39,0.38,0.37,0.36,0.35,0.34,0.33,0.32,0.31,0.30,0.29,0.28,0.27,0.26,0.25,0.24,0.23,0.22,0.21,0.20,0.19,0.18,0.17,0.16,0.15,0.14,0.13,0.12,0.11,0.10,0.09,0.08,0.07,0.06,0.05,0.04,0.03,0.02,0.01,0];                   % Fitness of RR in refuge with natural enemies
WErs_ref_arr = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5];                         % Fitness of RS in refuge with natural enemies
WEss_ref_arr = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5];                     % Fitness of SS in refuge with natural enemies
WErr_toxic_arr = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];                                      % Fitness of RR in toxic with natural enemies
WErs_toxic_arr = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];                                      % Fitness of RS in toxic with natural enemies
WEss_toxic_arr = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];                                      % Fitness of SS in toxic with natural enemies
gen_num_arr = [2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500];           % Number of generations models are run];      

%% INITIALIZE
gen2thresh_sto_median_arr = zeros(1,length(q_freq_arr));
gen2thresh_det_arr = zeros(1,length(q_freq_arr));
diff_s_d_arr = zeros(1,length(q_freq_arr));
runs_thresh_arr = zeros(1, length(q_freq_arr));
std_stoch_arr = zeros(1, length(q_freq_arr));

%% CALCULATIONS
% Using the parameters above call function that runs both models yielding 
% generations to threshold for the stochastic model (median) and 
% deterministic 
for mm = 1:Num_sim
    display(mm)
    [gen2thresh_sto_median, gen2thresh_det, gen2thresh_s_d,std_stochastic,runs_threshold]= ...
        Run_single_batch(q_freq_arr(mm), Pref_arr(mm), K_arr(mm), ...
        WErr_ref_arr(mm), WErs_ref_arr(mm), WEss_ref_arr(mm), ...
        WErr_toxic_arr(mm), WErs_toxic_arr(mm), WEss_toxic_arr(mm),...
        gen_num_arr(mm));
    gen2thresh_sto_median_arr(mm) = gen2thresh_sto_median;
    gen2thresh_det_arr(mm) = gen2thresh_det;
    diff_s_d_arr(mm) = gen2thresh_s_d;
    std_stoch_arr(mm) = std_stochastic;
    runs_thresh_arr(mm)= runs_threshold;
end

%% COMPARE RESULTS
% Assign parameter that is varied on x-axis
%x = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100];
x = [0,2,4,6,8,10,12,14,16,18, 20,22,24,26,28, 30,32,34,36,38, 40,42,44,46,48, 50,52,54,56,58, 60,62,64,66,68, 70,72,74,76,78, 80,82,84,86,88, 90,92,94,96,98, 100];

% Plot calculated stochastic over deterministic ratio of generations to 
% threshold (y_axis) and the parameter that is varied - x (x_axis)
figure

%%%% plot without error bars
%plot(x, gen2thresh_sto_median_arr, x, gen2thresh_det_arr)

%%%%% plot with error bars
%errorbar(x, gen2thresh_sto_median_arr,std_stoch_arr, 'r');
%hold on;
%plot(x, gen2thresh_det_arr)

%%%%% use scatterplot and trendline
scatter(x, gen2thresh_sto_median_arr);
hold on
plot(x, gen2thresh_det_arr, 'r')
lsline

%%%%% add labels and legend
xlabel('Differential predation (%)', 'FontSize', 12);
ylabel('Generations to threshold','FontSize', 12);
legend('Stochastic', 'Deterministic');

%%%%% bar chart of runs that don't cross the threshold
figure
bar(runs_thresh_arr);
xlabel('Runs below threshold')

%% SAVE RESULTS
filename = strcat('Compare_NEref_RR',num2str(WErr_ref_arr(1)), '_', num2str(WErr_ref_arr(1) - WErr_ref_arr(2)), '_', num2str(WErr_ref_arr(end)), '.mat');
save(filename, 'gen2thresh_sto_median_arr','gen2thresh_det_arr','std_stoch_arr','runs_thresh_arr', '-mat');

% %% INPUTS - NE present in both
% Num_sim_b = 10;                                                                         % Number of simulations that comparison is run
% q_freq_arr_b = [0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005];  % Initial frequency of resistant alleles
% K_arr_b = [42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000];       % Carrying capacity
% Pref_arr_b = [0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3];                        % Proportion of area that is refuge
% WErr_ref_arr_b = [0, 0.05, 0.1, 0.15, 0.20,0.25,0.30,0.35,0.40,0.45];                % Fitness of RR in refuge with natural enemies
% WErs_ref_arr_b = [0.5, 0.5, 0.5, 0.5, 0.5,0.5,0.5,0.5,0.5,0.5];                    % Fitness of RS in refuge with natural enemies
% WEss_ref_arr_b = [0.5, 0.5, 0.5, 0.5, 0.5,0.5,0.5,0.5,0.5,0.5];                      % Fitness of SS in refuge with natural enemies
% WErr_toxic_arr_b = [0, 0.05, 0.1, 0.15, 0.20,0.25,0.30,0.35,0.40,0.45];                                      % Fitness of RR in toxic with natural enemies
% WErs_toxic_arr_b = [0.5, 0.5, 0.5, 0.5, 0.5,0.5,0.5,0.5,0.5,0.5];                                      % Fitness of RS in toxic with natural enemies
% WEss_toxic_arr_b = [0.5, 0.5, 0.5, 0.5, 0.5,0.5,0.5,0.5,0.5,0.5];                                        % Fitness of SS in toxic with natural enemies
% gen_num_arr_b = [2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500];           % Number of generations models are run];      
% 
% %% INITIALIZE
% median_stochastic_b = zeros(1,length(q_freq_arr_b));
% gen2thresh_mean_stochastic_b = zeros(1,length(q_freq_arr_b));
% runs_below_thresh_array_b = zeros(1, length(q_freq_arr_b));
% 
% median_det50_b = zeros(1,length(q_freq_arr_b));
% gen2thresh_det_mean_b = zeros(1, length(q_freq_arr_b));
% 
% %% CALCULATIONS
% % Using the parameters above call function that runs both models yielding 
% % generations to threshold for the stochastic model (median) and 
% % deterministic 
% for mm = 1:Num_sim_b
%     display(mm)
%     [gen2thresh_sto_median_b, gen2thresh_det_b, gen2thresh_s_d_b]= ...
%         Run_single_batch(q_freq_arr_b(mm), Pref_arr_b(mm), K_arr_b(mm), ...
%         WErr_ref_arr_b(mm), WErs_ref_arr_b(mm), WEss_ref_arr_b(mm), ...
%         WErr_toxic_arr_b(mm), WErs_toxic_arr_b(mm), WEss_toxic_arr_b(mm),...
%         gen_num_arr_b(mm));
%     gen2thresh_sto_median_arr_b(mm) = gen2thresh_sto_median_b;
%     gen2thresh_det_arr_b(mm) = gen2thresh_det_b;
%     diff_s_d_arr_b(mm) = gen2thresh_s_d_b;
% %end
% 
% %% NORMALIZE
% %norm_stochastic_b = ((gen2thresh_sto_median_arr_b./gen2thresh_sto_median_arr_b(1,1))-1)*100;
% %norm_difference_b = ((diff_s_d_arr_b./diff_s_d_arr_b(1,1))-1)*100;
% %norm_deterministic_b = ((gen2thresh_det_arr_b./gen2thresh_det_arr_b(1,1))-1)*100;
% 
% %% COMPARE RESULTS
% % Assign parameter that is varied on x-axis
% %x_b = [0, 0.05, 0.1, 0.15, 0.20,0.25,0.30,0.35,0.40,0.45];
% 
% % Plot calculated stochastic over deterministic ratio of generations to 
% % threshold (y_axis) and the parameter that is varied - x (x_axis)
% %subplot(2,1,2)
% %plot(x_b, gen2thresh_sto_median_arr_b, x_b, gen2thresh_det_arr_b, x, diff_s_d_arr_b)
% %xlabel('Mortality (%)', 'FontSize', 12);
% %legend('Stochastic', 'Deterministic','Difference');
% 
% %% SAVE RESULTS
% %%%%filename = strcat('Compare_NEboth_fitnessRR',num2str(WErr_ref_arr_b(1)), '_', num2str(WErr_ref_arr_b(end)));
% %%%%save(filename, 'gen2thresh_sto_median_arr_b','gen2thresh_det_arr_b','diff_s_d_arr_b');