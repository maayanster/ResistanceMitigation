%% DESCRIPTION
% Compares deterministic and stochastic model to find difference which is
% assumed to be the abundance effect for a range of values.
%% INPUTS - NE present only in refuge
Num_sim = 10;                                                                         % Number of simulations that comparison is run
q_freq_arr = [0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005];  % Initial frequency of resistant alleles
K_arr = [42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000];       % Carrying capacity
Pref_arr = [0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3];                        % Proportion of area that is refuge
WErr_ref_arr = [0, 0.05, 0.1, 0.15, 0.20,0.25,0.30,0.35,0.40,0.45];                % Fitness of RR in refuge with natural enemies
WErs_ref_arr = [0.5, 0.5, 0.5, 0.5, 0.5,0.5,0.5,0.5,0.5,0.5];                    % Fitness of RS in refuge with natural enemies
WEss_ref_arr = [0.5, 0.5, 0.5, 0.5, 0.5,0.5,0.5,0.5,0.5,0.5];                      % Fitness of SS in refuge with natural enemies
WErr_toxic_arr = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1];                                      % Fitness of RR in toxic with natural enemies
WErs_toxic_arr = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1];                                      % Fitness of RS in toxic with natural enemies
WEss_toxic_arr = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1];                                        % Fitness of SS in toxic with natural enemies
gen_num_arr = [2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500];           % Number of generations models are run];      

%% INITIALIZE
median_stochastic = zeros(1,length(q_freq_arr));
gen2thresh_mean_stochastic = zeros(1,length(q_freq_arr));
runs_below_thresh_array = zeros(1, length(q_freq_arr));

median_det50 = zeros(1,length(q_freq_arr));
gen2thresh_det_mean = zeros(1, length(q_freq_arr));

%% CALCULATIONS
% Using the parameters above call function that runs both models yielding 
% generations to threshold for the stochastic model (median) and 
% deterministic 
for mm = 1:Num_sim
    display(mm)
    [gen2thresh_sto_median, gen2thresh_det, gen2thresh_s_d]= ...
        Runs(q_freq_arr(mm), Pref_arr(mm), K_arr(mm), ...
        WErr_ref_arr(mm), WErs_ref_arr(mm), WEss_ref_arr(mm), ...
        WErr_toxic_arr(mm), WErs_toxic_arr(mm), WEss_toxic_arr(mm),...
        gen_num_arr(mm));
    gen2thresh_sto_median_arr(mm) = gen2thresh_sto_median;
    gen2thresh_det_arr(mm) = gen2thresh_det;
    diff_s_d_arr(mm) = gen2thresh_s_d;
end

%% NORMALIZE
%norm_stochastic = ((gen2thresh_sto_median_arr./gen2thresh_sto_median_arr(1,1))-1)*100;
%norm_difference = ((diff_s_d_arr./diff_s_d_arr(1,1))-1)*100;
%norm_deterministic = ((gen2thresh_det_arr./gen2thresh_det_arr(1,1))-1)*100;


%% COMPARE RESULTS
% Assign parameter that is varied on x-axis
x = [0, 0.025, 0.05, 0.075,0.1,0.125, 0.15, 0.175, 0.2, 0.225, 0.25,0.275,0.3,0.325,0.35,0.375,0.4,0.425,0.45,0.475,0.5];

% Plot calculated stochastic over deterministic ratio of generations to 
% threshold (y_axis) and the parameter that is varied - x (x_axis)
figure
subplot(2,1,1)
plot(x, gen2thresh_sto_median_arr, x, gen2thresh_det_arr)
ylabel('Generations to threshold','FontSize', 12);
legend('Stochastic', 'Deterministic');

%% INPUTS - NE present in both
Num_sim_b = 10;                                                                         % Number of simulations that comparison is run
q_freq_arr_b = [0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005];  % Initial frequency of resistant alleles
K_arr_b = [42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000, 42000];       % Carrying capacity
Pref_arr_b = [0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3];                        % Proportion of area that is refuge
WErr_ref_arr_b = [0, 0.05, 0.1, 0.15, 0.20,0.25,0.30,0.35,0.40,0.45];                % Fitness of RR in refuge with natural enemies
WErs_ref_arr_b = [0.5, 0.5, 0.5, 0.5, 0.5,0.5,0.5,0.5,0.5,0.5];                    % Fitness of RS in refuge with natural enemies
WEss_ref_arr_b = [0.5, 0.5, 0.5, 0.5, 0.5,0.5,0.5,0.5,0.5,0.5];                      % Fitness of SS in refuge with natural enemies
WErr_toxic_arr_b = [0, 0.05, 0.1, 0.15, 0.20,0.25,0.30,0.35,0.40,0.45];                                      % Fitness of RR in toxic with natural enemies
WErs_toxic_arr_b = [0.5, 0.5, 0.5, 0.5, 0.5,0.5,0.5,0.5,0.5,0.5];                                      % Fitness of RS in toxic with natural enemies
WEss_toxic_arr_b = [0.5, 0.5, 0.5, 0.5, 0.5,0.5,0.5,0.5,0.5,0.5];                                        % Fitness of SS in toxic with natural enemies
gen_num_arr_b = [2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500];           % Number of generations models are run];      

%% INITIALIZE
median_stochastic_b = zeros(1,length(q_freq_arr_b));
gen2thresh_mean_stochastic_b = zeros(1,length(q_freq_arr_b));
runs_below_thresh_array_b = zeros(1, length(q_freq_arr_b));

median_det50_b = zeros(1,length(q_freq_arr_b));
gen2thresh_det_mean_b = zeros(1, length(q_freq_arr_b));

%% CALCULATIONS
% Using the parameters above call function that runs both models yielding 
% generations to threshold for the stochastic model (median) and 
% deterministic 
for mm = 1:Num_sim_b
    display(mm)
    [gen2thresh_sto_median_b, gen2thresh_det_b, gen2thresh_s_d_b]= ...
        Runs(q_freq_arr_b(mm), Pref_arr_b(mm), K_arr_b(mm), ...
        WErr_ref_arr_b(mm), WErs_ref_arr_b(mm), WEss_ref_arr_b(mm), ...
        WErr_toxic_arr_b(mm), WErs_toxic_arr_b(mm), WEss_toxic_arr_b(mm),...
        gen_num_arr_b(mm));
    gen2thresh_sto_median_arr_b(mm) = gen2thresh_sto_median_b;
    gen2thresh_det_arr_b(mm) = gen2thresh_det_b;
    diff_s_d_arr_b(mm) = gen2thresh_s_d_b;
end

%% NORMALIZE
%norm_stochastic_b = ((gen2thresh_sto_median_arr_b./gen2thresh_sto_median_arr_b(1,1))-1)*100;
%norm_difference_b = ((diff_s_d_arr_b./diff_s_d_arr_b(1,1))-1)*100;
%norm_deterministic_b = ((gen2thresh_det_arr_b./gen2thresh_det_arr_b(1,1))-1)*100;

%% COMPARE RESULTS
% Assign parameter that is varied on x-axis
x_b = [50, 55, 60, 65, 70, 75, 80, 85, 90, 95];

% Plot calculated stochastic over deterministic ratio of generations to 
% threshold (y_axis) and the parameter that is varied - x (x_axis)
subplot(2,1,2)
plot(x_b, gen2thresh_sto_median_arr_b, x_b, gen2thresh_det_arr_b)
xlabel('Mortality (%)', 'FontSize', 12);
legend('stochastic', 'deterministic');

