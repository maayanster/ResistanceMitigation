function [gen2thresh] = Deterministic(q_freq, Pref, ...
    WErr_ref, WErs_ref, WEss_ref, WErr_toxic, WErs_toxic, WEss_toxic, ...
    gen_num)
% Deterministic_13March genetic model simulating insecticide resistance 
% ------------------------------------------------------
% function [q_freq, q_fixation] = Deterministic_13March(q_freq, Pref, ...
%                                 WErr_ref, WErs_ref,WErr_toxic, ...
%                                 WErs_toxic, WEss_toxic, WEss_ref,...
%                                 gen_num)
% ------------------------------------------------------
% Description: Model simulating pest inseciticide resistance assuming an 
%              infinite pest population, Hardy-weinberg, with selection
%              from natural enemies and other sources
% Input:       {q_freq} 
%              {Pref} Proportion of land in refuge (no pesticide/bt)
%              {WErr_ref} Fitness of RR in refuge with natural enemies
%              {WErs_ref} Fitness of RS in refuge with natural enemies
%              {WEss_ref} Fitness of SS in refuge with natural enemies
%              {WErr_toxic} Fitness of RR in refuge with natural enemies
%              {WErs_toxic} Fitness of RS in refuge with natural enemies
%              {WEss_toxic} Fitness of SS in refuge with natural enemies
%              {gen_num} Number of generations to model to run for
% Output:      {q_freq} Frequency of resistant alleles (not used in all
%              runs)
%              {gen2thresh} Number of generations to resistance threshold
% Adrian Semmelink
% Classification: Honours project
% Last revision date: 20-Feb-2015

%% INITIALIZE
p_freq = 1 - q_freq;                  % Initial frequency of S allele
Wxrr_toxic = 0.207;                   % Fitness of RR in field w/o natural enemies
Wxss_toxic = 0;                       % Fitness of SS in field w/o natural enemies
Wxrs_toxic = 0;                       % Fitness of RS in field w/o natural enemies
Wxrr_ref = 0.207;                     % Fitness of RR in refuge w/o natural enemies
Wxrs_ref = 0.207;                     % Fitness of RS in refuge w/o natural enemies
Wxss_ref = 0.208;                     % Fitness of SS in refuge w/o natural enemies
Ptoxic = 1-Pref;                      % Proportion of total field that is toxic
gen2thresh = 0;                       % Generation when q_freq is larger than 0.1
MutationR = 0.00005;                  % Mutation rate of R to S or S to R 
i = 0;                                % Initialize generation count
q_threshold = 0.1;                    % q frequency at which program stops

%% OUTPUT
p_array=[];
q_array=[];

%% CALCULATIONS
while i <= gen_num;         % while generations less than gen_num run loop      
    % Count generations
    i = i+1;                      
    
    % Apply mutation rate 
    q_freq = q_freq + p_freq*MutationR;
    p_freq = p_freq + q_freq*MutationR;
    
    % Hardy-weinberg ratios applied to find expected relative frequency  
    RR = q_freq^2;
    SS = p_freq^2;
    RS = 2*q_freq*p_freq;
    
    % Oviposition distribution across field according to proportion
    % of toxic and refuge zones
    RRref = RR*Pref;
    SSref = SS*Pref;
    RSref = RS*Pref;
    RRtoxic = RR*Ptoxic;
    SStoxic = SS*Ptoxic;
    RStoxic = RS*Ptoxic;
    
    % Selection - fitness without natural enemies
    RRref = RRref*Wxrr_ref;
    SSref = SSref*Wxss_ref;
    RSref = RSref*Wxrs_ref;
    RRtoxic = RRtoxic*Wxrr_toxic;
    SStoxic = SStoxic*Wxss_toxic;
    RStoxic = RStoxic*Wxrs_toxic;
    
    % Selection by natural enemies
    RRref = RRref*WErr_ref;
    SSref = SSref*WEss_ref;
    RSref = RSref*WErs_ref;
    RRtoxic = RRtoxic*WErr_toxic;
    SStoxic = SStoxic*WEss_toxic;
    RStoxic = RStoxic*WErs_toxic;
    
    % Finding p and q 
    q_sum = 2*RRref + RSref + 2*RRtoxic + RStoxic;
    p_sum = 2*SSref + RSref + 2*SStoxic + RStoxic;
    
    q_freq = q_sum/(q_sum+p_sum);
    p_freq = p_sum/(p_sum+q_sum);
    
    % Displaying generation at which fixation occurs
    if  q_freq >= q_threshold && gen2thresh == 0;
        gen2thresh = i;
        break
    end
    
    % Saving fequencies for plots (not used in all simulations)
    q_array = [q_array, q_freq];
    p_array = [p_array, p_freq];
    
end

% If no resistance developed by the end of the run the number of
% generations to threshold is recorded as NaN
if q_freq <= q_threshold;
    gen2thresh = NaN;
end

%% Display outputs - not used in comparison runs
% gens = 1:length(q_array);
% figure
% plot(gens, q_array, gens, p_array)
% xlabel('generations', 'FontSize', 12);
% ylabel('frequency', 'FontSize', 12);
%legend('q', 'p');