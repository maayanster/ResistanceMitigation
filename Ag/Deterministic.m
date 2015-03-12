function [q_freq, q_fixation] = Deterministic(q_freq, Pref, WErr_ref, WErs_ref, WEss_ref, gen_num)

% Deterministic_4Feb simulates pesticide resistance starting with relatively large
% pop and large 
% ------------------------------------------------------
% function [q_freq, population] = Deterministic_20Feb(q_freq, Pref, Wxrr_ref, Wxrs_ref, Wxss_ref, gen_num)
% ------------------------------------------------------
% Description: This function is meant to include random pop fluctuations in
%              small initial frequencies of a mutant allele, with selection 
%              from various sources
% Input:       {q_freq} 
%              {Pref} proportion of land in refuge (no pesticide/bt)
%              {WErr_ref} fitness of RR in refuge w/o natural enemies
%              {WErs_ref} fitness of RS in refuge w/o natural enemies
%              {WEss_ref} fitness of SS in refuge w/o natural enemies
%              {gen_num} number of generations to model to run for
% Output:      {q_freq} frequency of resistant alleles
%              a plot
%
% Adrian Semmelink
% Classification: honours project
% Last revision date: 20-Feb-2015

%% INITIALIZE
%q_freq = 0.00001;                      % initial frequency of R allele
p_freq = 1 - q_freq;                  % initial frequency of S allele
Wxrr_bt = 1;                        % fitness of RR in field w/o natural enemies
Wxss_bt = 0;                          % fitness of SS in field w/o natural enemies
Wxrs_bt = 0;                          % fitness of RS in field w/o natural enemies
WErr_bt = 1;                          % fitness of RR in field due to natural enemies
WErs_bt = 1;                          % fitness of RS in field due to natural enemies
WEss_bt = 1;                          % fitness of SS in field due to natural enemies
%WErr_ref = 1;                         % fitness of RR in refuge due to natural enemies
%WErs_ref = 1;                         % fitness of RS in refuge due to natural enemies
%WEss_ref = 1;                         % fitness of SS in refuge due to natural enemies
Wxrr_ref = 0.207;                      % fitness of RR in refuge w/o natural enemies
Wxrs_ref = 0.207;                      % fitness of RS in refuge w/o natural enemies
Wxss_ref = 0.207;                      % fitness of SS in refuge w/o natural enemies

RR = 0;                               % RR genotype frequency total              
RS = 0;                               % RS genotype frequency total
SS = 0;                               % SS genotype frequency total
RRref = 0;                            % RR genotype frequency refuge             
RSref = 0;                            % RS genotype frequency refuge
SSref = 0;                            % SS genotype frequency refuge
RRbt = 0;                             % RR genotype frequency Bt               
RSbt = 0;                             % RS genotype frequency Bt
SSbt = 0;                             % SS genotype frequency Bt
%Pref = 0.1;                           % Proportion of total field that is refuge
Pbt = 1-Pref;                         % Proportion of total field that is Bt
q_fixation = 0;

i = 0;                                % initialize generation count
%% OUTPUT
p_array=[];
q_array=[];

%% CALCULATIONS

while i <= gen_num;       %<= 100000                  % while frequency of R allele less than 0.5 (dominant resistance) calculate the following
    i = i+1;                      % Count generations
    
    %Mating Hardy-Weinberg - assume random mating - take the proportions of the allele and output
    %frequencies of three genotypes
    RR = q_freq^2;
    SS = p_freq^2;
    RS = 2*q_freq*p_freq;
    
    %Oviposition assume random
    RRref = RR*Pref;
    SSref = SS*Pref;
    RSref = RS*Pref;
    RRbt = RR*Pbt;
    SSbt = SS*Pbt;
    RSbt = RS*Pbt;
    
    %Selection - fitness w/o natural enemies
    RRref = RRref*Wxrr_ref;
    SSref = SSref*Wxss_ref;
    RSref = RSref*Wxrs_ref;
    RRbt = RRbt*Wxrr_bt;
    SSbt = SSbt*Wxss_bt;
    RSbt = RSbt*Wxrs_bt;
    
    %Selection with natural enemies
    RRref = RRref*WErr_ref;
    SSref = SSref*WEss_ref;
    RSref = RSref*WErs_ref;
    RRbt = RRbt*WErr_bt;
    SSbt = SSbt*WEss_bt;
    RSbt = RSbt*WErs_bt;
    
    % Finding p and q 
    q_sum = 2*RRref + RSref + 2*RRbt + RSbt;
    p_sum = 2*SSref + RSref + 2*SSbt + RSbt;
    
    q_freq = q_sum/(q_sum+p_sum);
    p_freq = p_sum/(p_sum+q_sum);
    
    % displaying generation at which fixation occurs
    if  q_freq >= 0.1 && q_fixation == 0;
        q_fixation = i;
    end
    
    
    %saving our fequencies
    q_array = [q_array, q_freq];
    p_array = [p_array, p_freq];
    
end

display(q_fixation)
gens = 1:length(q_array);
figure
plot(gens, q_array, gens, p_array)
xlabel('generations');
ylabel('frequency');
legend('q', 'p');

%display('i', 'output');