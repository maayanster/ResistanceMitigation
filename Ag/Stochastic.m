function [q_freq, population] = stochastic_model_18Feb(population, Pref, Wxrr_ref, Wxrs_ref, Wxss_ref, gen_num)

% STOCHASTIC_MODEL simulates pesticide resistance starting with small pop
% ------------------------------------------------------
% [q_freq, population] = stochastic_model(population, Pref, Wxrr_ref, Wxrs_ref, Wxss_ref, gen_num)
% ------------------------------------------------------
% Description: This function is meant to include random pop fluctuations in
%              small initial frequencies of a mutant allele, with selection 
%              from various sources
% Input:       {population} initial pop size (this will vary in the model)
%              {Pref} proportion of land in refuge (no pesticide/bt)
%              {Wxrr_ref} fitness of RR in refuge w/o natural enemies
%              {Wxrs_ref} fitness of RS in refuge w/o natural enemies
%              {Wxss_ref} fitness of SS in refuge w/o natural enemies
%              {gen_num} numbe rof generation to run model for
% Output:      {q_freq} frequency of resistant alleles
%              {population} total pop soze at end of simulation
%              a plot
%
% Adrian Semmelink
% Classification: honours project
% Last revision date: 18-Feb-2015

%% Initialize

%%% Function inputs don't need to initialize
% Population = 100000;                % Initial pest population - not used
% Wxrr_ref = 0.79;                    % fitness of RR in refuge w/o natural enemies
% Wxrs_ref = 0.79;                    % fitness of RS in refuge w/o natural enemies
% Wxss_ref = 0.79;                    % fitness of SS in refuge w/o natural enemies
% Pref = 0.4;                         % Proportion of total field that is refuge

q_freq = 10/population;                % initial frequency of R allele, starting with 1 individual
K = 420000;                          % Carrying capacity
p_freq = 1 - q_freq;                  % initial frequency of S allele
mean_no_progeny = 70;                 % mean progeny produced from one mating (fecundity)
std_dev_progeny = 1;                  % standard dev of progeny distribution ---> should be standard error
Wxrr_bt = 0.207;                          % fitness of RR in field w/o natural enemies
Wxss_bt = 0;                          % fitness of SS in field w/o natural enemies
Wxrs_bt = 0;                          % fitness of RS in field w/o natural enemies
WErr_bt = 1;                          % fitness of RR in field due to natural enemies
WErs_bt = 1;                          % fitness of RS in field due to natural enemies
WEss_bt = 1;                          % fitness of SS in field due to natural enemies
WErr_ref = 1;                         % fitness of RR in refuge due to natural enemies
WErs_ref = 1;                         % fitness of RS in refuge due to natural enemies
WEss_ref = 1;                         % fitness of SS in refuge due to natural enemies
RR = 0;                               % RR genotype frequency total              
RS = 0;                               % RS genotype frequency total
SS = 0;                               % SS genotype frequency total
RRref = 0;                            % RR genotype frequency refuge             
RSref = 0;                            % RS genotype frequency refuge
SSref = 0;                            % SS genotype frequency refuge
RRbt = 0;                             % RR genotype frequency Bt               
RSbt = 0;                             % RS genotype frequency Bt
SSbt = 0;                             % SS genotype frequency Bt
Pbt = 1-Pref;                         % Proportion of total field that is Bt

i = 0;                                % Initialize generation count

%% OUTPUT
p_array=[];
q_array=[];
pop_array=[];

%% CALCULATIONS

while i <= gen_num %q_freq <= 0.1 && i <= 20;   % while frequency of R allele less than 0.1 or number of generations less than 20 calculate loop
    i = i+1;                      % Count generations
    disp(i)
    % Round population each loop to prevent binornd returning NaN
    population = round(population);
    
    % Determine how many males/females as integers - use binomial - assume
    % replacement 
    % distribution to choose a number of females, assume remainder males)  
    Fem = binornd(population,0.5,1);
    Male = population - Fem;   
    
    % Convert fem / male pop to type double - poissrnd doesn't work with
    % integers
    %Fem = double(Fem);
    %Male = double(Male);
               
    % Hardy-weinberg ratios are used to calculate expected relative frequency
    % Is it necessary to insert stochasticity into this step?  
    SS = p_freq^2;
    RS = 2*q_freq*p_freq;
    RR = q_freq^2;
    
    %%% If error statement to check that RR/RS/SS total 1
    if ((RR+RS+SS)-1) > 0.000001 && ((RR+RS+SS)-1) < -0.000001; 
        error('RR+RS+SS do not equal 1');
    end

    % Sample the number of individuals from each genotype in this generation
    % using poisson distribution
    
    MaleRR = poissrnd(RR*Male, 1);    
    FemRR = poissrnd(RR*Fem, 1);
    MaleRS = poissrnd(RS*Male, 1);
    FemRS = poissrnd(RS*Fem, 1);
    
    % Calculate remainders for SS genotypes
    MaleSS = Male - (MaleRS + MaleRR);
    FemSS = Fem - (FemRS + FemRR);   
    if (MaleRS + MaleRR) > Male;
        MaleSS = 0;
    end
    if (FemRS + FemRR) > Fem
        FemSS = 0;
    end
    
    %%% Calculate pairings between genotypes
    % Tracking Female RR genotype
    FemRRxMaleRR = poissrnd(FemRR*MaleRR/Male, 1); 
    FemRRxMaleRS = poissrnd(FemRR*MaleRS/Male, 1);
    FemRRxMaleSS = FemRR - (FemRRxMaleRR + FemRRxMaleRS);   
    if FemRRxMaleSS < 0 
        FemRRxMaleSS = 0;
        FemRRxMaleRS = FemRRxMaleRS/(FemRRxMaleRR+FemRRxMaleRS)*FemRR;
        FemRRxMaleRR = FemRRxMaleRR/(FemRRxMaleRR+FemRRxMaleRS)*FemRR;
    end
        
    
    % Tracking Female RS genotype
    FemRSxMaleRR = poissrnd(FemRS*MaleRR/Male, 1);
    FemRSxMaleRS = poissrnd(FemRS*MaleRS/Male, 1);
    FemRSxMaleSS = FemRS - (FemRSxMaleRR + FemRSxMaleRS);
    if FemRSxMaleSS < 0
        FemRSxMaleSS = 0;
        FemRSxMaleRR = FemRSxMaleRR/(FemRSxMaleRR + FemRSxMaleRS)*FemRS;
        FemRSxMaleRS = FemRSxMaleRS/(FemRSxMaleRR + FemRSxMaleRS)*FemRS;
    end
    
    % Tracking Female SS genotype
    FemSSxMaleRR = poissrnd(FemRS*MaleRR/Male, 1);
    FemSSxMaleRS = poissrnd(FemRS*MaleRS/Male, 1);
    FemSSxMaleSS = FemSS - (FemSSxMaleRR + FemSSxMaleRS);
    if FemSSxMaleSS < 0;
        FemSSxMaleRR = FemSSxMaleRR/(FemSSxMaleRR+FemSSxMaleRS)*FemSS;
        FemSSxMaleRS = FemSSxMaleRS/(FemSSxMaleRR+FemSSxMaleRS)*FemSS;
        FemSSxMaleSS = 0;
    end
   
    
    % Calculate 6 different types of pairings:
    RRxRR = FemRRxMaleRR; 
    RRxRS = FemRRxMaleRS + FemRSxMaleRR;
    RRxSS = FemRRxMaleSS + FemSSxMaleRR;
    SSxSS = FemSSxMaleSS;
    SSxRS = FemSSxMaleRS + FemRSxMaleSS;
    RSxRS = FemRSxMaleRS;
    
    % For each type of pairing find progeny number according to normal distribution 
    % still need to change std_dev to std error 
    RRxRR_progeny = round(RRxRR*normrnd(mean_no_progeny, std_dev_progeny, 1));
    RRxRS_progeny = round(RRxRS*normrnd(mean_no_progeny, std_dev_progeny, 1));
    RRxSS_progeny = round(RRxSS*normrnd(mean_no_progeny, std_dev_progeny, 1));
    SSxSS_progeny = round(SSxSS*normrnd(mean_no_progeny, std_dev_progeny, 1));
    SSxRS_progeny = round(SSxRS*normrnd(mean_no_progeny, std_dev_progeny, 1));
    RSxRS_progeny = round(RSxRS*normrnd(mean_no_progeny, std_dev_progeny, 1));
    
    % For each type of pairing assigning number of progeny for each
    % genotype (birth rate per pairing)
    % Find number of RR progeny produced by RRxRR pairing
    RRxRR_progeny_RR = RRxRR_progeny;
    
    % Find number of RR and RS progeny produced by RRxRS pairing
    RRxRS_progeny_RR = binornd(RRxRS_progeny, 0.5);
    RRxRS_progeny_RS = RRxRS_progeny - RRxRS_progeny_RR;
    
    % Find number of RS progeny produced by RRxSS pairing
    RRxSS_progeny_RS = RRxSS_progeny;
    
    % Find number of SS progeny produced by SSxSS pairing
    SSxSS_progeny_SS = SSxSS_progeny;
    
    % Find number of SS and RS progeny produced by SSxRS pairing
    SSxRS_progeny_RS = binornd(SSxRS_progeny, 0.5);
    SSxRS_progeny_SS = SSxRS_progeny - SSxRS_progeny_RS;
    
    % Find number of SS, RS, and RR progeny produced by RSxRS pairing
    RSxRS_progeny_SS = binornd(RSxRS_progeny, 0.25);
    RSxRS_progeny_RR = binornd(RSxRS_progeny, 0.25);
    RSxRS_progeny_RS = RSxRS_progeny - (RSxRS_progeny_SS+RSxRS_progeny_RR);
    
    % Find total amount of progeny for RR, SS, and RS
    tot_progeny_RR = RRxRR_progeny_RR + RRxRS_progeny_RR + RSxRS_progeny_RR;
    tot_progeny_RS = RRxRS_progeny_RS + RRxSS_progeny_RS + SSxRS_progeny_RS + RSxRS_progeny_RS;
    tot_progeny_SS = SSxSS_progeny_SS + SSxRS_progeny_SS + RSxRS_progeny_SS;
    
    %Oviposition assume random 
    RRref = tot_progeny_RR*Pref;
    SSref = tot_progeny_SS*Pref;
    RSref = tot_progeny_RS*Pref;
    RRbt = tot_progeny_RR*Pbt;
    SSbt = tot_progeny_SS*Pbt;
    RSbt = tot_progeny_RS*Pbt;   
    
    %Selection - fitness w/o natural enemies 
    RRref = RRref*Wxrr_ref;
    SSref = SSref*Wxss_ref;
    RSref = RSref*Wxrs_ref;
    RRbt = RRbt*Wxrr_bt;
    SSbt = SSbt*Wxss_bt;
    RSbt = RSbt*Wxrs_bt;
    
    %Selection due to natural enemies 
    RRref = RRref*WErr_ref;
    SSref = SSref*WEss_ref;
    RSref = RSref*WErs_ref;
    RRbt = RRbt*WErr_bt;
    SSbt = SSbt*WEss_bt;
    RSbt = RSbt*WErs_bt;
    
    % Selection on all genotypes, (carrying capacity situation to create a 
    % roughly logistic growth curve)
    %pop_total = RRref + SSref + RSref + RRbt + SSbt + RSbt;
    
    Growth_rate = (RRref + SSref + RSref + RRbt + SSbt + RSbt)/population;
    
    pop_new = population + Growth_rate*population*(1-population/K);
    
    if pop_new <= 0
        pop_new = K;
    end
   
    RRref_x = RRref/population*pop_new;
    SSref_x = SSref/population*pop_new;
    RSref_x = SSref/population*pop_new;
    RRbt_x = RRbt/population*pop_new;
    SSbt_x = SSbt/population*pop_new;
    RSbt_x = SSbt/population*pop_new;
    
    
    %Normalize to new population after carrying capacity change in
    %population
    %RRref = RRref/pop_total*pop_new;
    %SSref = SSref/pop_total*pop_new;
    %RSref = SSref/pop_total*pop_new;
    %RRbt = RRbt/pop_total*pop_new;
    %SSbt = SSbt/pop_total*pop_new;
    %RSbt = SSbt/pop_total*pop_new;
    
    % Prevent overshoot of carrying capacity
    %if pop_total >= K;
    %   pop = pop_total*((K - pop_total)/K);
    % elseif pop_total < K
    %     pop = pop_total*(1 - pop_total/K); 
    % end
    
    %hold = pop_total*(1 - pop_total/K);
    %if hold < 0
    %    K_survival = 0.1;
    %else
    %K_survival = hold/pop_total;
    %end
    K_survival = 1;
   
    %find number of each allele after all selection stages
    q_sum = 2*RRref_x*K_survival + RSref_x*K_survival + 2*RRbt_x*K_survival + RSbt_x*K_survival;
    p_sum = 2*SSref_x*K_survival + RSref_x*K_survival + 2*SSbt_x*K_survival + RSbt_x*K_survival;
    
    % Finding p and q
    q_freq = round(q_sum/(q_sum+p_sum));
    p_freq = round(p_sum/(p_sum+q_sum));
    
    % displaying generation at which fixation occurs
    if  q_freq >= 0.1 && q_fixation ==0;
        q_fixation = i;
    end
    
    population = (round(p_sum+q_sum))/2;
    
    %saving our population and fequencies info
    q_array = [q_array, q_freq];
    p_array = [p_array, p_freq];
    pop_array = [pop_array, population];
    
end

gens = 1:length(q_array);
%figure
subplot(2,1,1)
plot(gens, q_array, gens, p_array)
xlabel('generations');
ylabel('frequency');
legend('q', 'p');
subplot(2,1,2)
plot(gens, pop_array)
xlabel('generations');
ylabel('population');
display('i', 'output');
