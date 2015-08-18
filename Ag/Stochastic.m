function [generations2thresh, q_array] = Stochastic(q_freq, Pref, K, ...
    WErr_ref, WErs_ref, WEss_ref, WErr_toxic, WErs_toxic, WEss_toxic,...
   gen_num)
% Stochastic simulates insecticide resistance starting 
% ------------------------------------------------------
% [generations2thresh, q_array] = Stochastic(q_freq, Pref, K, WErr_ref,...
%  WErs_ref, WEss_ref, WErr_toxic,WErs_toxic, WEss_toxic,gen_num)
% ------------------------------------------------------
% Description: Model insecticie resistance in a finite population, 
%              panmictic mating, with selection from various sources,
%              and density dependent effects
% Input:       {q_freq} initial frequency of resistant alleles
%              {Pref} proportion of area that is refuge (no Bt)
%              {K} carrying capacity
%              {WErr_ref} fitness of RR in refuge with natural enemies
%              {WErs_ref} fitness of RS in refuge with natural enemies
%              {WEss_ref} fitness of SS in refuge with natural enemies
%              {WErr_toxic} Fitness of RR in field due to natural enemies
%              {WErs_toxic} Fitness of RS in field due to natural enemies
%              {WEss_toxic} Fitness of SS in field due to natural enemies
%              {gen_num} number of generation that model runs
% Output:      {generations2thresh} Generation when q_freq is larger than 0.1
%              {q_array} array of resistant allele frequencies over number
%              of generations

% Adrian Semmelink
% Classification: Insecticide resistance project
% Last revision date: 11-July-2015

%% Initialize
population = 2900;                     % Initial pest population  
p_freq = 1 - q_freq;                  % Initial frequency of S allele
progeny = 100;                         % Number of progeny produced per female 
Wxrr_toxic = 0.207;                   % Fitness of RR in field w/o natural enemies
Wxss_toxic = 0;                       % Fitness of SS in field w/o natural enemies
Wxrs_toxic = 0;                       % Fitness of RS in field w/o natural enemies
Wxrr_ref = 0.207;                     % Fitness of RR in refuge w/o natural enemies
Wxrs_ref = 0.207;                     % Fitness of RS in refuge w/o natural enemies
Wxss_ref = 0.208;                     % Fitness of SS in refuge w/o natural enemies
MutationR = 0.00005;                  % Mutation rate of R to S or S to R (Sisterson, 2004)
i = 0;                                % Initialize generation count
q_threshold= 0.1;                     % q frequency at which program stops 
generations2thresh = 0;               % Initialize number of generations to threshold
intrinsicR = 0.15;                     % Intrinsic growth rate (Gryspeit, 2012)

%% OUTPUT
p_array=[];
q_array=[];
pop_array = [];

%% CALCULATIONS
% Calculate number of generations for loop
gen_number = gen_num - 1;

% Run stochastic model
while i <= gen_number 
    i = i+1;                      % Count generations
    %%% STEP 1: Calculate proportions of adults by genotype and sex
    % Round population each loop to convert into whole number
    population = round(population);

    % Mutation rate applied
    q_freq = q_freq + (p_freq*MutationR);
    p_freq = p_freq + (q_freq*MutationR);

    % Determine proportion of males and females by applying binomial 
    % distribution to choose number of females, assume remainder males 
    Fem = binornd(population,0.5,1);
    Male = population - Fem;
    
    % Hardy-weinberg ratios applied to find expected relative frequency  
    SS = p_freq^2;
    RS = 2*(q_freq*p_freq);
    RR = q_freq^2;

    % If error statement to check that RR/RS/SS total 1
    if ((RR+RS+SS)-1) > 0.000001 && ((RR+RS+SS)-1) < -0.000001; 
        error('RR+RS+SS do not equal 1');
    end

    % Sample the number of individuals from each genotype in this generation
    % using poisson distribution - sampling "small" genotypes first
    MaleRR = poissrnd(RR*Male, 1);    
    FemRR = poissrnd(RR*Fem, 1);
    MaleRS = poissrnd(RS*Male, 1);
    FemRS = poissrnd(RS*Fem, 1);

    % Calculate remainders for SS genotypes
    MaleSS = Male - (MaleRS + MaleRR);   
    if (MaleRS + MaleRR) > Male;
        MaleSS = 0;
    end

    FemSS = Fem - (FemRS + FemRR);
    if (FemRS + FemRR) > Fem
        FemSS = 0;
    end

    MaleRR = round(RR*Male);    
    FemRR = round(RR*Fem);
    MaleRS = round(RS*Male);
    FemRS = round(RS*Fem);
    MaleSS = round(SS*Male);
    FemSS = round(SS*Fem);

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
    FemSSxMaleRR = poissrnd(FemSS*MaleRR/Male, 1);
    FemSSxMaleRS = poissrnd(FemSS*MaleRS/Male, 1);
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

    %%% STEP 2: For each type of pairing assign random number of progeny
    % for each genotype (birth rate per pairing)
    
    % Find number of RR progeny produced by RRxRR pairing
    RRxRR_progeny_RR = round(RRxRR*progeny);
    RRxRR_progeny_RR_random = poissrnd(RRxRR_progeny_RR);
    
    % Find number of RS progeny produced by RRxSS pairing
    RRxSS_progeny_RS = round(RRxSS*progeny);
    RRxSS_progeny_RS_random = poissrnd(RRxSS_progeny_RS);

    % Find number of SS progeny produced by SSxSS pairing
    SSxSS_progeny_SS = round(SSxSS*progeny);
    SSxSS_progeny_SS_random = poissrnd(SSxSS_progeny_SS);

    % Find number of RR and RS progeny produced by RRxRS pairing
    % Define mean number of progeny that would be RR
    Probability_RR = 0.5;
    % Create vectors of Probability and progeny with a length of RRxRS
    Prob_arr_RRxRS = Probability_RR(:,ones(1,RRxRS));
    RRxRS_progeny_arr = progeny(:,ones(1,RRxRS));
    % Calculate 'mean' number of progeny produced by pairing by genotype
    RRxRS_progeny_RR = round(binornd(RRxRS_progeny_arr, Prob_arr_RRxRS));
    RRxRS_progeny_RS = RRxRS_progeny_arr - RRxRS_progeny_arr;
    % Use mean number of progeny to calculate stochastic number of progeny
    RRxRS_progeny_RR_random = poissrnd(sum(RRxRS_progeny_RR));
    RRxRS_progeny_RS_random = poissrnd(sum(RRxRS_progeny_RS));
    
    % Find number of SS and RS progeny produced by SSxRS pairing
    % Define mean number of progeny that would be SS
    Probability_SS = 0.5;
    % Create vectors of Probability and progeny with a length of SSxRS
    Prob_arr_SSxRS = Probability_SS(:,ones(1,SSxRS));
    SSxRS_progeny_arr = progeny.*(ones(1,SSxRS));
    % Calculate 'mean' number of progeny produced by pairing by genotype
    SSxRS_progeny_SS = round(binornd(SSxRS_progeny_arr, Prob_arr_SSxRS));
    SSxRS_progeny_RS = SSxRS_progeny_arr - SSxRS_progeny_SS;
    % Use mean number of progeny to calculate stochastic number of progeny
    SSxRS_progeny_SS_random = poissrnd(sum(SSxRS_progeny_SS)); % add up progeny for each genotype, and randomized number
    SSxRS_progeny_RS_random = poissrnd(sum(SSxRS_progeny_RS));
    
    % Find number of SS, RS, and RR progeny produced by RSxRS pairing
    % Calculate total number of progeny produced for pairing without randomness 
    prog = progeny*RSxRS;
    % Calculate total number of progeny for pairing with stochasticity
    RSxRS_progeny = round(poissrnd(prog));
    % Create multinomial probability distribution with the probabilities
    % of each progeny being a specific genotype (RS =.5, SS & RR =.25)
    pr = makedist('Multinomial','Probabilities',[.5,.25,.25]);
    % Create vector with length of total number of progeny for pairing and
    % for each element pick genotype using multinomial probability
    % distribution
    genotype_vector = random(pr, 1, RSxRS_progeny);
    % Count number of progeny by genotype
    RSxRS_progeny_RR_random = sum(genotype_vector==3);
    RSxRS_progeny_SS_random = sum(genotype_vector==2);
    RSxRS_progeny_RS_random = sum(genotype_vector==1);
    
    % Find total amount of progeny for RR, SS, and RS
    tot_progeny_RR = RRxRR_progeny_RR_random + RRxRS_progeny_RR_random + RSxRS_progeny_RR_random;
    tot_progeny_RS = RRxRS_progeny_RS_random + RRxSS_progeny_RS_random + SSxRS_progeny_RS_random + RSxRS_progeny_RS_random;
    tot_progeny_SS = SSxSS_progeny_SS_random + SSxRS_progeny_SS_random + RSxRS_progeny_SS_random;

    %%% STEP 3: Divide progeny between toxic and refuge zone   
    RRref = tot_progeny_RR*Pref;
    SSref = tot_progeny_SS*Pref;
    RSref = tot_progeny_RS*Pref;
    RRtoxic = tot_progeny_RR - RRref;
    SStoxic = tot_progeny_SS - SSref;
    RStoxic = tot_progeny_RS - RSref; 
    
    %%% STEP 4: Selection by zones  
    % Selection - fitness not due to natural enemies 
    RRref_x = RRref* Wxrr_ref;
    RSref_x = RSref*Wxrs_ref;
    SSref_x = SSref*Wxss_ref;
    RRbt_x = RRtoxic*Wxrr_toxic;
    SSbt_x = SStoxic*Wxss_toxic;
    RSbt_x = RStoxic*Wxrs_toxic;

    % Selection due to natural enemies 
    RRref = RRref_x*WErr_ref;
    SSref = SSref_x*WEss_ref;
    RSref = RSref_x*WErs_ref;
    RRtoxic = RRbt_x*WErr_toxic;
    SStoxic = SSbt_x*WEss_toxic;
    RStoxic = RSbt_x*WErs_toxic;

    % Calculate total number of RR, RS, and SS in both field types
    RR = RRtoxic + RRref;
    RS = RStoxic + RSref;
    SS = SStoxic + SSref;

    % Calculate new population 
    N_new = RR + RS + SS;
    
    % Calculate population if change in population only followed
    % logistic growth (include density dependent mortality from carrying
    % capacity and intrinsic growth rate)
    N_logistic = population + intrinsicR * population * (1 - population/K);

    % Recalculate new population including density dependent effects
    % leaving population below carrying capacity at the same proportions of
    % RR, RS, and SS genotypes
    RR = RR*(N_logistic/N_new);
    RS = RS*(N_logistic/N_new);
    SS = SS*(N_logistic/N_new); 
    
    % Calculate total population
    population = RR+RS+SS;
    
    %%% STEP 5: Calculate frequency of resistant and suceptible alleles
    % Find number of each allele after all selection stages
    q_sum = 2*RR + RS;
    p_sum = 2*SS + RS;

    % Finding p and q
    q_freq = (q_sum/(q_sum+p_sum));
    p_freq = (p_sum/(p_sum+q_sum));

    % Displaying generation at which fixation occurs
    if  q_freq >= q_threshold && generations2thresh ==0;
        generations2thresh = i;
        break
    end

    % Saving our population and fequencies info for plots (not used in 
    % all runs)
    q_array = [q_array, q_freq];
end

% If no resistance developed by the end of the run the number of
% generations to threshold is recorded as the number of generations
% simulated
if q_freq <= q_threshold;
    generations2thresh = gen_num;
end


