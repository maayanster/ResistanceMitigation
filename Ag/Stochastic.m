%% INPUT
q_freq = 0.09;                        % initial frequency of R allele
p_freq = 1 - q_freq;                  % initial frequency of S allele
mean_no_progeny = int64(10);          % mean progeny produced from one mating
std_dev_progeny = int64(1);           % standard dev of progeny distribution
Wxrr_bt = 0.8;                        % fitness of RR in field w/o natural enemies
Wxss_bt = 0;                          % fitness of SS in field w/o natural enemies
Wxrs_bt = 0;                          % fitness of RS in field w/o natural enemies
WErr_bt = 1;                          % fitness of RR in field due to natural enemies
WErs_bt = 1;                          % fitness of RS in field due to natural enemies
WEss_bt = 1;                          % fitness of SS in field due to natural enemies
WErr_ref = 1;                         % fitness of RR in refuge due to natural enemies
WErs_ref = 1;                         % fitness of RS in refuge due to natural enemies
WEss_ref = 1;                         % fitness of SS in refuge due to natural enemies
Wxrr_ref = 0.79;                      % fitness of RR in refuge w/o natural enemies
Wxrs_ref = 0.79;                      % fitness of RS in refuge w/o natural enemies
Wxss_ref = 0.79;                      % fitness of SS in refuge w/o natural enemies

RR = 0;                               % RR genotype frequency total              
RS = 0;                               % RS genotype frequency total
SS = 0;                               % SS genotype frequency total
RRref = 0;                            % RR genotype frequency refuge             
RSref = 0;                            % RS genotype frequency refuge
SSref = 0;                            % SS genotype frequency refuge
RRbt = 0;                             % RR genotype frequency Bt               
RSbt = 0;                             % RS genotype frequency Bt
SSbt = 0;                             % SS genotype frequency Bt
Pref = 0.2;                           % Proportion of total field that is refuge
Pbt = 1-Pref;                         % Proportion of total field that is Bt

population = int64(10000);           % Initial pest population
K = int64(100000000);                % carrying capacity

i = 0;                                % initialize generation count
%% OUTPUT
p_array=[];
q_array=[];
pop_array=[];

%% CALCULATIONS

while q_freq <= 0.1 && i <= 20;   % while frequency of R allele less than 0.1 or number of generations less than 20 calculate loop
    i = i+1;                      % Count generations
    
    % Determine how many males/females as integers - use binomial 
    % distribution to choose a number of females, assume remainder males)  
    Fem = binornd(population,0.5,1);
    Male = population - Fem;      
               
    % hardy-weinberg ratios are used to calculate expected relative frequency
    % Is it necessary to insert stochasticity into this step? 
    SS = poissrnd(p_freq^2, 1);
    RS = poissrnd(2*q_freq*p_freq, 1);
    RR = poissrnd(q_freq^2, 1);
    
    %%% If error statement to check that RR/RS/SS total 1
    if ((RR+RS+SS)-1) > 0.000001 && ((RR+RS+SS)-1) < -0.000001; 
        error('RR/RS/SS do not equal 1');
    end

    % Sample the number of individuals from each genotype in this generation
    % using poisson distribution
    
    MaleRR = poissrnd(RR*Male, 1);    
    FemRR = poissrnd (RR*Fem, 1);
    MaleRS = poissrnd(RS*Male, 1);
    FemRS = poissrnd (RS*Fem, 1);
    
    %%%% Do a remainder samplemaleSS = males - (malesRS + malesRR)
    MaleSS = poissrnd(SS*Male, 1);
    FemSS = poissrnd (SS*Fem, 1);    %%%% some be the remainders
    
    % using the numbers generated above, calculate how many of each pairing 
    % there would be
    RRxRR = round(FemRR * MaleRR / Males); %% poisson distributionRR prob of each trial by number of trials
    RRxSS = round(FemRR * MaleSS / Males); 
    RRxRS = round(FemRR * MaleRS / Males); %% same as RRxRS
    SSxSS = round(FemSS * MaleSS / Males);
    SSxRS = round(FemSS * MaleRS / Males);
    RSxRS = round(FemRS * MaleRS / Males);
    
    %%% losing stochasticity by rounding replace with poisson distribution
   
  
    % prepare punnet square info on genotypes of progeny for each pairing
    pairings = [RRxRR RRxSS RRxRS SSxSS SSxRS RSxRS];
    Punnet_RR = [1, 0, .5, 0, 0, .25];
    Punnet_RS = [0, 1, .5, 0, .5, .5];
    Punnet_SS = [0, 0, 0, 1, .5, .25];
   
    tot_progeny_RR = 0;
    tot_progeny_RS = 0; 
    tot_progeny_SS = 0;
    
    %for each pairing 1) assign a number of progeny according to normal
    %distribution 2)calculate number of progeny of each genotype
    for jj=1:6
        Progeny_num = normrnd(mean_no_progeny, std_dev_progeny, 1, pairings(jj));   % normrnd(mean, sd, n x m array)
        
        Progeny_RR_vec = round(Punnet_RR(jj)*Progeny_num);  
        Progeny_RR_sum = sum(Progeny_RR_vec);
        Progeny_RS_vec = round(Punnet_RS(jj)*Progeny_num);
        Progeny_RS_sum = sum(Progeny_RS_vec);
        Progeny_SS_vec = round(Punnet_SS(jj)*Progeny_num);
        Progeny_SS_sum = sum(Progeny_SS_vec);
        
        tot_progeny_RR = tot_progeny_RR + Progeny_RR_sum;
        tot_progeny_RS = tot_progeny_RS + Progeny_RS_sum; 
        tot_progeny_SS = tot_progeny_SS + Progeny_SS_sum;
    end    
    
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
    pop_total = RRref + SSref + RSref + RRbt + SSbt + RSbt;
    population = pop_total*(1 - pop_total/K);
    K_survival = population/pop_total;
    
    %find number of each allele after all selection stages
    q_sum = 2*RRref*K_survival + RSref*K_survival + 2*RRbt*K_survival + RSbt*K_survival;
    p_sum = 2*SSref*K_survival + RSref*K_survival + 2*SSbt*K_survival + RSbt*K_survival;
    
    % Finding p and q
    q_freq = q_sum/(q_sum+p_sum);
    p_freq = p_sum/(p_sum+q_sum);
    
    population = (p_sum+q_sum)/2;
    
    %saving our population and fequencies info
    q_array = [q_array, q_freq];
    p_array = [p_array, p_freq];
    pop_array = [pop_array, population];
    
end

gens = 1:length(q_array);
figure
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
