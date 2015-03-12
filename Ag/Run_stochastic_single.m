%Run_stochastic_single
% q_freq ----> frequency of resistant alleles
% Pref ------> proportion of area that is refuge (no Bt)
% WErr_ref --> fitness of RR in refuge with natural enemies
% WErs_ref --> fitness of RS in refuge with natural enemies
% WEss_ref --> fitness of SS in refuge with natural enemies
% gen_num ---> number of generations (minus 1)

q_s = zeros(1,100);
pops = zeros(1,100);

for nn=1:10                                              %number of trials
    [q_freq, pop] = stochastic_model_05March(0.001,0.3, 0.5, 0.5, 0.5, 199);
    q_s(nn) = q_freq;
    pops(nn) = pop;
end

resistant_not_extinct = sum(q_s>0); % how many trials did the R allele not go to 0
display(resistant_not_extinct)

resistant_fixation = sum(q_s>0.1); % how many trial did the R allele get to 0.1 of pop.
display(resistant_fixation)
