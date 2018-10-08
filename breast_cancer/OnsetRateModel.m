clear;
clc;

%% NOTE : CODE till 13/06/17
% code is not consistent in names with the previous code, have to make
% thode changes.

%% Desription:
% following code is based on the two-step markov model developed by Prof.
% Chaitra G et. al. for parameterization of cancer state transition for
% LMIC (low and middle income countries).

% The whole model is split into two sub-models, first for calculation of
% onset rate and second for calculating the rate of diagnosis

%% MODEL 1:
% model 1 is a stochastic process {Yt: t>=0}, defined over state space as
% follows,
% s = { Ha, Ua, Da, M }, H, U, D and M stands for Healthy state,
% Undiagnosed state, Diagnosed state and Mortal state respectively; 'a'
% stands for age. TPM cam be referred from the paper

% Aim of model 1: Calculation of onset rate for each value of 'a' (age)

% Overview of algorithm of model 1:
% STEP 1: initialization step
% STEP 2: calculation of P_haua and onset rate from that value
% STEP 3: calculation of P_steady state_ua-1
% STEP 4: calculation of P_steady state_ha
% STEP 5: increment of a by 1, if a<a_max go to step 2, else stop

%% inputs needed for the algorithm
% 1. incidence estimates, inci
% 2. fraction of population in age 'a', c
% 3. fraction of population not living with disease diagnosed, i.e either
% healthy or i undiagnosed state, A
% 4. age matrix
% 5. p_steady_u at age zero
% 6. p_steady_h at age zero
% 7. A at age zero
% 8. transition rate from H to D
% 9. mortality rate from D to M
% 10. Probability of transition from H to U at age zero

% maximum age under consideration
a_max=100;

% defining incidence matrix
inci=zeros(1,a_max+1);

% defining fraction of population in age 'a'
C=zeros(1,a_max+1);

% defining fraction of population not living with disease diagnosed, A
A=zeros(1,a_max+1);
A(1,1)=1; % i.e. at age zero, probability of person being healthy is one, also A0=1

% assume that p_steady_u0=0
p_steady_u=zeros(1,a_max+1);

% assume that p_steady_h0=A0
p_steady_h=zeros(1,a_max+1);
p_steady_h(1,1)=1; % i.e. at age zero, probability of person being healthy is one, also A0=1

% defining transition probability matrix
p_trans_hu=zeros(1,a_max+1);

% lambda (transition rate from H o D)
l=zeros(1,a_max);

% mu (mortality rate from D to M)
mu=zeros(1,a_max+1);

% initializing algorithm
onset=zeros(1,a_max+1);
a=1;

while a<a_max
    
    % Remark 3 from the paper
    par1=inci(1,a+1)*C(1,a+1); % 'a+1' because every array is starting from a=0 (indexing adjustments)
    
    % defining parameter 2 for P(s>=a-k)
    par2_1=0;   % this parameter does summation of mu till a+1
    par2_3=0;   % this parameter does summation of mu till a (will be used in calculation of p_steady_h)
    for i=(0:a+1)
        par2_1=par2_1+mu(1,i+1);
        if i==a
            par2_3=par2_1;
        else
        end
    end
    par2_2=1-expcdf(1,1/par2_1);
    
    % from Lemma 2 of paper
    par3=0;
    for i=(0:a-1)
        par3=par3+(p_steady_h(1,i+1)*p_trans_hu(1,i+1)*(expcdf(a-i+1,1/l(1,i+1))-expcdf(a-i,1/l(1,i+1)))*par2_2);
    end
    
    par4=A(1,a+1)*expcdf(1,1/l(1,a+1))*(1-expcdf(1,1/mu(1,a+1)));
    
    % calculation of transition probability
    p_trans_hu(1,a+1)=par1*par3/(par4-par1);
    
    % calculation of onset rate
    onset(1,a+1)=-1*log(1-p_trans_hu(1,a+1));
    
    % caculation of steady state probability of being in state 'U' at age
    % 'a-1'
    p_steady_u(1,a)=A(1,a)-p_steady_h(1,a);
    
    % caculation of steady state probability of being in state 'H' at age
    % 'a'
    par1=0;
    for i=(0:a-1)
        % using par2_3 here to calculate numerator, par2_3 is summation of
        % mu from 0 to a
        par1=par1+(p_steady_h(1,i+1)*p_trans_hu(1,i+1)*(1-expcdf(a-i,1/l(1,i+1)))*(1-expcdf(1,1/par2_3)));
    end
    p_steady_h(1,a+1)=par1/(1+p_trans_hu(1,a+1));
    
    % increasing a by 1
    a=a+1;
end




