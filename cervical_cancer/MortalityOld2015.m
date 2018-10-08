function [Mu_Tx, Impact, BaseMortality, TxCoverage, StgDist, RMortalityTx, RMortalityNoTx, ModDwellTime] = Mortality(DwellTime, MortalityArray, Prevalence, StageDist, CancerMort,PAR)
global CancerMortality 
global CancerMortalityCalc 
%global StageDistribution
MaxAge = 101;


 
%StageDistribution = StageDist;

%Prashant: variables for proportion of people in corresponding stage 
%surviving after 10 years OR
%probability of people to be in corresponding stage surviving after 10
%years
Healthy = zeros(1,MaxAge);
Insitu = zeros(1,MaxAge);
Local = zeros(1,MaxAge);
Regional = zeros(1,MaxAge);
Distant = zeros(1,MaxAge);
StageIV = zeros(1,MaxAge);

%Prashant: variables for relative survival in corresponding stage. Ratio of
%proability of people in that stage to probability of people in healthy
%stage
RRInsitu = zeros(1,MaxAge);
RRLocal = zeros(1,MaxAge);
RRRegional = zeros(1,MaxAge);
RRDistant = zeros(1,MaxAge);
RRStageIV = zeros(1,MaxAge);

global Mu_TxVal 
Mu_TxVal = zeros(5,1);
Mu_Tx = zeros(5,MaxAge);
global Mu_NoTxVal 
Mu_NoTxVal = zeros(5,1);
Mu_NoTx = zeros(5,MaxAge);
Coverage = zeros(MaxAge,10);

% RInsitu_Tx = .99;
% RLocal_Tx = 0.94;
% RRegional_Tx = 0.60;
% RDistant_Tx = 0.12;
% RStageIV_Tx = 1;
% 
% value = mean( [0.042, 0.093]);
% RTx = [0.006, 0.006, value, 0.275 1];
% value = mean([0.063,.15 1]);
% RNoTx = [0.02, 0.02, value, 0.3 1];
% 
% RTxByNoTx = RTx./RNoTx;


ModDwellTime = DwellTime;


RMortalityTx = zeros(5,MaxAge);
RMortalityNoTx = zeros(5,MaxAge);
 t = 10;
 % Calculating 10 year survival rate%
for i = 1:MaxAge - 10
    
   DwellRate(1,i:i+9)= 1 ./ DwellTime(1,i:i+9);
   DwellRate(2,i:i+9)= 1 ./ DwellTime(2,i:i+9);
   DwellRate(3,i:i+9)= 1 ./ DwellTime(3,i:i+9);
   DwellRate(4,i:i+9)= 1 ./ DwellTime(4,i:i+9);
   if StageDist(1,5) == 0
       DwellRate(5,i:i+9) = 0;
   else
       DwellRate(5,i:i+9)= 1 ./ DwellTime(5,i:i+9);
   end
   
% Ten-year survival without treatment
    if StageDist(1,5) == 0
        stg = 4;
    else
        stg = 5;
    end
   Healthy(1,i) = TenYearSurvival(MortalityArray(1,i:i+9),zeros(1,10));
   Insitu(1,i) = TenYearSurvival(DwellRate(1:stg,i:i+9),MortalityArray(1,i:i+9));
   Local(1,i) = TenYearSurvival(DwellRate(2:stg,i:i+9),MortalityArray(1,i:i+9));
   Regional(1,i) = TenYearSurvival(DwellRate(3:stg,i:i+9),MortalityArray(1,i:i+9));
   Distant(1,i) = TenYearSurvival(DwellRate(4:stg,i:i+9),MortalityArray(1,i:i+9));
   StageIV(1,i) = TenYearSurvival(DwellRate(5:5,i:i+9),MortalityArray(1,i:i+9));
 
% Ten-year RELATIVE survival without treatment   
   Insitu(1,i) = max(Insitu(1,i),Local(1,i));% very high rate for in-situ in older age
   
   RRInsitu(1,i) = Insitu(1,i) / Healthy(1,i);
   RRLocal(1,i) = Local(1,i) / Healthy(1,i);
   RRRegional(1,i) = Regional(1,i) / Healthy(1,i);
   RRDistant(1,i) = Distant(1,i) / Healthy(1,i);
   if StageDist(1,5) == 0
       RRStageIV(1,i) = 0;
   else
       RRStageIV(1,i) = StageIV(1,i) / Healthy(1,i);
   end
 
  %Yearly mortality rate without treatment 
   Mu_NoTx(1,i) = (MortalityArray(1,i)*t - log(RRInsitu(1,i)))/t;
   Mu_NoTx(2,i) = (MortalityArray(1,i)*t - log(RRLocal(1,i)))/t;
   Mu_NoTx(3,i) = (MortalityArray(1,i)*t - log(RRRegional(1,i)))/t;
   Mu_NoTx(4,i) = (MortalityArray(1,i)*t - log(RRDistant(1,i)))/t;
   if StageDist(1,5) == 0
      Mu_NoTx(5,i) = 0;
   else
    Mu_NoTx(5,i) = (MortalityArray(1,i)*t - log(RRStageIV(1,i)))/t;
   end
     
  
end

%Prashant: Assuming that Mortality rate and Relative survival rate is same
%for highest 10 ages.
for i = MaxAge - 10 : MaxAge
   RRInsitu(1,i) =  RRInsitu(1,i - 1); 
   RRLocal(1,i) = RRLocal(1,i - 1);
   RRRegional(1,i) = RRRegional(1,i - 1);
   RRDistant(1,i) = RRDistant(1,i - 1);
   RRStageIV(1,i) = RRStageIV(1,i - 1);
   
   Mu_NoTx(1,i) = Mu_NoTx(1,i-1) ;
   Mu_NoTx(2,i) = Mu_NoTx(2,i-1) ;
   Mu_NoTx(3,i) = Mu_NoTx(3,i-1) ;
   Mu_NoTx(4,i) = Mu_NoTx(4,i-1) ;
   Mu_NoTx(5,i) = Mu_NoTx(5,i-1) ;
end

'Ten year'
%%IF relative survival is available
% InsituCF = mean(RRInsitu);%(1,30:80));
% LocalCF = mean(RRLocal);%(1,30:80));
% RegionalCF = mean(RRRegional);%(1,30:80));
% DistantCF = mean(RRDistant);%(1,30:80));
% StageIVCF = mean(RRStageIV);%(1,30:80));

MuTxAll = [0 0.027 0.062 0.167 0.316];% Goldie et a., 2003
for i = 2:5
MuRatio(i) = (sum(Mu_NoTx(i,:).*Prevalence)/sum(Prevalence)) / MuTxAll(i);
Mu_Tx(i,:) = Mu_NoTx(i,:)./MuRatio(i);
end



% RI_F = RInsitu_Tx / InsituCF;
% RL_F = RLocal_Tx / LocalCF;
% RR_F = RRegional_Tx / RegionalCF;
% RD_F = RDistant_Tx / DistantCF; 
% RIV_F = RStageIV_Tx / StageIVCF; 

% Mu_Tx(2,:) =[0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0.033447184	0.033447184	0.033447184	0.033447184	0.033447184	0.033447184	0.033447184	0.033447184	0.033447184	0.033447184	0.033447184	0.033447184	0.033447184	0.033447184	0.033447184	0.023306763	0.023306763	0.023306763	0.023306763	0.023306763	0.023306763	0.023306763	0.023306763	0.023306763	0.023306763	0.01820388	0.01820388	0.01820388	0.01820388	0.01820388	0.01820388	0.01820388	0.01820388	0.01820388	0.01820388	0.028082431	0.028082431	0.028082431	0.028082431	0.028082431	0.028082431	0.028082431	0.028082431	0.028082431	0.028082431	0.031330762	0.031330762	0.031330762	0.031330762	0.031330762	0.031330762	0.031330762	0.031330762	0.031330762	0.031330762	0.042886322	0.042886322	0.042886322	0.042886322	0.042886322	0.042886322	0.042886322	0.042886322	0.042886322	0.042886322	0.113479195	0.113479195	0.113479195	0.113479195	0.113479195	0.113479195	0.113479195	0.113479195	0.113479195	0.113479195	0.113479195	0.113479195	0.113479195	0.113479195	0.113479195	0.113479195	0.113479195	0.113479195	0.113479195	0.113479195	0.113479195
% ];
% Mu_Tx(3,:) =[0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0.089883399	0.089883399	0.089883399	0.089883399	0.089883399	0.089883399	0.089883399	0.089883399	0.089883399	0.089883399	0.089883399	0.089883399	0.089883399	0.089883399	0.089883399	0.071907235	0.071907235	0.071907235	0.071907235	0.071907235	0.071907235	0.071907235	0.071907235	0.071907235	0.071907235	0.073345056	0.073345056	0.073345056	0.073345056	0.073345056	0.073345056	0.073345056	0.073345056	0.073345056	0.073345056	0.067935474	0.067935474	0.067935474	0.067935474	0.067935474	0.067935474	0.067935474	0.067935474	0.067935474	0.067935474	0.071049478	0.071049478	0.071049478	0.071049478	0.071049478	0.071049478	0.071049478	0.071049478	0.071049478	0.071049478	0.093043023	0.093043023	0.093043023	0.093043023	0.093043023	0.093043023	0.093043023	0.093043023	0.093043023	0.093043023	0.198310643	0.198310643	0.198310643	0.198310643	0.198310643	0.198310643	0.198310643	0.198310643	0.198310643	0.198310643	0.198310643	0.198310643	0.198310643	0.198310643	0.198310643	0.198310643	0.198310643	0.198310643	0.198310643	0.198310643	0.198310643
% ];
% Mu_Tx(4,:) =[0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0.254593135	0.254593135	0.254593135	0.254593135	0.254593135	0.254593135	0.254593135	0.254593135	0.254593135	0.254593135	0.254593135	0.254593135	0.254593135	0.254593135	0.254593135	0.184260655	0.184260655	0.184260655	0.184260655	0.184260655	0.184260655	0.184260655	0.184260655	0.184260655	0.184260655	0.165564417	0.165564417	0.165564417	0.165564417	0.165564417	0.165564417	0.165564417	0.165564417	0.165564417	0.165564417	0.155741014	0.155741014	0.155741014	0.155741014	0.155741014	0.155741014	0.155741014	0.155741014	0.155741014	0.155741014	0.160146478	0.160146478	0.160146478	0.160146478	0.160146478	0.160146478	0.160146478	0.160146478	0.160146478	0.160146478	0.190383582	0.190383582	0.190383582	0.190383582	0.190383582	0.190383582	0.190383582	0.190383582	0.190383582	0.190383582	0.386204307	0.386204307	0.386204307	0.386204307	0.386204307	0.386204307	0.386204307	0.386204307	0.386204307	0.386204307	0.386204307	0.386204307	0.386204307	0.386204307	0.386204307	0.386204307	0.386204307	0.386204307	0.386204307	0.386204307	0.386204307
% ];
% Mu_Tx(5,:) =[0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0.591302312	0.591302312	0.591302312	0.591302312	0.591302312	0.591302312	0.591302312	0.591302312	0.591302312	0.591302312	0.591302312	0.591302312	0.591302312	0.591302312	0.591302312	0.591302312	0.591302312	0.591302312	0.591302312	0.591302312	0.591302312	0.591302312	0.591302312	0.591302312	0.591302312	0.366516293	0.366516293	0.366516293	0.366516293	0.366516293	0.366516293	0.366516293	0.366516293	0.366516293	0.366516293	0.366516293	0.366516293	0.366516293	0.366516293	0.366516293	0.366516293	0.366516293	0.366516293	0.366516293	0.366516293	0.349739996	0.349739996	0.349739996	0.349739996	0.349739996	0.349739996	0.349739996	0.349739996	0.349739996	0.349739996	0.294806655	0.294806655	0.294806655	0.294806655	0.294806655	0.294806655	0.294806655	0.294806655	0.294806655	0.294806655	0.294806655	0.294806655	0.294806655	0.294806655	0.294806655	0.294806655	0.294806655	0.294806655	0.294806655	0.294806655	0.294806655	0.294806655	0.294806655	0.294806655	0.294806655	0.294806655	0.294806655	0.294806655	0.294806655	0.294806655	0.294806655
% ];

% Mu_Tx(2,:) = ones(1,MaxAge) *0.027 ;
% Mu_Tx(3,:) = ones(1,MaxAge) *0.062 ;
% Mu_Tx(4,:) = ones(1,MaxAge) *0.167 ;
% Mu_Tx(5,:) = ones(1,MaxAge) *0.316 ;

Mu_Tx(1,:) = Mu_NoTx(1,:).*Mu_Tx(2,:)./Mu_NoTx(2,:);
for i = 1:MaxAge
        if StageDist(1,5) == 0
           Mu_Tx(5,i) = 0;
        end
  
    
     
%         R1 = RTx(2) / RTx(3) ;
%         R2 = RTx(2) / RTx(4) ;
%         R3 = RTx(2) / RTx(1) ;
%         R4 = RTx(2) / RTx(5) ;
%     if ((CancerMort(i,1) + MortalityArray(1,i)) / (StageDist(1,2)+ (StageDist(1,3)/R1)+ (StageDist(1,4)/R2 ) + (StageDist(1,1)/R3 ) + (StageDist(1,5)/R4 ))) > Mu_Tx(2,i)
%         Mu_Tx(2,i) = ((CancerMort(i,1) + MortalityArray(1,i)) / (StageDist(1,2)+ (StageDist(1,3)/R1)+ (StageDist(1,4)/R2 ) + (StageDist(1,1)/R3 ) + (StageDist(1,5)/R4 )));
%         Mu_Tx(3,i) = Mu_Tx(2,i)/ R1;
%         Mu_Tx(4,i) = Mu_Tx(2,i)/ R2;
%         Mu_Tx(1,i) = Mu_Tx(2,i)/ R3;
%         if StageDist(1,5) == 0
%              Mu_Tx(5,i) = 0;
%         else
%              Mu_Tx(5,i) = Mu_Tx(2,i)/ R4;
%         end
%          
%         ValueF(1:5,1) = (Mu_Tx(:,i) ./ RTxByNoTx(:))./ Mu_NoTx(:,i);
%         
%         Mu_NoTx(:,i) = Mu_Tx(:,i) ./ RTxByNoTx(:);
%         
%         ModDwellTime(2:5,i) =  ModDwellTime(2:5,i)./ValueF(2:5,1);
%        
%    
%     end
%         

     % 10-year Relative surival with_Tx
   
        RI(1,i) =  min(1,exp(MortalityArray(1,i)*t- Mu_Tx(1,i)*t));
        RL(1,i) =  min(1,exp(MortalityArray(1,i)*t- Mu_Tx(2,i)*t));
        RR(1,i) =  min(1,exp(MortalityArray(1,i)*t- Mu_Tx(3,i)*t));
        RD(1,i) =  min(1,exp(MortalityArray(1,i)*t- Mu_Tx(4,i)*t));
        RIV(1,i) =  min(1,exp(MortalityArray(1,i)*t- Mu_Tx(5,i)*t));
  
    
    CancerMortality = CancerMort(i,1) + MortalityArray(1,i);
    Mu_TxVal(:,1) = Mu_Tx(:,i);
    Mu_NoTxVal(:,1) = Mu_NoTx(:,i);
    
    if StageDist(1,1) == 0
        uVal = 0;
    else
        uVal = 1;
    end
    
    if StageDist(1,5) == 0
        uVal5 = 0;
    else
        uVal5 = 1;
    end
    % Assuming Tx coverage is 100% 
    LowerBound = [uVal 1 1 1 uVal5 0 0 0 0 0 ];
    UpperBound = [uVal 1 1 1 uVal5 uVal 1 1 1 uVal5 ];
    options = optimoptions(@fmincon,'MaxIter', 1000, 'TolX', 0.001);
    %Covergae: 1:4= coverage of Tx in  stage insitu, local, regional, and distant; 4:6 =
    %Stage distribution
    [Coverage(i, 1:10), error(i)] = fmincon(@FuncCancerMort, [0 1 1 1 0 0 0.3 0.4 .3 0], [], [], [0 0 0 0 0 1 1 1 1 1], [1],LowerBound,UpperBound,[] );
    %Coverage(6,i) = 1- Coverage(5,i) - Coverage (4,i);
    Cancer(1,i) = CancerMortalityCalc;

    for j = 1: 10
    if Coverage(i, j) < 0.01
       Coverage(i, j) = 0; 
    end
    end
end
    %to view wrtiting ten year survival probability
    dlmwrite('TenYearSurvival.csv',RI);
    dlmwrite('TenYearSurvival.csv',RL, '-append');
    dlmwrite('TenYearSurvival.csv',RR, '-append');
    dlmwrite('TenYearSurvival.csv',RD, '-append');
    dlmwrite('TenYearSurvival.csv',RIV, '-append');
    dlmwrite('TenYearSurvival.csv',Prevalence, '-append');        
    
    %not sure if this formula is correct.
     dlmwrite('TenYearSurvival.csv',dot(RI,Prevalence)/sum(Prevalence), '-append');
     dlmwrite('TenYearSurvival.csv',dot(RL,Prevalence)/sum(Prevalence), '-append');
     dlmwrite('TenYearSurvival.csv',dot(RR,Prevalence)/sum(Prevalence), '-append');
     dlmwrite('TenYearSurvival.csv',dot(RD,Prevalence)/sum(Prevalence), '-append');
     dlmwrite('TenYearSurvival.csv',dot(RIV,Prevalence)/sum(Prevalence), '-append');

      %Subtracting disease free mortality from cancer mortality as it is added
      %back in Spectrum


    Mu_Tx(1,:) = max(0,Mu_Tx(1,:) - MortalityArray(1,:));
    Mu_Tx(2,:) = max(0,Mu_Tx(2,:) - MortalityArray(1,:));
    Mu_Tx(3,:) = max(0,Mu_Tx(3,:) - MortalityArray(1,:));
    Mu_Tx(4,:) = max(0,Mu_Tx(4,:) - MortalityArray(1,:));
    Mu_Tx(5,:) = max(0,Mu_Tx(5,:) - MortalityArray(1,:));

    Mu_NoTx(1,:) = max(0,Mu_NoTx(1,:) - MortalityArray(1,:));
    Mu_NoTx(2,:) = max(0,Mu_NoTx(2,:) - MortalityArray(1,:));
    Mu_NoTx(3,:) = max(0,Mu_NoTx(3,:) - MortalityArray(1,:));
    Mu_NoTx(4,:) = max(0,Mu_NoTx(4,:) - MortalityArray(1,:));
    Mu_NoTx(5,:) = max(0,Mu_NoTx(5,:) - MortalityArray(1,:));

    Mu_Tx(1,80:end) = Mu_Tx(1,80);
    Mu_Tx(2,80:end) = Mu_Tx(2,80);
    Mu_Tx(3,80:end) = Mu_Tx(3,80);
    Mu_Tx(4,80:end) = Mu_Tx(4,80);
    Mu_Tx(5,80:end) = Mu_Tx(5,80);
    
    Mu_NoTx(1,80:end) = Mu_NoTx(1,80);
    Mu_NoTx(2,80:end) = Mu_NoTx(2,80);
    Mu_NoTx(3,80:end) = Mu_NoTx(3,80);
    Mu_NoTx(4,80:end) = Mu_NoTx(4,80);
    Mu_NoTx(5,80:end) = Mu_NoTx(5,80);
   
    RMortalityTx(1,:)= Mu_Tx(1,:)./ MortalityArray(1,:);
    RMortalityTx(2,:)= Mu_Tx(2,:)./ MortalityArray(1,:);
    RMortalityTx(3,:)= Mu_Tx(3,:)./ MortalityArray(1,:);
    RMortalityTx(4,:)= Mu_Tx(4,:)./ MortalityArray(1,:);
    RMortalityTx(5,:)= Mu_Tx(5,:)./ MortalityArray(1,:);
    
    RMortalityNoTx(1,:)= Mu_NoTx(1,:)./ MortalityArray(1,:);
    RMortalityNoTx(2,:)= Mu_NoTx(2,:)./ MortalityArray(1,:);
    RMortalityNoTx(3,:)= Mu_NoTx(3,:)./ MortalityArray(1,:);
    RMortalityNoTx(4,:)= Mu_NoTx(4,:)./ MortalityArray(1,:); 
    RMortalityNoTx(5,:)= Mu_NoTx(5,:)./ MortalityArray(1,:); 

    

% for i = 1: MaxAge
%     for j = 1: 4
%         if ((MortalityArray(1,i)*t - log(RRInsitu(1,i)))/t) < Mu_NoTx(1,i) 
%             
%                LowerBound = [0 ];
%                UpperBound = [40];
%                 options = optimoptions(@fmincon,'MaxIter', 1000, 'TolX', 0.001);
%                [ModDwellTime(j,i)] = fmincon(@FuncCancerMort, [], [], [], [0 0 0 0 1 1 1 1], [1],LowerBound,UpperBound,[] );
% 
% end

% Coverage(80:100,1) = Coverage(80,1);
% Coverage(80:100,2) = Coverage(80,2);
% Coverage(80:100,3) = Coverage(80,3);
% Coverage(80:100,4) = Coverage(80,4);
% Calculating basecase mortality, i.e., the mortality weighted by
% treatment coverage
% 

BaseMortality = zeros(size(Mu_Tx,1), size(Mu_Tx,2));
BaseMortality(1,:) = (Coverage(:,1)' .*  Mu_Tx(1,:) + (1- Coverage(:,1)').* Mu_NoTx(1,:))  ;
BaseMortality(2,:) = (Coverage(:,2)' .*  Mu_Tx(2,:) + (1- Coverage(:,2)').* Mu_NoTx(2,:))  ;
BaseMortality(3,:) = (Coverage(:,3)' .*  Mu_Tx(3,:) + (1- Coverage(:,3)').* Mu_NoTx(3,:))  ;
BaseMortality(4,:) = (Coverage(:,4)' .*  Mu_Tx(4,:) + (1- Coverage(:,4)').* Mu_NoTx(4,:)) ;
BaseMortality(5,:) = (Coverage(:,5)' .*  Mu_Tx(5,:) + (1- Coverage(:,5)').* Mu_NoTx(5,:)) ;

% estimating increase in weighted mortality for unit change in coverage
Impact (1,:) = (Mu_Tx(1,:) - Mu_NoTx(1,:)) ./ BaseMortality(1,:);
Impact (2,:) = (Mu_Tx(2,:) - Mu_NoTx(2,:)) ./ BaseMortality(2,:);
Impact (3,:) = (Mu_Tx(3,:) - Mu_NoTx(3,:)) ./ BaseMortality(3,:);
Impact (4,:) = (Mu_Tx(4,:) - Mu_NoTx(4,:)) ./ BaseMortality(4,:);
Impact (5,:) = (Mu_Tx(5,:) - Mu_NoTx(5,:)) ./ BaseMortality(5,:);

TxCoverage(:,:) = max(0,Coverage(:,1:5));
StgDist(:,:) =  max(0,Coverage(:,6:10));
close all
figure
hold on
plot(1:MaxAge, Healthy(1,1:MaxAge), 'r');
plot(1:MaxAge, Insitu(1,1:MaxAge), 'y');
plot(1:MaxAge, Local(1,1:MaxAge), 'b');
plot(1:MaxAge, Regional(1,1:MaxAge), 'g');
plot(1:MaxAge, Distant(1,1:MaxAge), 'c');
plot(1:MaxAge, StageIV(1,1:MaxAge), 'm');
legend('healthy','in-situ', 'local', 'regional', 'distant', 'stageIV');

% plot(30:80, T_Local(1,30:80), 'b-');
% plot(30:80, T_Regional(1,30:80), 'g-');
% plot(30:80, T_Distant(1,30:80), 'c-');
title('5-year Survival')


figure
hold on
plot(1:MaxAge, RMortalityNoTx(1,1:MaxAge), 'y');
plot(1:MaxAge, RMortalityNoTx(2,1:MaxAge), 'b');
plot(1:MaxAge, RMortalityNoTx(3,1:MaxAge), 'g');
plot(1:MaxAge,RMortalityNoTx(4,1:MaxAge), 'c');
plot(1:MaxAge,RMortalityNoTx(5,1:MaxAge), 'm');

legend('in-situ', 'local', 'regional', 'distant', 'stageIV');
plot(1:MaxAge, RMortalityTx(1,1:MaxAge), 'y+');
plot(1:MaxAge, RMortalityTx(2,1:MaxAge), 'b+');
plot(1:MaxAge, RMortalityTx(3,1:MaxAge), 'g+');
plot(1:MaxAge, RMortalityTx(4,1:MaxAge), 'c+');
plot(1:MaxAge, RMortalityTx(5,1:MaxAge), 'm+');

title('yearly relative mortality in cancer compared to non-cancer persons ("-"without Tx "+" with Tx)')

figure
hold on
plot(1:MaxAge, RRInsitu(1,1:MaxAge), 'y');
plot(1:MaxAge, RRLocal(1,1:MaxAge), 'b');
plot(1:MaxAge, RRRegional(1,1:MaxAge), 'g');
plot(1:MaxAge,RRDistant(1,1:MaxAge), 'c');
plot(1:MaxAge,RRStageIV(1,1:MaxAge), 'm');

legend('in-situ', 'local', 'regional', 'distant', 'stageIV');
plot(1:MaxAge, RI(1,1:MaxAge), 'y+');
plot(1:MaxAge, RL(1,1:MaxAge), 'b+');
plot(1:MaxAge, RR(1,1:MaxAge), 'g+');
plot(1:MaxAge, RD(1,1:MaxAge), 'c+');
plot(1:MaxAge, RD(1,1:MaxAge), 'm+');
title('10-year relative survival in cancer compared to non-cancer persons ("-"without Tx "+" with Tx)')

figure
hold on
plot(1:MaxAge, Coverage(1:MaxAge,1), 'b');
plot(1:MaxAge, Coverage(1:MaxAge,2), 'g');
plot(1:MaxAge,Coverage(1:MaxAge,3), 'c');
plot(1:MaxAge,Coverage(1:MaxAge,4), 'r');
plot(1:MaxAge,Coverage(1:MaxAge,5), 'm');
legend('in-situ', 'local', 'regional', 'distant', 'stageIV');

title(' Treatment Coverage')

figure
hold on
plot(1:MaxAge, Coverage(1:MaxAge,6), 'b');
plot(1:MaxAge, Coverage(1:MaxAge,7), 'g');
plot(1:MaxAge,Coverage(1:MaxAge,8), 'c');
plot(1:MaxAge,Coverage(1:MaxAge,9), 'r');
plot(1:MaxAge,Coverage(1:MaxAge,10), 'm');
legend('in-situ', 'local', 'regional', 'distant', 'stageIV');

title(' Stage Dist')

figure
hold on
plot(1:MaxAge, BaseMortality(1,1:MaxAge), 'b');
plot(1:MaxAge, BaseMortality(2,1:MaxAge), 'g');
plot(1:MaxAge,BaseMortality(3,1:MaxAge), 'c');
plot(1:MaxAge,BaseMortality(4,1:MaxAge), 'y');
plot(1:MaxAge,BaseMortality(5,1:MaxAge), 'm');
legend('in-situ', 'local', 'regional', 'distant', 'stageIV');

% plot(1:MaxAge,BaseMortality(4,1:MaxAge), 'r');
% plot(1:MaxAge, CancerMort(1:MaxAge,1), 'c+');

title(' Cancert mortality by stage')
figure
hold on
plot(1:MaxAge, Mu_Tx(1,1:MaxAge), 'b+');
plot(1:MaxAge, Mu_Tx(2,1:MaxAge), 'g+');
plot(1:MaxAge, Mu_Tx(3,1:MaxAge), 'y+');
plot(1:MaxAge, Mu_Tx(4,1:MaxAge), 'r+');
plot(1:MaxAge, Mu_Tx(5,1:MaxAge), 'm+');
legend('in-situ', 'local', 'regional', 'distant', 'stageIV');

plot(1:MaxAge, Mu_NoTx(1,1:MaxAge), 'b');
plot(1:MaxAge, Mu_NoTx(2,1:MaxAge), 'g');
plot(1:MaxAge, Mu_NoTx(3,1:MaxAge), 'y');
plot(1:MaxAge, Mu_NoTx(4,1:MaxAge), 'r');
plot(1:MaxAge, Mu_NoTx(5,1:MaxAge), 'm');

title(' Mu-Treatment(+) and No Treatment(-)')

figure
hold on
plot(1:80, CancerMort(1:80,1) - MortalityArray(1,1:80)' , 'c');
plot(1:80, CancerMort(1:80,1), 'c+');
% title(' Mu-Cancer')
% 
% 
% 
% 
% figure
% hold on
%plot(1:MaxAge, MortalityArray(1,1:MaxAge), 'r');
title(' Mortalities')

[~,~,TEMPMORT1]=xlsread('Mort Results.xlsx','10-year Survial');
[~,~,TEMPMORT2]=xlsread('Mort Results.xlsx','yrly reltv mort comp non cancer');
[~,~,TEMPMORT3]=xlsread('Mort Results.xlsx','10-yr reltv mort comp non cance');
[~,~,TEMPMORT4]=xlsread('Mort Results.xlsx','Treatment Coverage');
[~,~,TEMPMORT5]=xlsread('Mort Results.xlsx','Stage Dist');
[~,~,TEMPMORT6]=xlsread('Mort Results.xlsx','Cancer Mort by stage');
[~,~,TEMPMORT7]=xlsread('Mort Results.xlsx','Mu Treat and No Treat');
[~,~,TEMPMORT8]=xlsread('Mort Results.xlsx','Mortalities');

CountryIndex = PAR{13};
CountryName = PAR{14};

 [A,B,C,D,E,F,G,H] = MikeWriteTemplateMort(TEMPMORT1,TEMPMORT2,TEMPMORT3,TEMPMORT4,TEMPMORT5,TEMPMORT6,TEMPMORT7,TEMPMORT8,Healthy,Insitu,Local,Regional,Distant,StageIV,RMortalityNoTx,RMortalityTx,RRInsitu,RRLocal,RRRegional,RRDistant,RRStageIV,RI,RL,RR,RD,RIV,Coverage,BaseMortality,Mu_Tx,Mu_NoTx,CancerMort,MortalityArray,CountryIndex,CountryName)


 xlswrite('Mort Results.xlsx',A,'10-year Survial');
 xlswrite('Mort Results.xlsx',B,'yrly reltv mort comp non cancer');
 xlswrite('Mort Results.xlsx',C,'10-yr reltv mort comp non cance');
 xlswrite('Mort Results.xlsx',D,'Treatment Coverage');
 xlswrite('Mort Results.xlsx',E,'Stage Dist');
 xlswrite('Mort Results.xlsx',F,'Cancer Mort by stage');
 xlswrite('Mort Results.xlsx',G,'Mu Treat and No Treat');
 xlswrite('Mort Results.xlsx',H,'Mortalities');
 

  'Mort Results Written'


%pause
end

function Error = FuncCancerMort(coverage)
global Mu_TxVal 
global Mu_NoTxVal
%global StageDistribution 
global CancerMortality
global CancerMortalityCalc
CancerMortalityCalc =  (coverage(1,1) *  Mu_TxVal(1,1) + (1- coverage(1,1))* Mu_NoTxVal(1,1)) * coverage(1,6) + ...
                       (coverage(1,2) *  Mu_TxVal(2,1) + (1- coverage(1,2))* Mu_NoTxVal(2,1)) * coverage(1,7) + ...
                       (coverage(1,3) *  Mu_TxVal(3,1) + (1- coverage(1,3))* Mu_NoTxVal(3,1)) * coverage(1,8) + ...
                       (coverage(1,4) *  Mu_TxVal(4,1) + (1- coverage(1,4))* Mu_NoTxVal(4,1)) * coverage(1,9) + ...
                       (coverage(1,5) *  Mu_TxVal(5,1) + (1- coverage(1,5))* Mu_NoTxVal(5,1)) * coverage(1,10);
       
Error = CancerMortality - CancerMortalityCalc;
Error = Error * Error;
end


function Proportion = TenYearSurvival(RatesArray, Mu)
dt =24;
RatesArray = RatesArray /dt;
Mu = Mu / dt;
N = zeros (1,size(RatesArray,1));
N(1,1) = 100000; %simulated population size
for t = 1:10 * dt
    %Prashant: no of people leaving from this stage due to moving to next stage + dying
    N(1,1) = N(1,1) - N(1,1) * (RatesArray(1,max(1,round(t/dt)))+ Mu(1,max(1,round(t/dt))));
    for i = size(RatesArray):-1:2
        %Prashant: no of people already in this stage + no of people coming from prev
        %Prashant: stage + no of people dying from this stage
    N(1,i)= N(1,i) + N(1,i-1) * RatesArray(i-1,max(1,round(t/dt))) - N(1,i) * (RatesArray(i, max(1,round(t/dt))) + Mu(1, max(1,round(t/dt))));
    end
end

%Prashant: returns proportion of people alive in any stage to the initial population
Proportion = sum(sum(N))/100000; 
end



