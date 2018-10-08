function [RES]=MarkovChainCountry2(PAR)

   
    RES=[]; %empty array to store results

    close all
    nsims=10000; %number of simulations

    %Country-specific parameters
    probAgeA=PAR{1}; %distribution of ages in the population (14 years age range groups)
    MortalityArray=PAR{2}; %mortality by age
    
    for i = 80:101
        MortalityArray(i)=MortalityArray(i-1)* MortalityArray(79)/MortalityArray(78);
    end
    
    Births=PAR{3};
    prevalenceActualF=PAR{4}; %prevalence of the country by age (14 years range)
    %prevN=prevalenceActual/sum(prevalenceActual);

    %Prevalence and case fatality rate: Age Group is different from Spectrum age group so converting
    AgeArrayL = PAR{7};
    AgeArrayU = PAR{8};
    AgeArrayM=(AgeArrayU+AgeArrayL)/2;


    AgeArrayLM =[1    5   15  30	45	60	70	80];
    AgeArrayUM =[4    14  29 	44	59	69	79	100];
    AgeArrayMM=(AgeArrayUM+AgeArrayLM)/2;

    CFR = interp1(AgeArrayM,PAR{5},[0:100],'PCHIP'); %interploting grouped age to 0:100
    %CFR(1,1:14)=0;
    cancerMortalityArray =-log(1-min(.99,CFR)); %mortality by age (14 years range) -> converting from exp 
    %cancerMortalityArray(1,80:end) = cancerMortalityArray(1,79);
    %cancerMortalityArray = max(cancerMortalityArray, MortalityArray);
    PAR{5} = cancerMortalityArray; %assigning modified data to par{5}
    StageDist=PAR{6}; %distribution of stages (insitu, local, regional, distant)
    StageDist(1,5)=0; 

    % CountryIndex=PAR{9};

    %Interpolation to single age years, flatline after 80 years
    PrevalenceF = interp1(AgeArrayM,prevalenceActualF,[0:100],'PCHIP');
    %PrevalenceF(80:end)=PrevalenceF(80);
    PrevalenceF= max(0,PrevalenceF);


    %Spectrum Age groups
    AgeArrayLower=[1    5	10	15  20	25	30	40	50	60	70	80];
    AgeArrayUpper=[4	9	14  19 	24	29	39	49	59	69	79	100];
    AgeArrayMid=((AgeArrayUpper+AgeArrayLower)/2);

    prevalenceActual = zeros(size(prevalenceActualF,1),size(AgeArrayLower,2));

    for i = 1 : size(prevalenceActual,2)
        indA=AgeArrayLower(i):AgeArrayUpper(i);
        prevalenceActual(1,i)= dot(PrevalenceF(indA),probAgeA(indA))/sum(probAgeA(indA)) ;%PrevalenceF(AgeArrayMid(i));%
    end
    %prevalenceActual(1,1:3) =0;
    PAR{4}=prevalenceActual;

    %Progression times in undiagnoses cases
    DwellTimes = zeros(size(AgeArrayLower,1),size(AgeArrayLower,2));
    for i = 1 : size(AgeArrayLower,2) 
      
%       %below commented by Prashant 12 July 2016 and using newer values
%       DwellTimes(1,i)=3.4;
%       DwellTimes(2,i)=5;
%       DwellTimes(3,i)=0.95;
%       DwellTimes(4,i)=0.0667;
        %Using below Dwell times from the cancer assumptions file ~\Dropbox\Cancer manuscripts\Breast and CRC_ cancer impact model assumptions - 7-11-2016.xlsx
         DwellTimes(1,i)=3.4;
         DwellTimes(2,i)=5.1600;
         DwellTimes(3,i)=1.5873;
         DwellTimes(4,i)=1.1111;
         
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Interpolation to single age years, flatline after 80 years

    DwellTimes1 = interp1(AgeArrayMid,DwellTimes(1,:),[0:100],'PCHIP');
    % DwellTimes1(80:end)=DwellTimes1(80);

    DwellTimes2 = interp1(AgeArrayMid,DwellTimes(2,:),[0:100],'PCHIP');
    % DwellTimes2(80:end)=DwellTimes2(80);

    DwellTimes3 = interp1(AgeArrayMid,DwellTimes(3,:),[0:100],'PCHIP');
    % DwellTimes3(80:end)=DwellTimes3(80);

    DwellTimes4 = interp1(AgeArrayMid,DwellTimes(4,:),[0:100],'PCHIP');
    % DwellTimes4(80:end)=DwellTimes4(80);

    DwellTimesF=[DwellTimes1;DwellTimes2;DwellTimes3;DwellTimes4];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Calculates 10-year relative survival with and without treatment, which are
    %used in Excel file to calculate  relative mortality avergaed over the 3
    %countries. 
    %Also calculates coverage and stage distributions, which would be good for
    %verfication purposes.
    [Mu_Tx, ImpactAge, BaseMortality, TxCoverage, StgDist, RelativeMortTx, RelativeMortNoTx] = Mortality(DwellTimesF, MortalityArray, PrevalenceF, StageDist, smooth(PAR{5},'lowess'),PAR);

  
    %average values into age groups
    MortA=[];
    MortC=[];
    for J=1:length(AgeArrayLower)
        indA=AgeArrayLower(J):AgeArrayUpper(J);
        MortA(J)=dot(MortalityArray(indA),probAgeA(indA))/sum(probAgeA(indA));
        MortC(J)=dot(cancerMortalityArray(indA),probAgeA(indA))/sum(probAgeA(indA));
        for i = 1:4 % if the denominator is equal to 0, all the values are 0.
            if sum(PrevalenceF(indA))==0
                BaseMortalityC(i,J)=0;
                Impact(i,J)=0;
                TxCoverage_A(i,J)=0;
                StgDistribution(i,J)=0;
                DwellTimes(i,J)=0;
                %%%%%%%
                MORT(i,J)=0;
                BASE(i,J)=0;
            else
            BaseMortalityC(i,J) = dot(RelativeMortNoTx(i,indA),PrevalenceF(indA))/sum(PrevalenceF(indA));
            Impact(i,J) =         dot(RelativeMortTx(i,indA),PrevalenceF(indA))/sum(PrevalenceF(indA));
            TxCoverage_A(i,J) = dot(TxCoverage(indA,i)',PrevalenceF(indA))/sum(PrevalenceF(indA));
            StgDistribution(i,J) = dot(StgDist(indA,i)',PrevalenceF(indA))/sum(PrevalenceF(indA));
            DwellTimes(i,J) = dot(DwellTimesF(i, indA),PrevalenceF(indA))/sum(PrevalenceF(indA)); 
            %%%%%%%
            MORT(i,J)= dot((MortalityArray(indA)+MortalityArray(indA).*RelativeMortTx(i,indA)),PrevalenceF(indA))/sum(PrevalenceF(indA));
            BASE(i,J) = dot((MortalityArray(indA)+MortalityArray(indA).*RelativeMortNoTx(i,indA)),PrevalenceF(indA))/sum(PrevalenceF(indA));
            end
        end
      end
    MortA(end)=MortA(end-1);                                                                                                              
    MortC(end)=MortC(end-1);


    %TimeToDiagnosis
    TimeToDiagnosis = zeros(size(DwellTimes,1),size(DwellTimes,2));
    for i = 1 : size(DwellTimes,2)

        sumTime = zeros(size(DwellTimes,1));
        for j = 1:nsims
           randNum = rand;
           Insitu = - log(randNum) * DwellTimes(1,i);
           Local = - log(randNum) * DwellTimes(2,i);
           Regional = - log(randNum) * DwellTimes(3,i);
           Distant = - log(randNum) * DwellTimes(4,i); 

           sumTime(1) = sumTime(1) + 0.5*Insitu;
           sumTime(2) = sumTime(2) + Insitu + 0.5*Local;
           sumTime(3) = sumTime(3) + Insitu + Local + 0.5*Regional;
           sumTime(4) = sumTime(4) + Insitu + Local + Regional + 0.5*Distant; 
        end
       TimeToDiagnosis(1,i) = sumTime(1) / nsims;
       TimeToDiagnosis(2,i) = sumTime(2) / nsims;
       TimeToDiagnosis(3,i) = sumTime(3) / nsims;
       TimeToDiagnosis(4,i) = sumTime(4) / nsims;
    end

    TimeToDiagnosis

   %Modified TimeToDiagnosis (using DCO to estimate event of first occurence )
ModTimeToDiagnosis = zeros(size(DwellTimes,1),size(DwellTimes,2));

CumStageDISTwithDCO = StageDist;
DCO = 0.1;% proportion DCO  
%http://ci5.iarc.fr/CI5-X/PDF/INDICES/I_10.pdf
%http://ci5.iarc.fr/CI5-X/Pages/Indices_sel.aspx
for i = 1: size (StageDist,2)
    CumStageDISTwithDCO(i) = StageDist(i) / (sum(StageDist(i:end))+ DCO);
end 
%CumStageDISTwithDCO = CumStageDISTwithDCO / (1+DCO);

% for i = 1 : size(DwellTimes,2)
%    ModTimeToDiagnosis(1,i) = DwellTimes(1,i) * (1 - CumStageDISTwithDCO (1))/CumStageDISTwithDCO(1);
%    ModTimeToDiagnosis(2,i) = DwellTimes(1,i) + DwellTimes(2,i) * (1 - CumStageDISTwithDCO (2))/CumStageDISTwithDCO(2);
%    ModTimeToDiagnosis(3,i) = DwellTimes(1,i) + DwellTimes(2,i) + DwellTimes(3,i) * (1 - CumStageDISTwithDCO (3))/CumStageDISTwithDCO(3);
%    ModTimeToDiagnosis(4,i) = DwellTimes(1,i) + DwellTimes(2,i) + DwellTimes(3,i) + DwellTimes(4,i) * (1 - CumStageDISTwithDCO (4))/CumStageDISTwithDCO(4);
%    ModTimeToDiagnosis(5,i) = DwellTimes(1,i) + DwellTimes(2,i) + DwellTimes(3,i) + DwellTimes(4,i)+ DwellTimes(5,i) * (1 - CumStageDISTwithDCO (5))/CumStageDISTwithDCO(5);
% end 

% for i = 1 : size(DwellTimes,2)
%    ModTimeToDiagnosis(1,i) = (1 - CumStageDISTwithDCO (1))/((1/DwellTimes(1,i) + MortA(1,i)) * CumStageDISTwithDCO(1));
%    ModTimeToDiagnosis(2,i) = (1 - CumStageDISTwithDCO (2))/((1/(DwellTimes(1,i) + DwellTimes(2,i)) + MortA(1,i))* CumStageDISTwithDCO(2));
%    ModTimeToDiagnosis(3,i) = (1 - CumStageDISTwithDCO (3))/((1/(DwellTimes(1,i) + DwellTimes(2,i) + DwellTimes(3,i))+ MortA(1,i))* CumStageDISTwithDCO(3));
%    ModTimeToDiagnosis(4,i) = (1 - CumStageDISTwithDCO (4))/((1/(DwellTimes(1,i) + DwellTimes(2,i) + DwellTimes(3,i) + DwellTimes(4,i))+ MortA(1,i))* CumStageDISTwithDCO(4));
%    ModTimeToDiagnosis(5,i) = (1 - CumStageDISTwithDCO (5))/((1/(DwellTimes(1,i) + DwellTimes(2,i) + DwellTimes(3,i) + DwellTimes(4,i)+ DwellTimes(5,i))+ MortA(1,i))* CumStageDISTwithDCO(5));
% end 

for i = 1 : size(DwellTimes,2)
   ModTimeToDiagnosis(1,i) = (1 - CumStageDISTwithDCO (1))/((1/DwellTimes(1,i) + MortA(1,i)) * CumStageDISTwithDCO(1));
   ModTimeToDiagnosis(2,i) =(DwellTimes(1,i) + (1 - CumStageDISTwithDCO (2))/((1/ DwellTimes(2,i)) + MortA(1,i))* CumStageDISTwithDCO(2));
   ModTimeToDiagnosis(3,i) = (DwellTimes(1,i) + DwellTimes(2,i) +(1 - CumStageDISTwithDCO (3))/((1/ DwellTimes(3,i))+ MortA(1,i))* CumStageDISTwithDCO(3));
   ModTimeToDiagnosis(4,i) = (DwellTimes(1,i) + DwellTimes(2,i) + DwellTimes(3,i) +(1 - CumStageDISTwithDCO (4))/((1/ DwellTimes(4,i))+ MortA(1,i))* CumStageDISTwithDCO(4));
 end 
    
    %SojournTime and SojournTimePlusInsitu
    SojournTime =  zeros(1,size(DwellTimes,2));
    SojournTimePlusInsitu = zeros(1,size(DwellTimes,2));
    for j = 1: size(DwellTimes,2)
        SojournTime(1,j) = StageDist(1) * TimeToDiagnosis(1,j) + StageDist(2) * TimeToDiagnosis(2,j) + StageDist(3) * TimeToDiagnosis(3,j) +StageDist(4) * TimeToDiagnosis(4,j);
        SojournTimePlusInsitu(1,j) = SojournTime(1,j);

    end

    sumprobAgeA = zeros(1, size(AgeArrayLower,2));
    for j = 1: size(AgeArrayLower,2)
        for  probAgeInd = AgeArrayLower(j) : AgeArrayUpper(j)
            sumprobAgeA(1,j) = sumprobAgeA(1,j) + probAgeA(probAgeInd);
        end
    end

    %Calculating incidence by solving the differential equations 
    Incidence = IncidenceFunc(probAgeA,MortalityArray,prevalenceActual,cancerMortalityArray, AgeArrayLower,AgeArrayUpper);   
    %Incidence = Incidence *1000;
    Incidence =PAR{4};%GLOBOCAN breast cancer data has incidence not prevalence
    %Incidence = transpose(smooth(Incidence));

    %Calculating incidence by solving the differential equations MArkov Chain of the Markov Process 
    OnsetRates = PrevStageRates (StageDist /(1+DCO), Incidence,ModTimeToDiagnosis,probAgeA,MortalityArray,prevalenceActual,cancerMortalityArray, AgeArrayLower,AgeArrayUpper);
    % OnsetRates = smooth(OnsetRates);
    % OnsetRates = OnsetRates';
    % OnsetRates(1,1:4) =0;

    OnsetMatrix = zeros(1, size(AgeArrayLower,2)+ 1); 

    OnsetMatrix(1,1:size(AgeArrayLower,2))= OnsetRates * 1000;
    OnsetMatrix(1,size(AgeArrayLower,2)+1)= (OnsetRates * 1000) * transpose(sumprobAgeA);
     %sum((prevalenceActual) * (transpose((SojournTimeMat(i,:)) * transpose(sumprobAgeA)))) / sum(prevalenceActual)
     prevalenceTotal = prevalenceActual(1: size(AgeArrayLower,2))* transpose(sumprobAgeA) *1000;
    %OnsetRates(1,1:6)= 0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Interpolation to single age years, flatline after 80 years

    IncidenceF = interp1(AgeArrayMid,Incidence,[0:100],'PCHIP');
    %IncidenceF(80:end)=IncidenceF(80);

    OnsetRatesF = interp1(AgeArrayMid,OnsetRates,[0:100],'PCHIP');
    %OnsetRatesF(80:end)=OnsetRatesF(80);
    %OnsetRatesF(1:20) = 0;

    InsituOnset=sum(OnsetRatesF(1,15:80).*probAgeA(1,15:80)/(sum(probAgeA(1,15:80))));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [InsituOnsetRate,polypOnsetRateAgeGroup, PrevPolys]=PolypOnsetCal(OnsetRatesF,AgeArrayLower,AgeArrayUpper,PAR);
    
%% calculation of diagnostic rates starts here

    trial_run = 1;
    betaRuns = ones(4,size(AgeArrayLower,2),trial_run);

    for randRuns = 1:trial_run
        beta0 = ones(1,size(AgeArrayLower,2)) * rand()*4;
        initial_point=beta0(1,1);
        [diagr]=ExpSourj(beta0, StageDist,AgeArrayLower, AgeArrayUpper);

        %test shape by displaying
        % figure
        % hold on
        % X=0:100;
        % plot(X,diagr(:,1),'b')
        % plot(X,diagr(:,2),'r')
        % plot(X,diagr(:,3),'g')
        % plot(X,diagr(:,4),'c')
        % pause

        tic
        X=AgeArrayMid;
        %Y=[prevalenceActual*100];%if fit to prevalence
        Y=[Incidence*100];%if fit to incidence
        %[Z,diagr]=PopSim(beta0,X,PAR,IncidenceF,OnsetRatesF,DwellTimesF,AgeArrayLower,AgeArrayUpper,1)
        toc

        %Fit the prevalence model, find the diagnostic rates and plot the results
        tic
        [Z,PrevS,PrevC,diagrT,beta]=FitSoJourn(1000,beta0,X,Y,PAR,Incidence,OnsetRatesF,DwellTimesF,AgeArrayLower,AgeArrayUpper, Mu_Tx, BaseMortality);
        toc
        x = ExpSourjX(beta,AgeArrayMid, StageDist);
        betaRuns(:,:,randRuns)= (x').^(-1);
        initial_point_array(:,randRuns)=initial_point;

    end

    trial_detail=zeros(4,size(beta,2)+1,trial_run);
    trial_detail(1,1,:) = initial_point_array;
    trial_detail(:,2:(size(beta,2)+1),:) = betaRuns;
    for i = 1:trial_run
        trial_details(:,:,i) = transpose(trial_detail(:,:,i));
    end
    [DiagnosticRate]=ExpSourjX(beta,AgeArrayMid, StageDist);
%%
    %save results to be mapped to associaitons file
    RES{1}=PrevS;       %Prevalence of PreClinical Cancer
    RES{2}=PrevC;       %Prevalance of Clinical Cancer
    RES{3}=Incidence;   %GLOBOCAN incidence
    RES{4}=OnsetRates';%onsetRateSmooth'; Insitu onset rate from Markov Chain. we may not be writing this value to the template.
    RES{5}=polypOnsetRateAgeGroup;   %PolypOnset Rate for 3 kinds of polyps & also incidence from Healthy to Insitu (q1, q2, q3, q4 of flowchart) To be confirmed with prof.
    RES{6}=DiagnosticRate;
    RES{7}=DwellTimes';
    RES{8}=BaseMortalityC; %Relative Mortality without treatmet
    RES{9}=MortA';% Disease free mortality
    RES{10} = Impact;% %Relative Mortality with treatment
    RES{11} = TxCoverage_A;% Current coverage of treatment
    RES{12} = PrevPolys; %Prevalance of Polyps

   

end%end of MarkvoChainCountry

function [PrevA,PrevS,PrevC,diagr,beta]=FitSoJourn(MaxFunEvals,beta0,X,Y,PAR,Incidence,OnsetRatesF,DwellTimesF,AgeArrayLower,AgeArrayUpper, Mu_Tx, BaseMortality,  InsituOnsetRate)
    LB = ones(1,size(AgeArrayLower,2))*0.01;
    UB = ones(1,size(AgeArrayLower,2))*4;
    
    LB(1,1:3)=0;
    UB(1,1:3)=0;
    %LB=[-10,-1, -10,-10];
    %UB=[100,0.5, 10, 100];%if exponential curve for non-US
    %UB = [2000, 0.5,10,100];%% if logistic US
    %UB = [10000, 0.5,10,1000];%% if logistic non-US

    %OPTIONS = optimset('largescale','off','LevenbergMarquardt','on','display','off','Diagnostics','on','MaxFunEvals',700,'TolFun',1e-8,'TolCon',1e-8);
    OPTIONS = optimset('MaxFunEvals',MaxFunEvals,'Diagnostics','on','TolFun',1e-6,'TolCon',1e-6);

%     [beta,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@(x,X)...
%     PopSim(x,X,PAR,Incidence,OnsetRatesF,DwellTimesF,AgeArrayLower,AgeArrayUpper,Mu_Tx, BaseMortality, 0),...
%     beta0,X,Y,LB,UB,OPTIONS);

    %FIT EACH AGE INDIVIDUALLY
    store_beta=beta0;
    for age = 1 : size(AgeArrayLower,2)
        individualX = X(1,age);

        [betafit_out,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@(x,individualX )...
        PopSim(x,individualX,PAR,Incidence,OnsetRatesF,DwellTimesF,AgeArrayLower,AgeArrayUpper,Mu_Tx, BaseMortality, 0, age,store_beta),...
        beta0(1,age),individualX,Y(1,age),LB(1,age),UB(1,age),OPTIONS);

        store_beta(1,age)=betafit_out;
    end
    beta=store_beta;

    % beta = [0,0,0,8.04350879264037,2.06419265914802,1.05243573362542,1.32322102488666,0.865267058630671,0.373372280464398,0.189700814075861,0.105110504312803,0.0734006218779887];

    %plot the results
    [PrevA,diagr,PrevS,PrevC]=PopSim(beta,X,PAR,Incidence,OnsetRatesF,DwellTimesF,AgeArrayLower,AgeArrayUpper,Mu_Tx, BaseMortality,1,1000000,store_beta);

end

function [diagr]=ExpSourj(b, StageDist, AgeArrayLower, AgeArrayUpper)
    diagr=[];
    X=0:100;
    %b(3)=-1;
    diagr1=ones(1,101);
    diagr2=ones(1,101);
    diagr3=ones(1,101);
    diagr4=ones(1,101);

    for i = 1:size(AgeArrayLower,2)
        value = (StageDist(1));
        if value == 0
            diagr1(AgeArrayLower(i):AgeArrayUpper(i)) = 100000000;
        else
         diagr1(AgeArrayLower(i):AgeArrayUpper(i))=(1/b(1,i))/value;
        end
        value = (StageDist(2)+StageDist(1));
        diagr2(AgeArrayLower(i):AgeArrayUpper(i))=(1/b(1,i))/value;
        value = (StageDist(3)+StageDist(2)+StageDist(1));
        diagr3(AgeArrayLower(i):AgeArrayUpper(i))=(1/b(1,i))/value;
        value = (StageDist(4)+StageDist(3)+StageDist(2)+StageDist(1));
        diagr4(AgeArrayLower(i):AgeArrayUpper(i))=(1/b(1,i))/value;
    end
    %%%US%%%%
    % for i = 1:size(AgeArrayLower,2)
    %     value = (StageDist(1));
    %     if value == 0
    %         diagr1(AgeArrayLower(i):AgeArrayUpper(i)) = 100000000;
    %     else
    %      diagr1(AgeArrayLower(i):AgeArrayUpper(i))=(1/b(1,i))/value;
    %     end
    %     value = (StageDist(2));
    %     diagr2(AgeArrayLower(i):AgeArrayUpper(i))=(1/b(1,i))/value;
    %     value = (StageDist(3));
    %     diagr3(AgeArrayLower(i):AgeArrayUpper(i))=(1/b(1,i))/value;
    %     value = (StageDist(4));
    %     diagr4(AgeArrayLower(i):AgeArrayUpper(i))=(1/b(1,i))/value;
    %     
    % end


    %%%%%%%%%%%%%%%%%
    %%%LOGISTIC CURVE%%%
    % %%%%%%%%%%%%%%%%%
    % diagr1=b(4)./(1+(b(1)*(exp(b(2)*X))));%*(exp(-b(3)*1));
    % diagr1(80:end)= diagr1(80);
    % 
    % value = (StageDist(1));
    % %value= value^2;
    % diagr2=(b(4)/value)./(1+(b(1)*(exp(-b(2)*X))));%*(exp(-b(3)*2));
    % diagr2(80:end)= diagr2(80);
    % 
    % value = (StageDist(2));%+StageDist(1));
    % %value= value^2;
    % diagr3=(b(4)/value)./(1+(b(1)*(exp(-b(2)*X))));%*(exp(-b(3)*3));
    % diagr3(80:end)= diagr3(80);
    % 
    % value = (StageDist(3));%+ StageDist(2)+ StageDist(1));
    % %value= value^2;
    % diagr4=(b(4)/value)./(1+(b(1)*(exp(-b(2)*X))));%*(exp(-b(3)*4));
    % diagr4(80:end)= diagr4(80);
    % % 

    %%%%%%%%%%%%%%%%%
    %%%EXPONENTIAL CURVE%%%
    %%%%%%%%%%%%%%%%%
    % diagr1=b(1)*(exp(b(2)*X));%*(exp(-b(3)*1));
    % diagr1(80:end)= diagr1(80);
    % 
    % value = (StageDist(1));
    % %value= value^2;
    % diagr2=(b(1)/value)*(exp(b(2)*X));%*(exp(-b(3)*2));
    % diagr2(80:end)= diagr2(80);
    % 
    % value = (StageDist(2)+StageDist(1));
    % %value= value^2;
    % diagr3=(b(1)/value)*(exp(b(2)*X));%*(exp(-b(3)*3));
    % diagr3(80:end)= diagr3(80);
    % 
    % value = (StageDist(3)+ StageDist(2)+ StageDist(1));
    % %value= value^2;
    % diagr4=(b(1)/value)*(exp(b(2)*X));%*(exp(-b(3)*4));
    % diagr4(80:end)= diagr4(80);

    diagr=[diagr1;diagr2;diagr3;diagr4]';
end

function [diagr]=ExpSourjX(b,X, StageDist);
    diagr=[];
    %b(3)=-0.3;

    %%%underdevloped countries where there is no screening diagnostic rate is likely to be highre in late stages than in early stages due to absence of screening 
    for i = 1: size(X,2)
        value = (StageDist(1));
        if value == 0
            diagr1(i) = 100000000;
        else
            diagr1(i)=(1/b(1,i))/value;
        end
    value = (StageDist(2)+StageDist(1));
    diagr2(i)=(1/b(1,i))/value;
    value = (StageDist(3)+StageDist(2)+StageDist(1));
    diagr3(i)=(1/b(1,i))/value;
    value = (StageDist(4)+StageDist(3)+StageDist(2)+StageDist(1));
    diagr4(i)=(1/b(1,i))/value;
    end

    %%%US
    % for i = 1: size(X,2)
    %     value = (StageDist(1));
    %     if value == 0
    %         diagr1(i) = 100000000;
    %     else
    %         diagr1(i)=(1/b(1,i))/value;
    %     end
    %     value = (StageDist(2));
    %     diagr2(i)=(1/b(1,i))/value;
    %     value = (StageDist(3));
    %     diagr3(i)=(1/b(1,i))/value;
    %     value = (StageDist(4));
    %     diagr4(i)=(1/b(1,i))/value;
    % end

    %%%%%%%%%%%%%%%%%
    %%%LOGISTIC CURVE%%%
    % %%%%%%%%%%%%%%%%%
    % diagr1=b(4)./(1+(b(1)*(exp(b(2)*X))));%*(exp(-b(3)*1));
    % diagr1(end)= diagr1(end-1);
    % 
    % value = (StageDist(1));
    % %value= value^2;
    % diagr2=(b(4)/value)./(1+(b(1)*(exp(-b(2)*X))));%*(exp(-b(3)*2));
    % diagr2(end)= diagr2(end-1);
    % 
    % value = (StageDist(2));%+StageDist(1));
    % %value= value^2;
    % diagr3=(b(4)/value)./(1+(b(1)*(exp(-b(2)*X))));%*(exp(-b(3)*3));
    % diagr3(end)= diagr3(end-1);
    % 
    % value = (StageDist(3));%+ StageDist(2)+ StageDist(1));
    % %value= value^2;
    % diagr4=(b(4)/value)./(1+(b(1)*(exp(-b(2)*X))));%*(exp(-b(3)*4));
    % diagr4(end)= diagr4(end-1);

    %%%%%%%%%%%%%%%%%
    %%%EXPONENTIAL CURVE%%%
    %%%%%%%%%%%%%%%%%
    % diagr1=b(1)*(exp(b(2)*X));%*(exp(-b(3)*1));
    % diagr1(end)= diagr1(end-1);
    % 
    % value = (StageDist(1));
    % %value= value^2;
    % diagr2=(b(1)/value)*(exp(b(2)*X));%*(exp(-b(3)*2));
    % diagr2(end)= diagr2(end-1);
    % 
    % value = (StageDist(2)+StageDist(1));
    % %value= value^2;
    % diagr3=(b(1)/value)*(exp(b(2)*X));%*(exp(-b(3)*3));
    % diagr3(end)= diagr3(end-1);
    % 
    % value = (StageDist(3)+ StageDist(2)+ StageDist(1));
    % %value= value^2;
    % diagr4=(b(1)/value)*(exp(b(2)*X));%*(exp(-b(3)*4));
    % diagr4(end)= diagr4(end-1);

    diagr=[diagr1;diagr2;diagr3;diagr4]';
    %diagr=max(diagr,0.5);

end

function [Z,diagr,PrevS,PrevC]=PopSim(fitPAR,X,PAR,Incidence,OnsetRatesF,DwellTimesF,AgeArrayLower,AgeArrayUpper, Mu_Tx, BaseMortality, PlotModel, age, store_beta)

    Z=[];
    PrevS=[];
    PrevC=[];

    AgeArrayMid=(AgeArrayUpper+AgeArrayLower)/2;

    StageDist=PAR{6};
    
    if age < 13
        if length(fitPAR)<2
            current_beta = fitPAR;
            fitPAR = store_beta;
            fitPAR(1,age) = current_beta;
        end
    end  

    [diagrT]=ExpSourj(fitPAR,StageDist, AgeArrayLower, AgeArrayUpper);

    diagr=1./diagrT;
    
    if StageDist(1) == 0
        diagr(:,1) = 0;
        diagrT(:,1) = 0;
    end
  
    PrevA=zeros(size(AgeArrayLower));
    %IncidenceF=IncidenceF';
    OnsetRatesF=OnsetRatesF';
    lambda=(1./DwellTimesF)';

    probAgeA=PAR{1};probAgeA(81:end)=probAgeA(81)/length(probAgeA(81:end));
    mu=PAR{2}';mu(81:end)=mu(80)+1/5;
    BirthsF=PAR{3};%female pop only
    prev=PAR{4};
    prevN=prev/sum(prev);
    muC =smooth(PAR{5}');
    StageDist=PAR{6};

    agemax=101;
    m_a=length(probAgeA);
    %age,state, duration of disease
    pop0=100000;%work with pop of standard size

    POP=zeros(m_a,9);
    POP(1:m_a,1)=probAgeA*pop0;

    t_max=500;
    CFR = zeros(m_a,t_max);


    dt=1/15;
    for t=1:t_max

        %fast time loop
        Diags=zeros(length(POP),4,t_max);
        IncidenceE = zeros(101,1);
        death_pre_clinical_4 = zeros(101,1);
        for t1=1:1/dt    
            %healhy pop    
            POP(:,1)=POP(:,1) + dt*(-mu-OnsetRatesF).*POP(:,1);

            %non-diagnosed cancer
            %POP(:,2)=POP(:,2)+dt*(-mu-lambda(:,1)).*POP(:,2) + dt*(OnsetRatesF).*POP(:,1);
            POP(:,2)=POP(:,2)+dt*(-mu-lambda(:,1)-diagr(:,1)).*POP(:,2) + dt*(OnsetRatesF).*POP(:,1);
            POP(:,3)=POP(:,3)+dt*(-mu-lambda(:,2)-diagr(:,2)).*POP(:,3) + dt*(lambda(:,1)).*POP(:,2);
            POP(:,4)=POP(:,4)+dt*(-mu-lambda(:,3)-diagr(:,3)).*POP(:,4) + dt*(lambda(:,2)).*POP(:,3);
            POP(:,5)=POP(:,5)+dt*(-mu-lambda(:,4)-diagr(:,4)).*POP(:,5) + dt*(lambda(:,3)).*POP(:,4);
            
            % Deaths in preclinical last stage
            death_pre_clinical_4(:,1) = death_pre_clinical_4(:,1) + POP(:,5).*(lambda(:,4))*dt;

            %diagnosed cancer
            IncidenceE(:,1) = IncidenceE(:,1) + dt*(diagr(:,1)).*POP(:,2);
            POP(:,6)=POP(:,6)+dt*(-mu ).*POP(:,6) - dt*BaseMortality(1,:)'.*POP(:,6) + dt*(diagr(:,1)).*POP(:,2);
            CFR(:,t)= CFR(:,t)+ dt*(mu ).*POP(:,6)+ dt*BaseMortality(1,:)'.*POP(:,6);

            IncidenceE(:,1) = IncidenceE(:,1) + dt*(diagr(:,2)).*POP(:,3);
            POP(:,7)=POP(:,7)+dt*(-mu ).*POP(:,7) - dt*BaseMortality(2,:)'.*POP(:,7)+ dt*(diagr(:,2)).*POP(:,3);
            CFR(:,t)= CFR(:,t)+ dt*(mu ).*POP(:,7)+ dt*BaseMortality(2,:)'.*POP(:,7);

            IncidenceE(:,1) = IncidenceE(:,1) + dt*(diagr(:,3)).*POP(:,4);
            POP(:,8)=POP(:,8)+dt*(-mu ).*POP(:,8) - dt*BaseMortality(3,:)'.*POP(:,8)+ dt*(diagr(:,3)).*POP(:,4);
            CFR(:,t)= CFR(:,t)+dt*(mu ).*POP(:,8) + dt*BaseMortality(3,:)'.*POP(:,8);

            IncidenceE(:,1) = IncidenceE(:,1) + dt*(diagr(:,4)).*POP(:,5);
            POP(:,9)=POP(:,9)+dt*(-mu ).*POP(:,9) - dt*BaseMortality(4,:)'.*POP(:,9)+ dt*(diagr(:,4)).*POP(:,5);
            CFR(:,t)= CFR(:,t)+dt*(mu ).*POP(:,9) + dt*BaseMortality(4,:)'.*POP(:,9);
            
            %if t > t_max - 100
            Diags(:,1,t)=Diags(:,1,t)+dt*( ((diagr(:,1)).*POP(:,2)) );
            Diags(:,2,t)=Diags(:,2,t)+dt*( ((diagr(:,2)).*POP(:,3)) );
            Diags(:,3,t)=Diags(:,3,t)+dt*( ((diagr(:,3)).*POP(:,4)) );
            Diags(:,4,t)=Diags(:,4,t)+dt*( ((diagr(:,4)).*POP(:,5)) );
            %end

        end

        DiagsAll=sum(Diags,2);

        %aging
        for I=1:9
            POP(2:end,I)=POP(1:end-1,I);
        end

        %must decide if births are to be based on population size
        Births=BirthsF*10000;%BirthsF*sum(POP(:));
        POP(1,1)=Births;

        DiagsA(1,t) = sum(Diags(:,1,t))/ sum(sum(Diags(:,:,t)));
        DiagsA(2,t) = sum(Diags(:,2,t))/ sum(sum(Diags(:,:,t)));
        DiagsA(3,t) = sum(Diags(:,3,t))/ sum(sum(Diags(:,:,t)));
        DiagsA(4,t) = sum(Diags(:,4,t))/ sum(sum(Diags(:,:,t)));
    end
    %prevalence of cancer by age
    %average mort in each age group
    for J=1:length(AgeArrayLower)
        indA=AgeArrayLower(J):AgeArrayUpper(J);

        inc=IncidenceE(indA,1); 

        p_a=POP(indA,:);
        p_c=POP(indA,[6:9]);

        %pre clinical
        p_s1=POP(indA,[2]);
        p_s2=POP(indA,[3]);
        p_s3=POP(indA,[4]);
        p_s4=POP(indA,[5]);

        %diagnosed/clinical
        p_c1=POP(indA,[6]);
        p_c2=POP(indA,[7]);
        p_c3=POP(indA,[8]);
        p_c4=POP(indA,[9]);

        incidenceCalc(J) = sum(inc(:))/sum(p_a(:));
        PrevA(J)=sum(p_c(:))/sum(p_a(:));

        %prev cancer pre clinical
        PrevS(J,1)=sum(p_s1(:))/sum(p_a(:));
        PrevS(J,2)=sum(p_s2(:))/sum(p_a(:));
        PrevS(J,3)=sum(p_s3(:))/sum(p_a(:));
        PrevS(J,4)=sum(p_s4(:))/sum(p_a(:));

        %prev cancer clinical
        PrevC(J,1)=sum(p_c1(:))/sum(p_a(:));
        PrevC(J,2)=sum(p_c2(:))/sum(p_a(:));
        PrevC(J,3)=sum(p_c3(:))/sum(p_a(:));
        PrevC(J,4)=sum(p_c4(:))/sum(p_a(:));

    
    end

    %Diags=Diags/sum(Diags)*100;
    %Z=[PrevA*100];If fitting to prevalence
    if age >1000
        Z=incidenceCalc;
        DCO =sum(death_pre_clinical_4) / (sum(death_pre_clinical_4)+sum(sum(IncidenceE)));
    else
        Z=incidenceCalc(age)*100;
        DCO =sum(death_pre_clinical_4) / (sum(death_pre_clinical_4)+sum(sum(IncidenceE)));
    end
    %PlotModel=1;
    if(PlotModel==1)

        close all
        figure
        subplot(4,1,1)
        plot(POP(:,1))
        title('total pop')

        subplot(4,1,2)
        plot(sum(POP(:,[2:5]),2) )
        title('cancer not diagnosed')

        subplot(4,1,3)
        plot(sum(POP(:,[6:9]),2) )
        title('total pop')
        title('cancer diagnosed')

        subplot(4,1,4)
        plot(OnsetRatesF )
        title('total pop')
        title('Onset Rate')


        %PrevA=Z(1:length(prev))
        figure
        hold on
        plot(AgeArrayMid,prev*100,'k*')
        plot(AgeArrayMid,PrevA,'r-')
        title('Prevalence of Diagnosed cases by age')

        figure
        hold on
        plot(AgeArrayMid,Incidence*100,'k*')
        plot(AgeArrayMid,incidenceCalc*100,'r-')
        title('Incidence of Diagnosed cases by age')

        figure
        hold on
        plot(0:100,OnsetRatesF )

        title('Onset Rate')

        figure
        hold on
        plot(0:100,CFR(:,t))
        plot(0:100,CFR(:,t-100))
        plot(0:100,CFR(:,t-200))
        plot(0:100,CFR(:,t-300))
        title('CFR')



        figure
        hold on
        X=0:100;
        plot(X,diagrT(:,1),'b')
        plot(X,diagrT(:,2),'r')
        plot(X,diagrT(:,3),'g')
        plot(X,diagrT(:,4),'c')
        legend('in-situ', 'local', 'regional', 'distant');
        title('Time to diagnosis')

        figure
        hold on
        X=0:100;
        plot(X,DwellTimesF(1,:),'b')
        plot(X,DwellTimesF(2,:),'r')
        plot(X,DwellTimesF(3,:),'g')
        plot(X,DwellTimesF(4,:),'c')
        legend('in-situ', 'local', 'regional', 'distant');

        title('Dwell Time')

        figure
        hold on
        plot(1:4,DiagsA(1:4,t),'r-')
        plot(1:4,StageDist(1:4),'k*')
        title('Stage at diagnosis')


        figure
        hold on
        subplot(4,1,1)
        plot(1:500,DiagsA(2,:),'r-')
        title('Local')

        subplot(4,1,2)
        plot(1:500,DiagsA(3,:),'r-')
        title('Regional')

        subplot(4,1,3)
        plot(1:500,DiagsA(4,:),'r-')
        title('Distant')

    end

end


%%ESTIMATE INCIDENCE
function [incidenceArray] = IncidenceFunc(probAgeA,MortalityArray,prevalenceActual,cancerMortalityArray, AgeArrayLower,AgeArrayUpper)

    incidenceArray = zeros(1, size(AgeArrayLower,2));
    sumAgeArray= zeros(1, size(AgeArrayLower,2));
    for j = 1 : size (AgeArrayLower,2)
        sumAge = 0;
        mortSum = 0;
        for  i = AgeArrayLower(j) : AgeArrayUpper(j)
            sumAge = sumAge + probAgeA(i);
            mortSum = mortSum  + cancerMortalityArray(i);
        end

        sumAgeArray(j) = sumAge;
        mortSum = mortSum / (AgeArrayUpper(j) - AgeArrayLower(j) + 1);

        if j == 1
            incidenceArray(j) = (prevalenceActual(j) * sumAgeArray(j) * (mortSum + (1 / (AgeArrayUpper(j) - AgeArrayLower(j) + 1))))/ (sumAgeArray(j) - prevalenceActual(j));
            % incidenceArray(j) = max(0,(prevalenceActual(j) * sumAgeArray(j) * (mortSum + (1 / (AgeArrayUpper(j) - AgeArrayLower(j) + 1))))/ (sumAgeArray(j) - prevalenceActual(j)));
        else
            %incidenceArray(j) = max(incidenceArray(j-1),((prevalenceActual(j) * sumAgeArray(j) * (mortSum + (1 / (AgeArrayUpper(j) - AgeArrayLower(j) + 1)))) -  (prevalenceActual(j - 1) * sumAgeArray(j-1) * (1 / (AgeArrayUpper(j - 1) - AgeArrayLower(j -1) + 1))))/ (sumAgeArray(j) - prevalenceActual(j)));
             incidenceArray(j) = ((prevalenceActual(j) * sumAgeArray(j) * (mortSum + (1 / (AgeArrayUpper(j) - AgeArrayLower(j) + 1)))) -  (prevalenceActual(j - 1) * sumAgeArray(j-1) * (1 / (AgeArrayUpper(j - 1) - AgeArrayLower(j -1) + 1))))/ (sumAgeArray(j) - prevalenceActual(j));
            if incidenceArray(j) < 0
             incidenceArray(j) = incidenceArray(j-1);
            end

        end%if J==1

     end
end%IncidenceFunc

%%ESTIMATE Previous Stage 
function [OnsetRates] = PrevStageRates(StageDist,CurrentStageRate, TimeToDiagnosis, probAgeA, MortalityArray,prevalenceActual,cancerMortalityArray, AgeArrayLower, AgeArrayUpper)

% % input list for model 1
% % mu = Mortalityarray
% % lamda = lamdaarray
% % I_D_a = current stage rates
% % C_a = A_a = ProbAgeA
% % agearraylower and upper are lower and upper bounds of the age groups
% 
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
% %%                    
%                     % For analysis by age group %     
%     
%     % Here, all the matrices should in dimension of 1 X 12
%     
%     % Dimension of matrix "probAgeA" is 1 X 100, now converting this matrix
%     % into 1 X 12 matrix. This has to done in accordance with predesigned
%     % age groups, "AgeArrayLower" and "AgeArrayUpper", these
%     % two matrices define age groups. With help of these matrices, elements
%     % of probAgeA are added together (Simple "OR" operation of
%     % probability). By doing this we'll get prob of a person being in that
%     % group.
% 
%     a_max_group=size(AgeArrayLower,2);
%     
%     % defiing matrix for storing values
%     C_a=zeros(1,size(AgeArrayLower,2));
%     
%     for i=(1:size(AgeArrayLower,2))
%         
%         % Reading values of lower and upper limit of "i th" age group
%         l_limit=AgeArrayLower(1,i);
%         u_limit=AgeArrayUpper(1,i);
%         
%         % with help l_limit and u_limit, we'll do summation over those
%         % indexes in probAgeA
%        C_a(1,i)=sum(probAgeA(1,l_limit:u_limit));
%         
%     end
%     
%     A_a=C_a;
%                                     % here, we have made an assumption that
%                                     % A_a = C_a, because probability of a
%                                     % person being diagnosed in age "a" is
%                                     % very low
%     
%     % Next matrix which is not in 1 X 12 dimension is MortalityArray. this
%     % matrix contains rates and not prbability values, hence we cannot
%     % simply add them together. We have to calculate expected value first
%     % and then add them together.
%     
%     % calculating expected value of mortality rate
%     mu_exp=probAgeA.*MortalityArray;
%     
%     % defining matrix for mu
%     mu=zeros(1,size(AgeArrayLower,2));
%     
%     % converting mu_exp in 1 X 12 matrix
%     for i=(1:size(AgeArrayLower,2))
%         
%         % Reading values of lower and upper limit of "i th" age group
%         l_limit=AgeArrayLower(1,i);
%         u_limit=AgeArrayUpper(1,i);
%         
%         % with help l_limit and u_limit, we'll do summation over those
%         % indexes in mu_exp
%         mu(1,i)=sum(mu_exp(1,l_limit:u_limit))/sum(probAgeA(l_limit:u_limit));
%         
%     end
%     
%     % defining lamda
%     lamda=zeros(4,size(AgeArrayLower,2));
%     for i=(1:4)
%         for n=(1:12)
%             lamda(i,n)=1/TimeToDiagnosis(i,n);
%         end
%     end
%     
%     % defining matix for long term prob of being in Healthy state and being
%     % in Undiagnosed state and prob of transition from H to U
%     steady_h=zeros(1,size(AgeArrayLower,2));
%     steady_h(1,1:3)=A_a(1,1:3);
%   
%     steady_u=zeros(1,size(AgeArrayLower,2));
%     probability_hu=zeros(1,size(AgeArrayLower,2));
%     
%     % defining incidence matrix
%     I_D_a=CurrentStageRate;
%     I_D_a(1,1:3)=0;
%     
%     % creating array for storing mid age values
%     mid_age_group=zeros(1,a_max_group);
%     for i=(1:12)
%         mid_age_group(1,i)=(AgeArrayLower(1,i)+AgeArrayUpper(1,i))/2;
%     end
%     
%     % defining prob of being in one of 4 cancer states
%     s_i=StageDist;
%     
%     % onset rate matrix
%     OnsetRates=zeros(1,a_max_group);
%     
%     %% Plotting some parameters
%     
%     figure
%     x=(1:12);
%     plot(x,I_D_a);
%     ylabel('Incidence');
%     figure
%     plot(x,A_a);
%     ylabel('A_a/C_a');
%     plot(x,mu);
%     ylabel('Mortalities');
%     figure
%     for i=(1:4)
%         plot(x,lamda(i,:));
%         hold on
%     end
%     ylabel('lamda');
%     figure
%     for i=(1:4)
%         plot(x,TimeToDiagnosis(i,:));
%         hold on
%     end
%     ylabel('Time to Diagnosis');
%     
%     %%
%     
%     for a = 4:12   % START ALGORITHM 
%       
%     %%STEP 1
%         IC = I_D_a(1,a) * C_a(1,a);
%         
%         sumKloop =0;
%         for k = 3:a-1
%                       
%            sumIloopNumerator = 0;
%            t = mid_age_group(1,a)-mid_age_group(1,k);
%             for i = 1:4
%                 sumIloopNumerator = sumIloopNumerator + s_i(1,i)* (  (1- exp(-t*lamda(i,k))) - (1- exp(-(t-1)*lamda(i,k))) );
%             end
%             muProduct=1;
%             for j = k:a
%              muProduct = muProduct * exp(-(mid_age_group(j)-mid_age_group(j-1)) *mu(1,j));
%             end
%             sumKloop = sumKloop + steady_h(1,k) * probability_hu(1,k) * sumIloopNumerator * muProduct;
%         end
%         
%         
%         sumIloopDenominator = 0;
%         t = (AgeArrayUpper(1,a)-AgeArrayLower(1,a))/2;
%         for i = 1:4
%                 sumIloopDenominator = sumIloopDenominator + s_i(1,i)* (1- exp(-t*lamda(i,a)));
%         end
%         denominator = A_a(1,a) * sumIloopDenominator * exp(-t * mu(1,a)) - IC;
%         probability_hu(1,a) = (IC - sumKloop)/ denominator;
%         
%     %% STEP 2: does not look like we use this anywhere
%     
%     
%     %% STEP 3
%      sumKloop =0;
%         for k = 3:a-1
%                       
%            sumIloopNumerator = 0;
%            t = mid_age_group(1,a)-mid_age_group(1,k);
%             for i = 1:4
%                 sumIloopNumerator = sumIloopNumerator + s_i(1,i)*  exp(-t*lamda(i,k)) ;
%             end
%             muProduct=1;
%             for j = k:a
%              muProduct = muProduct * exp(-(mid_age_group(j)-mid_age_group(j-1)) *mu(1,j));
%             end
%             sumKloop = sumKloop + steady_h(1,k) * probability_hu(1,k) * sumIloopNumerator * muProduct;
%         end
%         
%       steady_h(1,a) = (A_a(1,a) - sumKloop) / (1+probability_hu(1,a));
%     
%     
%     
%     end %% END ALGORITHM
%    
    


numStage = size(TimeToDiagnosis,1);
numAgeGroups = size(AgeArrayLower,2);
minAgeCancer = 3;

%CurrentStageRate = incidence =(I_(D_a ))
%lamdaArray= lambda
%probAgeA = age distribution =(C_a in step 2 of algo)
%MortalityArray= mu
%AgeArrayLower, AgeArrayUpper =lower and upper bounds of age group

%prevalenceActual= do not need this anymore, but lets not delete it as yet 
%cancerMortalityArray = do not need this anymore, but lets not delete it as yet 

    
%%                    
                    % For analysis by age group %     
    
    % Here, all the matrices should in dimension of 1 X 12
    
    % Dimension of matrix "probAgeA" is 1 X 100, now converting this matrix
    % into 1 X 12 matrix. This has to done in accordance with predesigned
    % age groups, "AgeArrayLower" and "AgeArrayUpper", these
    % two matrices define age groups. With help of these matrices, elements
    % of probAgeA are added together (Simple "OR" operation of
    % probability). By doing this we'll get prob of a person being in that
    % group.

    a_max_group=size(AgeArrayLower,2);
    
    % defiing matrix for storing values
    C_a=zeros(1,size(AgeArrayLower,2));
    
    for i=(1:size(AgeArrayLower,2))
        
        % Reading values of lower and upper limit of "i th" age group
        l_limit=AgeArrayLower(1,i);
        u_limit=AgeArrayUpper(1,i);
        
        % with help l_limit and u_limit, we'll do summation over those
        % indexes in probAgeA
       C_a(1,i)=sum(probAgeA(1,l_limit:u_limit));
        
    end
    
    A_a=C_a;
                                    % here, we have made an assumption that
                                    % A_a = C_a, because probability of a
                                    % person being diagnosed in age "a" is
                                    % very low
    
    % Next matrix which is not in 1 X 12 dimension is MortalityArray. this
    % matrix contains rates and not prbability values, hence we cannot
    % simply add them together. We have to calculate expected value first
    % and then add them together.
    
    % calculating expected value of mortality rate
    mu_exp=probAgeA.*MortalityArray;
    
    % defining matrix for storing values
    mu=zeros(1,size(AgeArrayLower,2));
    
  
    
    % converting mu_exp in 1 X 12 matrix
    for i=(1:size(AgeArrayLower,2))
        
        % Reading values of lower and upper limit of "i th" age group
        l_limit=AgeArrayLower(1,i);
        u_limit=AgeArrayUpper(1,i);
        
        % with help l_limit and u_limit, we'll do summation over those
        % indexes in mu_exp
        mu(1,i)=sum(mu_exp(1,l_limit:u_limit))/sum(probAgeA(l_limit:u_limit));
        
    end
    
    % defining lamda
    lamda=zeros(numStage,size(AgeArrayLower,2));
    for i=(1:numStage)
        for n=(1:numAgeGroups)
            lamda(i,n)=1/TimeToDiagnosis(i,n);
        end
    end
    
    % defining matix for long term prob of being in Healthy state and being
    % in Undiagnosed state and prob of transition from H to U
  steady_h=zeros(1,size(AgeArrayLower,2));
  steady_h(1,1:minAgeCancer-1)=A_a(1,1:minAgeCancer-1);
  
   steady_u=zeros(1,size(AgeArrayLower,2));
  probability_hu=zeros(1,size(AgeArrayLower,2));
    
    % defining incidence matrix
    I_D_a=CurrentStageRate;
    I_D_a(1,1:minAgeCancer-1)=0;
    
    % creating array for storing mid age values
    mid_age_group=zeros(1,a_max_group);
    for i=(1:numAgeGroups)
        mid_age_group(1,i)=(AgeArrayLower(1,i)+AgeArrayUpper(1,i))/2;
    end
    
    % defining prob of being in one of 4 cancer states
  s_i=StageDist;
    
    % onset rate matrix
    OnsetRates=zeros(1,a_max_group);
  %% Plotting some parameters
    
    figure
    x=(1:numAgeGroups);
    plot(x,I_D_a);
    ylabel('Incidence');
    figure
    plot(x,A_a);
    ylabel('A_a/C_a');
    plot(x,mu);
    ylabel('Mortalities');
    figure
    for i=(1:numStage)
        plot(x,lamda(i,:));
        hold on
    end
    ylabel('lamda');
    figure
    for i=(1:numStage)
        plot(x,TimeToDiagnosis(i,:));
        hold on
    end
    ylabel('Time to Diagnosis');
    
    for a = minAgeCancer+1:numAgeGroups   % START ALGORITHM 
      
   %%STEP 1
        IC = I_D_a(1,a) * C_a(1,a);
        
        sumKloop =0;
        for k = minAgeCancer:a-1
                      
           sumIloopNumerator = 0;
%           t = mid_age_group(1,a)-mid_age_group(1,k);
%            for i = 1:4
%                 sumIloopNumerator = sumIloopNumerator + s_i(1,i)* (  (1- exp(-t*lamda(i,k))) - (1- exp(-(t-1)*lamda(i,k))) );
%            end

           t = (AgeArrayUpper(1,a)-mid_age_group(1,k));
           t_minus_1 = (AgeArrayLower(1,a)-mid_age_group(1,k));
            for i = 1:numStage
                sumIloopNumerator = sumIloopNumerator + s_i(1,i)* (  (1- exp(-t*lamda(i,k))) - (1- exp(-(t_minus_1)*lamda(i,k))) );
            end
        
            muProduct=1;
%             for j = k:a
%              muProduct = muProduct * exp(-(mid_age_group(j)-mid_age_group(j-1)) *mu(1,j));
%             end
            for j = round(mid_age_group(1,k)):AgeArrayUpper(1,a)
             muProduct = muProduct * exp(-MortalityArray(1,j));
            end
            sumKloop = sumKloop + steady_h(1,k) * probability_hu(1,k) * sumIloopNumerator * muProduct;
        end
        
        
        sumIloopDenominator = 0;
        t = (AgeArrayUpper(1,a)-AgeArrayLower(1,a))/2;
        for i = 1:numStage
                sumIloopDenominator = sumIloopDenominator + s_i(1,i)* (1- exp(-t*lamda(i,a)));
        end
        denominator = A_a(1,a) * sumIloopDenominator * exp(-t * mu(1,a)) - IC;
        probability_hu(1,a) = (IC - sumKloop)/ denominator;
    %% STEP 2: does not look like we use this anywhere
    
    
    %% STEP 3
     sumKloop =0;
        for k = minAgeCancer:a-1
                      
           sumIloopNumerator = 0;
%           t = mid_age_group(1,a)-mid_age_group(1,k);
             t = AgeArrayUpper(1,a)-mid_age_group(1,k);
            for i = 1:numStage
                sumIloopNumerator = sumIloopNumerator + s_i(1,i)*  exp(-t*lamda(i,k)) ;
            end
       
            muProduct=1;
%             for j = k:a
%              muProduct = muProduct * exp(-(mid_age_group(j)-mid_age_group(j-1)) *mu(1,j));
%             end
            for j = round(mid_age_group(1,k)):AgeArrayUpper(1,a)
             muProduct = muProduct * exp(-MortalityArray(1,j));
            end
            sumKloop = sumKloop + steady_h(1,k) * probability_hu(1,k) * sumIloopNumerator * muProduct;
        end
        
      steady_h(1,a) = (A_a(1,a) - sumKloop) / (1+probability_hu(1,a));
    
    
    
    end %% END ALGORITHM
 
   
    OnsetRates = -log(1- probability_hu);

    onset_ref=[0	0	0	0	4.45E-05	7.77E-05	8.36E-05	0.00020033	0.000463829	0.000835422	0.001020914	0.000977113];
    
    x=(1:101);
    figure
    plot(mid_age_group(1:12),onset_ref,'--',mid_age_group,OnsetRates);
    legend('Reference values of Onset','Estimated');
    title('Colorectal Cancer, AFRE');

    y=1;

end %PrevStageRates


% apply stage at diagnosis to incidence to get persons by stage
% Use dwell times to estimate onset rate in Insitu
