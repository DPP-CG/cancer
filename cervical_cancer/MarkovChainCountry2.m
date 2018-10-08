
function [RES,RES2]=MarkovChainCountry2(PAR)
distHPV = PAR{10};
RES=[];

close all
nsims=10000;
%Country-specific parameters
probAgeA=PAR{1};
MortalityArray=PAR{2};
for i = 80:101
    MortalityArray(i)=MortalityArray(i-1)* MortalityArray(79)/MortalityArray(78);
end
Births=PAR{3};
prevalenceActualF=PAR{4};
%prevN=prevalenceActual/sum(prevalenceActual);

%Prevalence and case fatality rate: Age Group is different from Spectrum age group so converting
AgeArrayL = PAR{7};
AgeArrayU = PAR{8};
AgeArrayM=(AgeArrayU+AgeArrayL)/2;

CFR = interp1(AgeArrayM,PAR{5},[0:100],'PCHIP');

cancerMortalityArray =-log(1-min(.99,CFR));
%cancerMortalityArray(1,80:end) = cancerMortalityArray(1,79);
%cancerMortalityArray = max(cancerMortalityArray, MortalityArray);
PAR{5} = cancerMortalityArray;
StageDist=PAR{6};

%Spectrum Age groups
AgeArrayLower=[1    5	10	15  20	25	30	40	50	60	70	80];
AgeArrayUpper=[4	9	14  19 	24	29	39	49	59	69	79	100];
AgeArrayMid=((AgeArrayUpper+AgeArrayLower)/2);

%Interpolation to single age years, flatline after 80 years
PrevalenceF = interp1(AgeArrayM,prevalenceActualF,[0:100],'PCHIP');
%PrevalenceF(80:end)=PrevalenceF(80);
PrevalenceF= max(0,PrevalenceF);

prevalenceActual = zeros(size(prevalenceActualF,1),size(AgeArrayLower,2));

for i = 1 : size(prevalenceActual,2)
    indA=AgeArrayLower(i):AgeArrayUpper(i);
    prevalenceActual(1,i)= dot(PrevalenceF(indA),probAgeA(indA))/sum(probAgeA(indA)) ;%PrevalenceF(AgeArrayMid(i));%
end
PAR{4}=prevalenceActual;

%Progression times in undiagnoses cases (Goldie et al., 2003)
DwellTimesMidAge =  [1 29 54 58 70];
%
DwellTimes(2,:) = 1./[.31 .31 .31 .31 .31];%1./[0.198450939 0.198450939 0.198450939 0.198450939];%
DwellTimes(3,:)= 1./[.332 .332 .332 .332 .332];%1./[0.105360516 0.105360516 0.105360516 0.105360516];%
DwellTimes(4,:)= 1./[.485 .485 .485 .485 .485];%1./[0.356674944 0.356674944 0.356674944 0.356674944];%
DwellTimes(5,:)= 1./[0.7 0.7 0.7 0.7 0.7];


%DwellTimes(1,:)= 1./(0.0302*exp(0.0601*DwellTimesMidAge));%[0.01 0.273 1.185 5.290];
DwellTimes(1,:)= 1./[0.03 0.03 0.273 1.185 5.290];%[0.01 0.04 .04 .04];

% DwellTimes(2,:) = 1./[0.198450939 0.198450939 0.198450939 0.198450939];%
% DwellTimes(3,:)= 1./[0.105360516 0.105360516 0.105360516 0.105360516];%
% DwellTimes(4,:)= 1./[0.356674944 0.356674944 0.356674944 0.356674944];%
% DwellTimes(5,:)= 1./[0.4 0.4 0.4 0.4];
%
%
% DwellTimes(1,:)= 1./[0.051293294 0.051293294 0.051293294 0.051293294];%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Interpolation to single age years, flatline after 80 years

DwellTimes1 = interp1(DwellTimesMidAge,DwellTimes(1,:),[0:100],'PCHIP');
% DwellTimes1(80:end)=DwellTimes1(80);

DwellTimes2 = interp1(DwellTimesMidAge,DwellTimes(2,:),[0:100],'PCHIP');
% DwellTimes2(80:end)=DwellTimes2(80);

DwellTimes3 = interp1(DwellTimesMidAge,DwellTimes(3,:),[0:100],'PCHIP');
% DwellTimes3(80:end)=DwellTimes3(80);

DwellTimes4 = interp1(DwellTimesMidAge,DwellTimes(4,:),[0:100],'PCHIP');
% DwellTimes4(80:end)=DwellTimes4(80);

DwellTimes5 = interp1(DwellTimesMidAge,DwellTimes(5,:),[0:100],'PCHIP');
% DwellTimes5(80:end)=DwellTimes5(80);

DwellTimesF=[DwellTimes1;DwellTimes2;DwellTimes3;DwellTimes4;DwellTimes5];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculates 10-year relative survival with and without treatment,
[Mu_Tx, ImpactAge, BaseMortality, TxCoverage, StgDist, RelativeMortTx, RelativeMortNoTx] = Mortality(DwellTimesF, MortalityArray, PrevalenceF, StageDist, smooth(PAR{5}),PAR);

%average values into age groups
MortA=[];
MortC=[];
for J=1:length(AgeArrayLower)
    indA=AgeArrayLower(J):AgeArrayUpper(J);
    MortA(J)=dot(MortalityArray(indA),probAgeA(indA))/sum(probAgeA(indA));
    MortC(J)=dot(cancerMortalityArray(indA),probAgeA(indA))/sum(probAgeA(indA));
    for i = 1:5
        BaseMortalityC(i,J) = dot(RelativeMortNoTx(i,indA),PrevalenceF(indA))/sum(PrevalenceF(indA));
        Impact(i,J) =         dot(RelativeMortTx(i,indA),PrevalenceF(indA))/sum(PrevalenceF(indA));
        TxCoverage_A(i,J) = dot(TxCoverage(indA,i)',PrevalenceF(indA))/sum(PrevalenceF(indA));
        StgDistribution(i,J) = dot(StgDist(indA,i)',PrevalenceF(indA))/sum(PrevalenceF(indA));
        DwellTimes(i,J) = dot(DwellTimesF(i, indA),PrevalenceF(indA))/sum(PrevalenceF(indA));
    end
end
MortA(end)=MortA(end-1);
MortC(end)=MortC(end-1);

dlmwrite('RelativeMortalityTxtoNoTx.csv',(Impact + 1) ./ (BaseMortalityC + 1));


%TimeToDiagnosis (exponential distribution)
TimeToDiagnosis = zeros(size(DwellTimes,1),size(DwellTimes,2));
for i = 1 : size(DwellTimes,2)
    
    sumTime = zeros(size(DwellTimes,1));
    for j = 1:nsims
        randNum = rand;
        Insitu = - log(randNum) * DwellTimes(1,i);
        Local = - log(randNum) * DwellTimes(2,i);
        Regional = - log(randNum) * DwellTimes(3,i);
        Distant = - log(randNum) * DwellTimes(4,i);
        StageIV = - log(randNum) * DwellTimes(5,i);
        
        sumTime(1) = sumTime(1) + 0.5*Insitu; % multiplying by 0.5 when we are considering average dwell time and by 1 when consideering full dwell time
        sumTime(2) = sumTime(2) + Insitu + 0.5*Local;
        sumTime(3) = sumTime(3) + Insitu + Local + 0.5*Regional;
        sumTime(4) = sumTime(4) + Insitu + Local + Regional + 0.5*Distant;
        sumTime(5) = sumTime(5) + Insitu + Local + Regional + Distant+ 0.5*StageIV;
    end
    TimeToDiagnosis(1,i) = sumTime(1) / nsims;
    TimeToDiagnosis(2,i) = sumTime(2) / nsims;
    TimeToDiagnosis(3,i) = sumTime(3) / nsims;
    TimeToDiagnosis(4,i) = sumTime(4) / nsims;
    TimeToDiagnosis(5,i) = sumTime(5) / nsims;
end


TimeToDiagnosis


%Modified TimeToDiagnosis (using DCO to estimate event of first occurence )
ModTimeToDiagnosis = zeros(size(DwellTimes,1),size(DwellTimes,2));

CumStageDISTwithDCO = StageDist;
DCO = 0.03;% proportion DCO  
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
   ModTimeToDiagnosis(5,i) = (DwellTimes(1,i) + DwellTimes(2,i) + DwellTimes(3,i) + DwellTimes(4,i)+(1 - CumStageDISTwithDCO (5))/((1/ DwellTimes(5,i))+ MortA(1,i))* CumStageDISTwithDCO(5));
end 

%SojournTime and SojournTimePlusInsitu
SojournTime =  zeros(1,size(DwellTimes,2));
SojournTimePlusInsitu = zeros(1,size(DwellTimes,2));
for j = 1: size(DwellTimes,2)
    SojournTime(1,j) = StageDist(1) * TimeToDiagnosis(1,j) + StageDist(2) * TimeToDiagnosis(2,j) + StageDist(3) * TimeToDiagnosis(3,j) +StageDist(4) * TimeToDiagnosis(4,j)+StageDist(5) * TimeToDiagnosis(5,j);
    SojournTimePlusInsitu(1,j) = SojournTime(1,j);
    
    %     sumTime = 0;
    %     for z = 1: nsims
    %         randNum = rand;
    %         sumTime = sumTime - log(randNum) * SojournTime(1,j) - log(randNum) * DwellTimes(1,j);
    %     end
    %     SojournTimePlusInsitu(1,j) = sumTime / nsims;
    
end%for j

sumprobAgeA = zeros(1, size(AgeArrayLower,2));
for j = 1: size(AgeArrayLower,2)
    for  probAgeInd = AgeArrayLower(j) : AgeArrayUpper(j)
        sumprobAgeA(1,j) = sumprobAgeA(1,j) + probAgeA(probAgeInd);
    end
end

%Calculating incidence by solving the differential equations
Incidence = IncidenceFunc(probAgeA,MortalityArray,prevalenceActual,cancerMortalityArray, AgeArrayLower,AgeArrayUpper);
%Incidence = Incidence *1000;
Incidence =PAR{4};%Cervical cancer has incidence data (breast cancer has prevalence data and IncidenceFunc is used for calculating Incidence)
%Incidence = transpose(smooth(Incidence));

IAll = interp1(AgeArrayMid,Incidence,[0:100],'PCHIP');

% X=15:100;
% JordanInc = [0.000007	0.000007	0.000007	0.000007	0.000007	0.000007	0.000007	0.000007	0.000007	0.000007	0.000007	0.000007	0.000007	0.000007	0.000007	0.000007	0.000007	0.000007	0.000007	0.000007	0.000007	0.000007	0.000007	0.000007	0.000007	0.00005	0.00005	0.00005	0.00005	0.00005	0.000087	0.000087	0.000087	0.000087	0.000087	0.000057	0.000057	0.000057	0.000057	0.000057	0.000099	0.000099	0.000099	0.000099	0.000099	0.000051	0.000051	0.000051	0.000051	0.000051	0.00006	0.00006	0.00006	0.00006	0.00006	0.000027	0.000027	0.000027	0.000027	0.000027	0.000075	0.000075	0.000075	0.000075	0.000075	0.000075	0.000075	0.000075	0.000075	0.000075	0.000075	0.000075	0.000075	0.000075	0.000075	0.000075	0.000075	0.000075	0.000075	0.000075	0.000075	0.000075	0.000075	0.000075	0.000075	0.000075
% ];
% close all
% figure
% hold on
%
% plot(X,IAll(1, 15:100),'b')
% plot(X,JordanInc(1,:),'b*')
% title('Incidence');
% pause
%Calculating onset of insitu by solving the MArkov Chain of the Markov Process
%[OnsetRates] = PrevStageRates (AgeArrayLower,AgeArrayUpper,PAR,TimeToDiagnosis,Impact)
%OnsetRates = PrevStageRatesMod (AgeArrayLower,AgeArrayUpper,PAR,TimeToDiagnosis,Impact);
OnsetRates = PrevStageRates (Incidence, ModTimeToDiagnosis,StageDist,probAgeA,MortalityArray,prevalenceActual,cancerMortalityArray, AgeArrayLower,AgeArrayUpper);
% function [P_theta] = PrevStageRates (AgeArrayLower,AgeArrayUpper,PAR,TimeToDiagnosis,Impact)
%OnsetRates(1,1:4) =0;

% converting to rate
OnsetRates = -log(1-OnsetRates);%
%OnsetMatrix = zeros(1, size(AgeArrayLower,2)+ 1);

%OnsetMatrix(1,1:size(AgeArrayLower,2))= OnsetRates * 1000;
%OnsetMatrix(1,size(AgeArrayLower,2)+1)= (OnsetRates * 1000) * transpose(sumprobAgeA);
%sum((prevalenceActual) * (transpose((SojournTimeMat(i,:)) * transpose(sumprobAgeA)))) / sum(prevalenceActual)
% prevalenceTotal = prevalenceActual(1: size(AgeArrayLower,2))* transpose(sumprobAgeA) *1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Interpolation to single age years, flatline after 80 years

%IncidenceF = interp1(AgeArrayMid,Incidence,[0:100],'PCHIP');
% IncidenceF(80:end)=IncidenceF(80);

%OnsetRatesF = OnsetRates';
OnsetRatesF=interp1(AgeArrayMid,OnsetRates,[0:100],'PCHIP');
% OnsetRatesF(80:end)=OnsetRatesF(80);
OnsetRatesF(1:15) = 0;

Y1=[distHPV(1)* OnsetRatesF(1, 1:101)*1000]';
Y2=[distHPV(2)* OnsetRatesF(1, 1:101)*1000]';

PrevHPV=zeros(12,13);
HPVincidence =zeros(12,30);
HPVprevalence =zeros(12,30);
betaMat=zeros(10,12);
for loopOpt = 0:0
Mark=1; %fit 1/18 + HR if MArk= 1
[InsituOnsetRate,HPVonsetRate,PrevG,HPVonsetG,Immunity,Prate,startPoints] = CalcHPVOnsetFeb2015new(Y1, AgeArrayLower,AgeArrayUpper,PAR,Mark);
ImmunityP(:,1)=PrevG(:,8);
PrevHPV=PrevHPV+PrevG;
HPVonsetH(1:12,1) = HPVonsetG(1:12,1);
HPVonsetH(1:12,2) = HPVonsetG(1:12,2);

Mark=0; %LR  if MArk= 0
[InsituOnsetRate,HPVonsetRate,PrevG,HPVonsetG,Immunity,startPoints] = CalcHPVOnsetFeb2015new(Y2, AgeArrayLower,AgeArrayUpper,PAR,Mark);
ImmunityP(:,2)=PrevG(:,8);
PrevHPV=PrevHPV+PrevG;
PrevHPV_H(:,1)=PrevHPV(:,1)+PrevHPV(:,2)+PrevHPV(:,3)+PrevHPV(:,7);
PrevHPV_H(:,2)=PrevHPV(:,10)+PrevHPV(:,11)+PrevHPV(:,12)+PrevHPV(:,13);
PrevHPV_L=PrevHPV(:,4)+PrevHPV(:,5)+PrevHPV(:,6)+PrevHPV(:,9);


HPVonsetL(1:12,1) = HPVonsetG(1:12,3);
% % 16/18+HR and LR are fit separately. So scale the rates. 
% HPVonsetH(:,1) = HPVonsetH(:,1).*(ones(12,1)-PrevHPV_H(:,1)-PrevHPV_H(:,2)-ImmunityP(:,1))./(ones(12,1)-PrevHPV_H(:,1)-PrevHPV_H(:,2)-PrevHPV_L(:,1)-PrevHPV(:,8));
% HPVonsetH(:,2) = HPVonsetH(:,2).*(ones(12,1)-PrevHPV_H(:,1)-PrevHPV_H(:,2)-ImmunityP(:,1))./(ones(12,1)-PrevHPV_H(:,1)-PrevHPV_H(:,2)-PrevHPV_L(:,1)-PrevHPV(:,8));
% HPVonsetL(:,1) = HPVonsetL(:,1).*(ones(12,1)-PrevHPV_L(:,1)-ImmunityP(:,2))./(ones(12,1)-PrevHPV_H(:,1)-PrevHPV_H(:,2)-PrevHPV_L(:,1)-PrevHPV(:,8));

HPVonsetG(:,1) = HPVonsetH(:,1);
HPVonsetG(:,2) = HPVonsetH(:,2);
HPVonsetG(:,3) = HPVonsetL(:,1);



PrevW(:,1)=PrevHPV(:,1)+PrevHPV(:,2)+PrevHPV(:,3)+PrevHPV(:,7);
PrevW(:,2)=PrevHPV(:,10)+PrevHPV(:,11)+PrevHPV(:,12)+PrevHPV(:,13);
PrevW(:,3)=PrevHPV(:,4)+PrevHPV(:,5)+PrevHPV(:,6)+PrevHPV(:,9);

HPVincidence(:,loopOpt *3+1:loopOpt *3+3)=HPVonsetG;
HPVprevalence(:,loopOpt *3+1:loopOpt *3+3)=PrevW;
%betaMat(loopOpt+1,:)=startPoints(1,:);
end
[OnsetW,RES2] =Transmission(HPVonsetG,PrevW,PrevHPV,PAR);
%HPVonsetG(:,3)=OnsetW(13:end)/1000; %Taking onset rates for low riks from model fit in transmission.
Mark=2;% run both together  
Y3=HPVonsetG;
[InsituOnsetRate,HPVonsetRate,PrevG,HPVonsetG,Immunity] = CalcHPVOnsetFeb2015new(Y3, AgeArrayLower,AgeArrayUpper,PAR,Mark);

HPVonsetG(:,1) = HPVonsetH(:,1);
HPVonsetG(:,2) = HPVonsetH(:,2);
HPVonsetG(:,3) = HPVonsetL(:,1);
HPVonsetG(1:3,:)=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i = 3 : size(OnsetRates,2)
%      indA=AgeArrayLower(i):AgeArrayUpper(i);
%     onsetRateSmooth(1,i)= dot(OnsetRatesF(indA),PrevalenceF(indA))/sum(PrevalenceF(indA)) ;%PrevalenceF(AgeArrayMid(i));%
% end

% OnsetRatesF=InSituOnset;


%beta0=[5,0.001, 5,5];

trial_run = 1;
betaRuns = ones(trial_run,size(AgeArrayLower,2));

for randRuns = 1:trial_run
    beta0 = ones(1,size(AgeArrayLower,2)) * rand()*20;
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
    [Z,PrevS,PrevC,diagrT,beta]=FitSoJourn(1000,beta0,X,Y,PAR,Incidence,OnsetRatesF,DwellTimesF,AgeArrayLower,AgeArrayUpper, Mu_Tx, BaseMortality,RelativeMortNoTx,RelativeMortTx,InsituOnsetRate);
    toc
    betaRuns(randRuns,:)= beta;
    initial_point_array(randRuns,1)=initial_point;
    
end

trial_detail=zeros(trial_run,size(beta,2)+1);
trial_detail(:,1) = initial_point_array;
trial_detail(:,2:(size(beta,2)+1)) = betaRuns;
[DiagnosticRate]=1./ExpSourjX(beta,AgeArrayMid, StageDist);

for i = 2: 5
    %DiagnosticRateMean(i) = (OnsetRates' * DiagnosticRate(:,i))/sum(OnsetRates');
end
%DiagnosticRateMean
%save results to be mapped to associaitons file
RES{1}=PrevS;%preclinical prevalence
RES{2}=PrevC;%clinical prevalence
RES{3}=Incidence;%transition to clincal states
RES{4}=OnsetRates;%onsetRateSmooth';
RES{5}=DiagnosticRate;
RES{6}=DwellTimes';
RES{7}=BaseMortalityC; %Relative Mortality without treatmet
RES{8}=MortA';% Disease free mortality
RES{9} = Impact;% %Relative Mortality with treatment
RES{10} = TxCoverage_A;% Current coverage of treatment
RES{11}=PrevHPV;
RES{12}=HPVonsetG;
RES{13}=Immunity;
end%end of MarkvoChainCountry

function [PrevA,PrevS,PrevC,diagr,beta]=FitSoJourn(MaxFunEvals,beta0,X,Y,PAR,Incidence,OnsetRatesF,DwellTimesF,AgeArrayLower,AgeArrayUpper, Mu_Tx, BaseMortality,RelativeMortNoTx,RelativeMortTx,InsituOnsetRate)
LB = ones(1,size(AgeArrayLower,2))*0.001;
UB = ones(1,size(AgeArrayLower,2));
%LB=[-10,-1, -10,-10];
%UB=[100,0.5, 10, 100];%if exponential curve for non-US
%UB = [2000, 0.5,10,100];%% if logistic US
%UB = [10000, 0.5,10,1000];%% if logistic non-US

%OPTIONS = optimset('largescale','off','LevenbergMarquardt','on','display','off','Diagnostics','on','MaxFunEvals',700,'TolFun',1e-8,'TolCon',1e-8);
OPTIONS = optimset('MaxFunEvals',MaxFunEvals,'Diagnostics','on','TolFun',1e-6,'TolCon',1e-6);

% [beta,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@(x,X)...
%     PopSim(x,X,PAR,Incidence,OnsetRatesF,DwellTimesF,AgeArrayLower,AgeArrayUpper,Mu_Tx, BaseMortality,RelativeMortNoTx,RelativeMortTx, 0),...
%     beta0,X,Y,LB,UB,OPTIONS);

store_beta=beta0;
for age = 1 : size(AgeArrayLower,2)
    individualX = X(1,age);

    [betafit_out,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@(x,individualX )...
    PopSim(x,individualX,PAR,Incidence,OnsetRatesF,DwellTimesF,AgeArrayLower,AgeArrayUpper,Mu_Tx, BaseMortality,RelativeMortNoTx,RelativeMortTx, 0, age,store_beta),...
    beta0(1,age),individualX,Y(1,age),LB(1,age),UB(1,age),OPTIONS);

    store_beta(1,age)=betafit_out;
end
beta=store_beta;

%beta = [0,0,0,8.04350879264037,2.06419265914802,1.05243573362542,1.32322102488666,0.865267058630671,0.373372280464398,0.189700814075861,0.105110504312803,0.0734006218779887];

%plot the results
[PrevA,diagr,PrevS,PrevC]=PopSim(beta,X,PAR,Incidence,OnsetRatesF,DwellTimesF,AgeArrayLower,AgeArrayUpper,Mu_Tx, BaseMortality,RelativeMortNoTx,RelativeMortTx,1,1000000,store_beta);

end

function [diagr]=ExpSourj(b, StageDist, AgeArrayLower, AgeArrayUpper)
diagr=[];
X=0:100;
%b(3)=-1;
diagr1=ones(1,101);
diagr2=ones(1,101);
diagr3=ones(1,101);
diagr4=ones(1,101);
diagr5=ones(1,101);

%%%If US dont take cumulative
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
    value =  (StageDist(5)+StageDist(4)+StageDist(3)+StageDist(2)+StageDist(1));
    if value == 0
        diagr5(AgeArrayLower(i):AgeArrayUpper(i)) = 100000000;
    else
        diagr5(AgeArrayLower(i):AgeArrayUpper(i))=(1/b(1,i))/value;
    end
end


diagr=[diagr1;diagr2;diagr3;diagr4;diagr5]';
end

function [diagr]=ExpSourjX(b,X, StageDist);
diagr=[];


%%%If US dont take cumulative
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
    value = (StageDist(5)+StageDist(4)+StageDist(3)+StageDist(2)+StageDist(1));
    if value == 0
        diagr5(i) = 100000000;
    else
        diagr5(i)=(1/b(1,i))/value;
    end
end



diagr=[diagr1;diagr2;diagr3;diagr4;diagr5]';
%diagr=max(diagr,0.5);

end

function [Z,diagr,PrevS,PrevC]=PopSim(fitPAR,X,PAR,Incidence,OnsetRatesF,DwellTimesF,AgeArrayLower,AgeArrayUpper, Mu_Tx, BaseMortality,RelativeMortNoTx,RelativeMortTx, PlotModel, age, store_beta)
%

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
if StageDist(5) == 0
    diagr(:,5) = 0;
    diagrT(:,5) = 0;
end
% diagr(:,1) =smooth(diagr(:,1) );
% diagr(:,2) =smooth(diagr(:,2) );
% diagr(:,3) =smooth(diagr(:,3) );
% diagr(:,4) =smooth(diagr(:,4) );
% diagr(:,5) =smooth(diagr(:,5) );
% diagr(2:4) = StageDist(1:3)./diagrT(2:4);

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

POP=zeros(m_a,11);
POP(1:m_a,1)=probAgeA*pop0;

t_max=500;
CFR = zeros(m_a,t_max);

%%
dt=1/15;
for t=1:t_max
    
    %fast time loop
    Diags=zeros(length(POP),5,t_max);
    IncidenceE = zeros(101,5);
     death_pre_clinical_5 = zeros(101,1);
    for t1=1:1/dt
        %healhy pop
        POP(:,1)=POP(:,1) + dt*(-mu-OnsetRatesF).*POP(:,1);
        
        %non-diagnosed cancer
        %POP(:,2)=POP(:,2)+dt*(-mu-lambda(:,1)).*POP(:,2) + dt*(OnsetRatesF).*POP(:,1);
        POP(:,2)=POP(:,2)+dt*(-mu-lambda(:,1)-diagr(:,1)).*POP(:,2) + dt*(OnsetRatesF).*POP(:,1);
        POP(:,3)=POP(:,3)+dt*(-mu-lambda(:,2)-diagr(:,2)).*POP(:,3) + dt*(lambda(:,1)).*POP(:,2);
        POP(:,4)=POP(:,4)+dt*(-mu-lambda(:,3)-diagr(:,3)).*POP(:,4) + dt*(lambda(:,2)).*POP(:,3);
        POP(:,5)=POP(:,5)+dt*(-mu-lambda(:,4)-diagr(:,4)).*POP(:,5) + dt*(lambda(:,3)).*POP(:,4);
        POP(:,6)=POP(:,6)+dt*(-mu-lambda(:,5)-diagr(:,5)).*POP(:,6) + dt*(lambda(:,4)).*POP(:,5);
        %POP(:,6)=POP(:,6)+dt*(-mu.* (1+ RelativeMortNoTx(1,:)')-diagr(:,5)).*POP(:,6) + dt*(lambda(:,4)).*POP(:,5);
        
        %% 
          % Deaths in preclinical last stage
             death_pre_clinical_5(:,1) = death_pre_clinical_5(:,1) + POP(:,6).*(lambda(:,5))*dt;
        %diagnosed cancer
        IncidenceE(:,1) = IncidenceE(:,1) + dt*(diagr(:,1)).*POP(:,2);
        %POP(:,7)=POP(:,7)+dt*(-mu).*POP(:,7) - dt*BaseMortality(1,:)'.*POP(:,7)+ dt*(diagr(:,1)).*POP(:,2);
        POP(:,7)=POP(:,7)+dt*(-mu .* (1+ RelativeMortTx(1,:)') ).*POP(:,7) + dt*(diagr(:,1)).*POP(:,2);
        CFR(:,t)= CFR(:,t)+ dt*(mu ).*POP(:,7)+ dt*BaseMortality(1,:)'.*POP(:,7);
        
        IncidenceE(:,1) = IncidenceE(:,1) + dt*(diagr(:,2)).*POP(:,3);
        %POP(:,8)=POP(:,8)+dt*(-mu ).*POP(:,8) - dt*BaseMortality(2,:)'.*POP(:,8)+ dt*(diagr(:,2)).*POP(:,3);
        POP(:,8)=POP(:,8)+dt*(-mu  .* (1+ RelativeMortTx(2,:)') ).*POP(:,8) + dt*(diagr(:,2)).*POP(:,3);
        CFR(:,t)= CFR(:,t)+ dt*(mu ).*POP(:,8)+ dt*BaseMortality(2,:)'.*POP(:,8);
        
        IncidenceE(:,1) = IncidenceE(:,1) + dt*(diagr(:,3)).*POP(:,4);
        %POP(:,9)=POP(:,9)+dt*(-mu ).*POP(:,9) - dt*BaseMortality(3,:)'.*POP(:,9)+ dt*(diagr(:,3)).*POP(:,4);
        POP(:,9)=POP(:,9)+dt*(-mu.* (1+ RelativeMortTx(3,:)') ).*POP(:,9) + dt*(diagr(:,3)).*POP(:,4);
        CFR(:,t)= CFR(:,t)+dt*(mu ).*POP(:,9) + dt*BaseMortality(3,:)'.*POP(:,9);
        
        IncidenceE(:,1) = IncidenceE(:,1) + dt*(diagr(:,4)).*POP(:,5);
        %POP(:,10)=POP(:,10)+dt*(-mu ).*POP(:,10) - dt*BaseMortality(4,:)'.*POP(:,10)+ dt*(diagr(:,4)).*POP(:,5);
        POP(:,10)=POP(:,10)+dt*(-mu .* (1+ RelativeMortTx(4,:)')).*POP(:,10) + dt*(diagr(:,4)).*POP(:,5);
        CFR(:,t)= CFR(:,t)+dt*(mu ).*POP(:,10) + dt*BaseMortality(4,:)'.*POP(:,10);
        
        IncidenceE(:,1) = IncidenceE(:,1) + dt*(diagr(:,5)).*POP(:,6);
        %POP(:,11)=POP(:,11)+dt*(-mu ).*POP(:,11) - dt*BaseMortality(5,:)'.*POP(:,11)+ dt*(diagr(:,5)).*POP(:,6);
        POP(:,11)=POP(:,11)+dt*(-mu .* (1+ RelativeMortTx(5,:)')).*POP(:,11) + dt*(diagr(:,5)).*POP(:,6);
        CFR(:,t)= CFR(:,t)+dt*(mu ).*POP(:,11) + dt*BaseMortality(5,:)'.*POP(:,11);
        
        % POP(:,7)=POP(:,7) - dt*BaseMortality(2,:)'.*POP(:,7)+ dt*(diagr(:,2)).*POP(:,3);
        % POP(:,8)=POP(:,8) - dt*BaseMortality(3,:)'.*POP(:,8)+ dt*(diagr(:,3)).*POP(:,4);
        % POP(:,9)=POP(:,9) - dt*BaseMortality(4,:)'.*POP(:,9)+ dt*(diagr(:,4)).*POP(:,5);
        
        
        
        %if t > t_max - 100
        Diags(:,1,t)=Diags(:,1,t)+dt*( ((diagr(:,1)).*POP(:,2)) );
        Diags(:,2,t)=Diags(:,2,t)+dt*( ((diagr(:,2)).*POP(:,3)) );
        Diags(:,3,t)=Diags(:,3,t)+dt*( ((diagr(:,3)).*POP(:,4)) );
        Diags(:,4,t)=Diags(:,4,t)+dt*( ((diagr(:,4)).*POP(:,5)) );
        Diags(:,5,t)=Diags(:,5,t)+dt*( ((diagr(:,5)).*POP(:,6)) );
        %end
        
    end %%%%%%
    
    DiagsAll=sum(Diags,2);
    
    %aging
    for I=1:11
        POP(2:end,I)=POP(1:end-1,I);
    end
    
    %must decide if births are to be based on population size
    Births=BirthsF*100000;
    POP(1,1)=Births;
    
    
    
    DiagsA(1,t) = sum(Diags(:,1,t))/ sum(sum(Diags(:,:,t)));
    DiagsA(2,t) = sum(Diags(:,2,t))/ sum(sum(Diags(:,:,t)));
    DiagsA(3,t) = sum(Diags(:,3,t))/ sum(sum(Diags(:,:,t)));
    DiagsA(4,t) = sum(Diags(:,4,t))/ sum(sum(Diags(:,:,t)));
    DiagsA(5,t) = sum(Diags(:,5,t))/ sum(sum(Diags(:,:,t)));
    
    
end
%prevalence of cancer by age
%average mort in each age group
prevAI=zeros(1,101);
for J=1:length(AgeArrayLower)
    indA=AgeArrayLower(J):AgeArrayUpper(J);
    
    inc=IncidenceE(indA,1);
    
    p_a=POP(indA,:);
    p_c=POP(indA,[7:11]);
    
    %pre clinical
    p_s1=POP(indA,[2]);
    p_s2=POP(indA,[3]);
    p_s3=POP(indA,[4]);
    p_s4=POP(indA,[5]);
    p_s5=POP(indA,[6]);
    
    %diagnosed/clinical
    p_c1=POP(indA,[7]);
    p_c2=POP(indA,[8]);
    p_c3=POP(indA,[9]);
    p_c4=POP(indA,[10]);
    p_c5=POP(indA,[11]);
    
    %%%%%%
      total_inc = inc';
        total_inc = sum(total_inc,1);
    prevAI(1,indA)=(inc(:)'./sum(p_a(:,:)'));
    %      for k=indA(1,1):indA(end)
    %          prevAI(1,k)=sum(inc(k-indA(1,1)+1,:))/sum(p_a(k-indA(1,1)+1,:))*100;
    %      end
    %%%%%%
    PrevA(J)=sum(p_c(:))/sum(p_a(:));
    incidenceCalc(J) = sum(inc(:))/sum(p_a(:));
    
    %prev cancer pre clinical
    PrevS(J,1)=sum(p_s1(:))/sum(p_a(:));
    PrevS(J,2)=sum(p_s2(:))/sum(p_a(:));
    PrevS(J,3)=sum(p_s3(:))/sum(p_a(:));
    PrevS(J,4)=sum(p_s4(:))/sum(p_a(:));
    PrevS(J,5)=sum(p_s5(:))/sum(p_a(:));
    
    %prev cancer clinical
    PrevC(J,1)=sum(p_c1(:))/sum(p_a(:));
    PrevC(J,2)=sum(p_c2(:))/sum(p_a(:));
    PrevC(J,3)=sum(p_c3(:))/sum(p_a(:));
    PrevC(J,4)=sum(p_c4(:))/sum(p_a(:));
    PrevC(J,5)=sum(p_c5(:))/sum(p_a(:));
    
end

%Diags=Diags/sum(Diags)*100;
%Z=[PrevA*100];
% if age >1000
%     Z=incidenceCalc;
% else
%     Z=incidenceCalc(age)*100;
% end
if age >1000
        Z=incidenceCalc;      
         DCO =sum(death_pre_clinical_5) / (sum(death_pre_clinical_5)+sum(sum(IncidenceE)));
    else
        % incidenceCalc (simulated incidence) already considers state of
        % cancer, see line ~713.
        Z=incidenceCalc(1,age)*100; 
         DCO =sum(death_pre_clinical_5) / (sum(death_pre_clinical_5)+sum(sum(IncidenceE)));
    end
    
    % plotting populations
    for i = 6:9
        plot(POP(:,i));
        hold on
        title('POP(:,2:9)');
    end
if(PlotModel==1)
    
    prevInd=xlsread('prevInd.xlsx','prevAI');
    prevInd(PAR{13},1:101)=prevAI*100;
    xlswrite('prevInd.xlsx',prevInd,'prevAI');
    
    close all
    figure
    subplot(4,1,1)
    plot(POP(:,1))
    title('total pop')
    
    subplot(4,1,2)
    plot(sum(POP(:,[2:6]),2) )
    title('cancer not diagnosed')
    
    subplot(4,1,3)
    plot(sum(POP(:,[7:11]),2) )
    title('total pop')
    title('cancer diagnosed')
    
    subplot(4,1,4)
    plot(OnsetRatesF )
    title('total pop')
    title('Onset Rate')
    
    
    PrevA=Z(1:length(prev))%cervical cancer has incidence data not prevalence. Validating against incidence
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
    
    
    
    DiagrT_weighted(:,1)=diagrT(:,1)'.*probAgeA;
    DiagrT_weighted(:,2)=diagrT(:,2)'.*probAgeA;
    DiagrT_weighted(:,3)=diagrT(:,3)'.*probAgeA;
    DiagrT_weighted(:,4)=diagrT(:,4)'.*probAgeA;
    DiagrT_weighted(:,5)=diagrT(:,5)'.*probAgeA;
    
    figure
    hold on
    X=0:100;
    plot(X,diagrT(:,1),'b')
    plot(X,diagrT(:,2),'r')
    plot(X,diagrT(:,3),'g')
    plot(X,diagrT(:,4),'c')
    plot(X,diagrT(:,5),'y')
    legend('in-situ/0', 'local/I', 'regional/II', 'distant/III', 'IV');
    title('Time to diagnosis')
    
    figure
    hold on
    X=0:100;
    plot(X,DwellTimesF(1,:),'b')
    plot(X,DwellTimesF(2,:),'r')
    plot(X,DwellTimesF(3,:),'g')
    plot(X,DwellTimesF(4,:),'c')
    plot(X,DwellTimesF(5,:),'y')
    legend('in-situ/0', 'local/I', 'regional/II', 'distant/III', 'IV');
    
    title('Dwell Time')
    
    figure
    hold on
    plot(1:5,DiagsA(1:5,t),'r-')
    plot(1:5,StageDist(1:5),'k*')
    title('Stage at diagnosis')
    
    
    figure
    hold on
    subplot(5,1,1)
    plot(1:500,DiagsA(1,:),'r-')
    title('In-situ')
    
    subplot(5,1,2)
    plot(1:500,DiagsA(2,:),'r-')
    title('Local/I')
    
    subplot(5,1,3)
    plot(1:500,DiagsA(3,:),'r-')
    title('Regional/II')
    
    subplot(5,1,4)
    plot(1:500,DiagsA(4,:),'r-')
    title('Distant/III')
    
    subplot(5,1,5)
    plot(1:500,DiagsA(5,:),'r-')
    title('IV')
    
    
    [~,~,TEMPMARKOV1]=xlsread('Markov Results.xlsx','Total Pop');
    [~,~,TEMPMARKOV2]=xlsread('Markov Results.xlsx','prev & prevA');
    [~,~,TEMPMARKOV3]=xlsread('Markov Results.xlsx','Onset Rate');
    [~,~,TEMPMARKOV4]=xlsread('Markov Results.xlsx','CFR');
    [~,~,TEMPMARKOV5]=xlsread('Markov Results.xlsx','Time to diagnosis');
    [~,~,TEMPMARKOV6]=xlsread('Markov Results.xlsx','Dwell Time');
    [~,~,TEMPMARKOV7]=xlsread('Markov Results.xlsx','Stage at diagnosis');
    [~,~,TEMPMARKOV8]=xlsread('Markov Results.xlsx','DiagsA');
    
    %pause
    CountryName = PAR{14};
    CountryIndex=PAR{13};
    [AA,BB,CC,DD,EE,FF,GG,HH] = MikeWriteTemplateMarkov(TEMPMARKOV1,TEMPMARKOV2,TEMPMARKOV3,TEMPMARKOV4,TEMPMARKOV5,TEMPMARKOV6,TEMPMARKOV7,TEMPMARKOV8,PAR{13},POP,OnsetRatesF,CFR,diagr,DwellTimesF,DiagsA,StageDist,AgeArrayMid,prev,PrevA,CountryName)
    
    
    
    xlswrite('Markov Results.xlsx',AA,'Total Pop')
    xlswrite('Markov Results.xlsx',BB,'prev & prevA')
    xlswrite('Markov Results.xlsx',CC,'Onset Rate')
    xlswrite('Markov Results.xlsx',DD,'CFR')
    xlswrite('Markov Results.xlsx',EE,'Time to diagnosis')
    xlswrite('Markov Results.xlsx',FF,'Dwell Time')
    xlswrite('Markov Results.xlsx',GG,'Stage at diagnosis')
    %   xlswrite('Markov Results.xlsx',HH,'DiagsA')
    
    
    'Markov Results Written'
    
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

function [Rates] = PrevStageRates (CurrentStageRate, TimeToDiagnosis, StageDist, probAgeA, MortalityArray, prevalenceActual, cancerMortalityArray, AgeArrayLower, AgeArrayUpper)

% %NOTE: cervical cancer has one  stage more than breast cancer
% 
% 
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
%    A_a=C_a;
%                                     % here, we have made an assumption that
%                                   % A_a = C_a, because probability of a
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
%     lamda=zeros(5,size(AgeArrayLower,2));
%     for i=(1:5)
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
%     
%     %%
%     startAge = 2;
%     for a = startAge+1:12   % START ALGORITHM 
%       
%     %%STEP 1
%         IC = I_D_a(1,a) * C_a(1,a);
%         
%         sumKloop =0;
%         for k = startAge:a-1
%                       
%            sumIloopNumerator = 0;
%            t = (AgeArrayUpper(1,a)-AgeArrayLower(1,k))/2;
%             for i = 1:5
%                 sumIloopNumerator = sumIloopNumerator + s_i(1,i)* (  (1- exp(-t*lamda(i,k))) - (1- exp(-(t-1)*lamda(i,k))) );
%             end
%             
%             
%             
%             muProduct=1;
%             for j = k:a
%              muProduct = muProduct * exp(-(AgeArrayUpper(j)-AgeArrayLower(j-1))/2 *mu(1,j));
%             end
%             sumKloop = sumKloop + steady_h(1,k) * probability_hu(1,k) * sumIloopNumerator * muProduct;
%         end
%         
%         
%         sumIloopDenominator = 0;
%         t = (AgeArrayUpper(1,a)-AgeArrayLower(1,a))/2;
%         for i = 1:5
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
%         for k = startAge:a-1
%                       
%            sumIloopNumerator = 0;
%           % t = mid_age_group(1,a)-mid_age_group(1,k);
%            t = (AgeArrayUpper(1,a)-AgeArrayLower(1,k))/2;
%             for i = 1:5
%                 sumIloopNumerator = sumIloopNumerator + s_i(1,i)*  exp(-t*lamda(i,k)) ;
%             end
%             muProduct=1;
%             for j = k:a
%              muProduct = muProduct * exp(-(AgeArrayUpper(j)-AgeArrayLower(j-1))/2 *mu(1,j));
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
    

numStage = size(StageDist,2);
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
           t = mid_age_group(1,a)-mid_age_group(1,k);
%             t = AgeArrayUpper(1,a)-mid_age_group(1,k);
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
    Rates=probability_hu;
    
    
    
    onset_ref=[0	0	0	0.000329083	0.000768594	0.00132553	0.001327089	0.002014344	0.001415079	0.001059426	0.001065192	0.000921512];


%     onset_ref=[1.77095E-16
% 5.01067E-11
% 1.05561E-08
% 7.7181E-07
% 0.000147044
% 0.001099142
% 0.002539701
% 0.004250965
% 0.003305252
% 0.003121092
% 0.003330814
% 0.003539146
% ];

% onset_ref = [0
% 0
% 0
% 0
% 0.02234
% 0.00963
% 0.00357
% 0.00027
% 0
% 0
% 0
% 0];

    x=(1:101);
    figure
    plot(mid_age_group(1:12),onset_ref,'--',mid_age_group,OnsetRates);
    legend('Reference values of Onset -SEARB','Estimated');
    title('Cervical Cancer, SEARB');
    
    y=1;
    
    
    
end %PrevStageRates


% %ESTIMATE Previous Stage
% function [Rates] = PrevStageRates (CurrentStageRate, TimeToDiagnosis, StageDist, probAgeA, MortalityArray, prevalenceActual, cancerMortalityArray, AgeArrayLower, AgeArrayUpper)
% 
% %NOTE: cervical cancer has one  stage more than breast cancer
% 
% 
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
%     lamda=zeros(5,size(AgeArrayLower,2));
%     for i=(1:5)
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
%             for i = 1:5
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
%         for i = 1:5
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
%             for i = 1:5
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
%     
%     OnsetRates = -log(1- probability_hu);
%     Rates=probability_hu;
%     
%     onset_ref=[0	0	0	0.000329083	0.000768594	0.00132553	0.001327089	0.002014344	0.001415079	0.001059426	0.001065192	0.000921512];
% 
% 
% %     onset_ref=[1.77095E-16
% % 5.01067E-11
% % 1.05561E-08
% % 7.7181E-07
% % 0.000147044
% % 0.001099142
% % 0.002539701
% % 0.004250965
% % 0.003305252
% % 0.003121092
% % 0.003330814
% % 0.003539146
% % ];
% 
% % onset_ref = [0
% % 0
% % 0
% % 0
% % 0.02234
% % 0.00963
% % 0.00357
% % 0.00027
% % 0
% % 0
% % 0
% % 0];
% 
%     x=(1:101);
%     figure
%     plot(mid_age_group(1:12),onset_ref,'--',mid_age_group,OnsetRates);
%     legend('Reference values of Onset','Estimated');
%     title('Cervical Cancer, SEARB');
%     
%     y=1;
%     
%     
%     
% end %PrevStageRates
% 

% apply stage at diagnosis to incidence to get persons by stage
% Use dwell times to estimate onset rate in Insitu


