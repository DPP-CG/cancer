
function [RES,RES2]=MarkovChainCountry2_StageAndAge(PAR)
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
AgeArrayM=AgeArrayL+(AgeArrayU-AgeArrayL)/2;

CFR = interp1(AgeArrayM,PAR{5},[0:100],'PCHIP');

cancerMortalityArray =-log(1-min(.99,CFR));
%cancerMortalityArray(1,80:end) = cancerMortalityArray(1,79);
%cancerMortalityArray = max(cancerMortalityArray, MortalityArray);
PAR{5} = cancerMortalityArray;
StageDist=PAR{6};

%Spectrum Age groups
AgeArrayLower=[1    5	10	15  20	25	30	40	50	60	70	80];
AgeArrayUpper=[4	9	14  19 	24	29	39	49	59	69	79	100];
AgeArrayMid=(AgeArrayLower+(AgeArrayUpper-AgeArrayLower)/2);

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
%DwellTimes(1,:)= 1./[0.03 0.03 0.273 0.273 0.273];%[0.01 0.04 .04 .04];

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
        
        sumTime(1) = sumTime(1) + 0.5*Insitu;   % multiplying by 0.5 when we are considering average dwell time and by 1 when consideering full dwell time
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
OnsetRates = PrevStageRates (Incidence, ModTimeToDiagnosis,StageDist /(1+DCO),probAgeA,MortalityArray,prevalenceActual,cancerMortalityArray, AgeArrayLower,AgeArrayUpper);
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

%% heavy loop

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
% 16/18+HR and LR are fit separately. So scale the rates. 
%     HPVonsetH(:,1) = HPVonsetH(:,1).*(ones(12,1)-PrevHPV_H(:,1)-PrevHPV_H(:,2)-ImmunityP(:,1))./(ones(12,1)-PrevHPV_H(:,1)-PrevHPV_H(:,2)-PrevHPV_L(:,1)-ImmunityP(:,1)- ImmunityP(:,2));
%     HPVonsetH(:,2) = HPVonsetH(:,2).*(ones(12,1)-PrevHPV_H(:,1)-PrevHPV_H(:,2)-ImmunityP(:,1))./(ones(12,1)-PrevHPV_H(:,1)-PrevHPV_H(:,2)-PrevHPV_L(:,1)-ImmunityP(:,1)- ImmunityP(:,2));
%     HPVonsetL(:,1) = HPVonsetL(:,1).*(ones(12,1)-PrevHPV_L(:,1)-ImmunityP(:,2))./(ones(12,1)-PrevHPV_H(:,1)-PrevHPV_H(:,2)-PrevHPV_L(:,1)-ImmunityP(:,1)- ImmunityP(:,2));

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

%% test by plotting
for i = 5:12
    hold on
    plot(HPVincidence(i,1:3:30));
end

%%

[OnsetW,RES2] =Transmission(HPVonsetG,PrevW,PrevHPV,PAR);
%HPVonsetG(:,3)=OnsetW(13:end)/1000; %Taking onset rates for low riks from model fit in transmission.
Mark=2;% run both together  
Y3=HPVonsetG;
[InsituOnsetRate,HPVonsetRate,PrevG,~,Immunity] = CalcHPVOnsetFeb2015new(Y3, AgeArrayLower,AgeArrayUpper,PAR,Mark);
HPVonsetG(1:3,:)=0;
%% skip loop

%% SEARB

% with intervention
% HPVonsetG = [0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0];
% HPVonsetH = [1.59872115546023e-14,6.21724893790088e-15;1.02117876643140e-13,3.97125075834435e-14;3.17544545634893e-14,1.23489545524681e-14;7.20090441947218e-06,2.80035171868363e-06;-0.0753582465883098,-0.0293059847843427;0.237948147751969,0.0925353907924325;0.178497857014148,0.0694158332832799;0.173941191843273,0.0676437968279395;0.141863914344123,0.0551693000227146;0.135882388628376,0.0528431511332572;0.0666246227164099,0.0259095755008261;0.164537544565826,0.0639868228867101];
% HPVonsetL = [2.22044604925031e-14;4.04796927246621e-14;9.03982359652431e-06;1.00008934531708e-05;-0.0679033753598050;-0.00104120821096358;-0.000285267659593809;-0.000345302733341982;-0.000181113202069411;-0.000250838108067644;-0.000240451599630176;-0.000333872291149686];
% Immunity = 10;
% res2_1 = [0,0,0,0;0,0,0,0;0,0,0,0;7.92706998378146e-17,3.08274943813723e-17,2.49559797876111e-17,2.20001125524335e-18;1.62280787516432e-05,6.31091951452789e-06,4.44621563667336e-06,4.39613139441838e-07;0.0208603202147878,0.00811234675019525,0.00545453470285062,0.000718968645816150;0.499328605966228,0.194183346764644,0.129715763500358,0.112626028763194;0.460743157794411,0.179177894697826,0.119644290837892,0.224674788713586;0.446293805543067,0.173558702155637,0.111962070465458,0.252286010429123;0.448013917119232,0.174227634435257,0.102220829652627,0.259698129348848;0.451909113030970,0.175742432845377,0.0888674839804717,0.260847618745033;0.461629904766852,0.179522740742665,0.0783090827496817,0.260972705636401];
% res2_2 = [0,2.22044604925031e-14;0,2.22044604925031e-14;0,2.22044604925031e-14;0.0244415938684172,2.22044604925031e-14;0.219470524017777,2.22044604925031e-14;0.755657143583482,0.00166721709385975;0.455851147827351,0.889812051545996;0.309659063474317,0.277482266935187;0.347412841001692,0.437907835822268;0.275303490574089,0.437850457851390;0.152607129551308,0.283850265275873;0.350661649554295,0.672713889597017];
% res2_3 = [0,0,0;0,0,0;0,0,0;0,0,0;1.01821237461871e-14,3.95971479018388e-15,3.20553338023210e-15;0.00189392999012623,0.000736528329493534,0.000518857873625894;2.08901801758512,0.812395895727549,0.542283014657587;0.679628220944397,0.264299863700599,0.177574304768651;1.06380480755994,0.413701869606643,0.215196188375810;1.07998086915820,0.419992560228190,0.151473490360207;0.685124966600415,0.266437487011272,0.0755783689817683;1.62274174472747,0.631066234060684,0.148913061677433];
% RES2 = [{res2_1},{res2_2},{res2_3}];

% w/o intervention
% HPVonsetG = [0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0];
% HPVonsetH = [1.82076576038526e-14,3.99680288865056e-15;2.66560171128178e-14,5.85132082964294e-15;8.19996980335702e-06,1.79999337146861e-06;8.20097691620650e-06,1.80021444502094e-06;-0.0720632046270772,-0.0158187522352121;0.0154118175713411,0.00338308190590415;0.00891003623478899,0.00195586161251466;0.00374116917843320,0.000821232258680459;0.00258842038342906,0.000568189840264915;0.00474502308168639,0.00104159043256531;0.00664279709938126,0.00145817497303491;0.0240204940182189,0.00527279136985294];
% HPVonsetL = [2.22044604925031e-14;3.69685292874745e-09;8.35029849130506e-06;1.00013308125076e-05;-0.0385971108783342;-0.000941326147060541;-0.000574881047198138;-0.000369399243256977;-0.000362889070933058;-0.000479244278847936;-0.000456397588814123;-0.000628424290365078];
% Immunity = 10;

%% AFRE (yet to be changed, copied from PERU for now)

% with intervention
% HPVonsetG = [0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0];
% HPVonsetH = [1.82076576038526e-14,3.99680288865056e-15;2.66560171128178e-14,5.85132082964294e-15;8.19996980335702e-06,1.79999337146861e-06;8.20097691620650e-06,1.80021444502094e-06;-0.0720632046270772,-0.0158187522352121;0.0154118175713411,0.00338308190590415;0.00891003623478899,0.00195586161251466;0.00374116917843320,0.000821232258680459;0.00258842038342906,0.000568189840264915;0.00474502308168639,0.00104159043256531;0.00664279709938126,0.00145817497303491;0.0240204940182189,0.00527279136985294];
% HPVonsetL = [2.22044604925031e-14;3.69685292874745e-09;8.35029849130506e-06;1.00013308125076e-05;-0.0385971108783342;-0.000941326147060541;-0.000574881047198138;-0.000369399243256977;-0.000362889070933058;-0.000479244278847936;-0.000456397588814123;-0.000628424290365078];
% Immunity = 10;
% res2_1 = [0,0,0,0;0,0,0,0;0,0,0,0;0.00351927010188690,0.00527890515283036,0.0330368264043774,0.000711404582864694;0.0720939012560607,0.108140851884091,0.664470584224971,0.0695347369582200;0.0716577486279779,0.107486622941967,0.648850965450262,0.165460834824832;0.0637356225482437,0.0956034338223656,0.557616866637596,0.222686503843265;0.0623031611470127,0.0934547417205190,0.472679491085724,0.230408650546094;0.0611052762354016,0.0916579143531025,0.410266255701008,0.224227815678473;0.0457647197354643,0.0686470796031965,0.298979812376601,0.190679256878754;0.0370240256836279,0.0555360385254418,0.228349195426778,0.152663058892401;0.0297404394023137,0.0446106591034705,0.174653977189639,0.121758320071171];
% res2_2 = [0,2.22044604925031e-14;0,2.22044604925031e-14;0,2.22044604925031e-14;0.701751047937015,2.22044604925031e-14;0.666526568492246,0.999973313390550;0.292772622997245,0.615159732231950;0.0287395461896347,3.31705897593952e-14;0.0326187562590745,0.0603227557815062;0.00247760028328365,4.30337491592369e-14;0.00352142375635848,0.00607351326697128;0.00154519048137304,0.00286144490070230;0.000823625761366586,0.00133690435051796];
% res2_3 = [0,0,0;0,0,0;0,0,0;0,0,0;0.219049986610874,0.328574979916311,1.83481151617194;0.265534669813290,0.398302004719935,1.83485302938324;2.19734951815118e-14,3.29602427722676e-14,1.10910358204504e-13;0.0301090582780300,0.0451635874170449,0.137998736107661;1.49051798955957e-14,2.23577698433936e-14,6.58328364626225e-14;0.00174550642732549,0.00261825964098824,0.00694555839053685;0.000843463158162018,0.00126519473724303,0.00297194503463120;0.000430081169810328,0.000645121754715491,0.00125461155167200];
% RES2 = [{res2_1},{res2_2},{res2_3}];

% w/o intervention
% HPVonsetG = [0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0];
% HPVonsetH = [1.82076576038526e-14,3.99680288865056e-15;2.66560171128178e-14,5.85132082964294e-15;8.19996980335702e-06,1.79999337146861e-06;8.20097691620650e-06,1.80021444502094e-06;-0.0720632046270772,-0.0158187522352121;0.0154118175713411,0.00338308190590415;0.00891003623478899,0.00195586161251466;0.00374116917843320,0.000821232258680459;0.00258842038342906,0.000568189840264915;0.00474502308168639,0.00104159043256531;0.00664279709938126,0.00145817497303491;0.0240204940182189,0.00527279136985294];
% HPVonsetL = [2.22044604925031e-14;3.69685292874745e-09;8.35029849130506e-06;1.00013308125076e-05;-0.0385971108783342;-0.000941326147060541;-0.000574881047198138;-0.000369399243256977;-0.000362889070933058;-0.000479244278847936;-0.000456397588814123;-0.000628424290365078];
% Immunity = 10;


%% PERU

% % with intervention
% HPVonsetG = [0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0];
% Immunity = 10;
% HPVonsetH = [8.88178419700125e-15,1.33226762955019e-14;4.33445460612014e-11,6.50168190918020e-11;1.48643515610105e-14,2.22965273415157e-14;1.57079361177360e-14,2.35619041766040e-14;-0.0201049487003780,-0.0301574230505670;0.0673322143343671,0.100998321501551;0.0577688555353866,0.0866532833030799;0.0669418904417652,0.100412835662648;0.0672256865177116,0.100838529776567;0.0689118204744028,0.103367730711604;0.0684430405533930,0.102664560830089;0.0726918988955634,0.109037848343345];
% HPVonsetL = [2.22044604925031e-14;4.44089013187371e-14;3.95821353701741e-14;2.22044609940273e-14;-0.0101359380137796;0.200732897481630;0.136563839454177;0.154767150963948;0.158753414365271;0.162160936197415;0.164076732504672;0.165845308868065];
% res2_1 = [0,0,0,0;0,0,0,0;0,0,0,0;4.60788044729090e-17,6.91182067093636e-17,1.80141694134700e-16,4.81096321857456e-18;6.31413651337382e-15,9.47120477006073e-15,2.36115227548585e-14,1.95022087314643e-15;0.00862523848460485,0.0129378577269073,0.0223845776390381,0.000749257986181553;0.165655818511349,0.248483727767024,0.429888594654271,0.118223223340933;0.150547313754471,0.225820970631706,0.387016218055812,0.227097159451951;0.147815494883572,0.221723242325359,0.369303452843013,0.254696875728964;0.147998693530265,0.221998040295398,0.360917113447403,0.262056830478880;0.149404553049905,0.224106829574857,0.357627256871374,0.264183083776963;0.150165414835639,0.225248122253458,0.354976661246763,0.265238762355024];
% res2_2 = [0,2.22044604925031e-14;0,2.22044604925031e-14;0,2.22044604925031e-14;0.0220554012712575,2.22044604925031e-14;0.207273810400778,2.22044604925031e-14;0.714249152881696,2.24477451585106e-14;0.437603992610239,0.839983368281933;0.376728449653175,0.278731950508790;0.479527079702899,0.675600632016029;0.443569307099780,0.508104180560315;0.471404760765494,0.891331084755137;0.475465244321430,0.772670871116584];
% res2_3 = [0,0,0;0,0,0;0,0,0;0,0,0;3.86357599035588e-15,5.79536398553382e-15,1.41891279716702e-14;1.13938139985910e-14,1.70907209978864e-14,3.53595882733800e-14;1.03358258383098,1.55037387574647,2.68239650826922;0.384651942738712,0.576977914108069,0.895854421878507;0.928341002928582,1.39251150439287,2.14195717405588;0.719223992330565,1.07883598849585,1.65051649138512;1.28586933432691,1.92880400149037,2.94401938001384;1.13072370791895,1.69608556187842,2.58416179982556];
% res2_4 = [0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0.357234012000000,0.522377791000000,0.119867991000000,0.000520165000000000,4.13429000000000e-08,6.01843000000000e-14,1.60468000000000e-21,3.86399000000000e-36;0,0,0,0,0.103761012000000,0.412439816000000,0.468759712000000,0.0150306330000000,8.82727000000000e-06,9.49505000000000e-11,1.87064000000000e-17,9.04736000000000e-31;0,0,0,0,0.00184802700000000,0.0363836160000000,0.534920730000000,0.420783550000000,0.00606247800000000,1.59979000000000e-06,7.73213000000000e-12,4.54406000000000e-23;0,0,0,0,2.43644000000000e-07,3.54440000000000e-05,0.0127840410000000,0.549054445000000,0.431901524000000,0.00622266100000000,1.64206000000000e-06,3.89315000000000e-15;0,0,0,0,5.73209000000000e-13,6.16153000000000e-10,5.45201000000000e-06,0.0127844490000000,0.549071944000000,0.431915289000000,0.00622285900000000,5.95208000000000e-09;0,0,0,0,2.48495000000000e-20,1.97370000000000e-16,4.28442000000000e-11,5.48523000000000e-06,0.0128623460000000,0.552417496000000,0.434546993000000,0.000167680000000000;0,0,0,0,3.01117000000000e-29,1.76721000000000e-24,9.41114000000000e-18,6.57844000000000e-11,8.42221000000000e-06,0.0197492840000000,0.848200665000000,0.132041629000000;0,0,0,0,4.34976000000000e-46,5.12746000000000e-40,3.31794000000000e-31,9.35658000000000e-22,4.83267000000000e-14,4.57171000000000e-08,0.000792124000000000,0.999207830000000];
% res2_5 = [0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0.538135176000000,0.368009879000000,0.0925833930000000,0.00127138900000000,1.63305000000000e-07,3.84188000000000e-13,1.65542000000000e-20,1.30646000000000e-29,1.60203000000000e-46;0,0,0,0.237150835000000,0.440846338000000,0.301477798000000,0.0205055660000000,1.94618000000000e-05,3.38310000000000e-10,1.07713000000000e-16,6.28125000000000e-25,1.54704000000000e-40;0,0,0,0.0143939950000000,0.132530273000000,0.448904914000000,0.394970506000000,0.00919639100000000,3.92186000000000e-06,3.06330000000000e-11,4.38236000000000e-18,1.31153000000000e-31;0,0,0,1.15825000000000e-05,0.000787995000000000,0.0197220380000000,0.425700768000000,0.541171740000000,0.0126005020000000,5.37357000000000e-06,4.19720000000000e-11,5.06753000000000e-22;0,0,0,1.26404000000000e-10,6.35437000000000e-08,1.17514000000000e-05,0.00622278600000000,0.431910189000000,0.549065461000000,0.0127842980000000,5.45195000000000e-06,2.65556000000000e-14;0,0,0,2.49035000000000e-17,9.25042000000000e-14,1.26406000000000e-10,1.64212000000000e-06,0.00622288300000000,0.431916927000000,0.549074026000000,0.0127844970000000,2.51220000000000e-08;0,0,0,9.09856000000000e-26,2.49725000000000e-21,2.52149000000000e-17,8.03595000000000e-12,1.66265000000000e-06,0.00630069200000000,0.437317466000000,0.555939457000000,0.000440722000000000;0,0,0,1.69731000000000e-41,9.35694000000000e-36,1.89763000000000e-30,7.34862000000000e-23,6.13390000000000e-15,9.37757000000000e-09,0.000262583000000000,0.134667680000000,0.865069728000000];
% RES2 = [{res2_1},{res2_2},{res2_3},{res2_4},{res2_5}];

% % w/o screen;
% HPVonsetG = [0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0];
% Immunity = 10;
% HPVonsetH = [8.88178419700125e-15,1.33226762955019e-14;2.01873974525104e-12,3.02810961787657e-12;1.36589329962913e-14,2.04883994944369e-14;1.05736359276578e-09,1.58604538914866e-09;-0.0358018753127187,-0.0537028129690781;0.0674666599456066,0.101199989918410;0.0579706257303483,0.0869559385955224;0.0668608700736604,0.100291305110491;0.0685234775782775,0.102785216367416;0.0702719131957931,0.105407869793690;0.0707763033607811,0.106164455041172;0.0727168954862932,0.109075343229440];
% HPVonsetL = [2.22044604925031e-14;1.47993346797773e-11;2.29439351132292e-14;8.09110035436526e-13;-0.0826555313750642;0.202275983044186;0.137036261114152;0.154593907118893;0.157491719927059;0.160928051469537;0.162106034904298;0.165741764559457];
% res2_1 = [0,0,0,0;0,0,0,0;0,0,0,0;0.00351927010188690,0.00527890515283036,0.0330368264043774,0.000711404582864694;0.0720939012560607,0.108140851884091,0.664470584224971,0.0695347369582200;0.0716577486279779,0.107486622941967,0.648850965450262,0.165460834824832;0.0637356225482437,0.0956034338223656,0.557616866637596,0.222686503843265;0.0623031611470127,0.0934547417205190,0.472679491085724,0.230408650546094;0.0611052762354016,0.0916579143531025,0.410266255701008,0.224227815678473;0.0457647197354643,0.0686470796031965,0.298979812376601,0.190679256878754;0.0370240256836279,0.0555360385254418,0.228349195426778,0.152663058892401;0.0297404394023137,0.0446106591034705,0.174653977189639,0.121758320071171];
% res2_2 = [0,2.22044604925031e-14;0,2.22044604925031e-14;0,2.22044604925031e-14;0.701751047937015,2.22044604925031e-14;0.666526568492246,0.999973313390550;0.292772622997245,0.615159732231950;0.0287395461896347,3.31705897593952e-14;0.0326187562590745,0.0603227557815062;0.00247760028328365,4.30337491592369e-14;0.00352142375635848,0.00607351326697128;0.00154519048137304,0.00286144490070230;0.000823625761366586,0.00133690435051796];
% res2_3 = [0,0,0;0,0,0;0,0,0;0,0,0;0.219049986610874,0.328574979916311,1.83481151617194;0.265534669813290,0.398302004719935,1.83485302938324;2.19734951815118e-14,3.29602427722676e-14,1.10910358204504e-13;0.0301090582780300,0.0451635874170449,0.137998736107661;1.49051798955957e-14,2.23577698433936e-14,6.58328364626225e-14;0.00174550642732549,0.00261825964098824,0.00694555839053685;0.000843463158162018,0.00126519473724303,0.00297194503463120;0.000430081169810328,0.000645121754715491,0.00125461155167200];
% RES2 = [{res2_1},{res2_2},{res2_3}];

%%


% HPVonsetG(:,1) = HPVonsetH(:,1);
% HPVonsetG(:,2) = HPVonsetH(:,2);
% HPVonsetG(:,3) = HPVonsetL(:,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i = 1 : size(OnsetRates,2)
%      indA=AgeArrayLower(i):AgeArrayUpper(i);
%     onsetRateSmooth(1,i)= dot(OnsetRatesF(indA),PrevalenceF(indA))/sum(PrevalenceF(indA)) ;%PrevalenceF(AgeArrayMid(i));%
% end

% OnsetRatesF=InSituOnset;

%% calculation of diagnostic rates starts here

    for upper_start = 1:1
        trial_run = 1;
        betaRuns = ones(5,size(AgeArrayLower,2),trial_run);

        for randRuns = 1:trial_run
            
            beta0 = ones(5,size(AgeArrayLower,2)) * rand() * (upper_start);
            initial_point = beta0(:,1);

            % ExpSourj function was needed when we were fitting the diagnostic
            % rates values only for distant state of cancer. Then by using the
            % fitted values, ExpSourj fuction gave output of whole diagnostic
            % rate values. As we are now fitting each age and stage, we don't
            % need this function.
            %[diagr] = ExpSourj(beta0, StageDist,AgeArrayLower, AgeArrayUpper);


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
            %Y=StageDist(1,1:4)'*[prevalenceActual*100];%if fit to prevalence
            Y=StageDist(1,1:5)'*(Incidence*100);%if fit to incidence
            %[Z,diagr]=PopSim(beta0,X,PAR,IncidenceF,OnsetRatesF,DwellTimesF,AgeArrayLower,AgeArrayUpper,1)
            toc

            %Fit the prevalence model, find the diagnostic rates and plot the results
            tic
            [Z,PrevS,PrevC,diagrT,beta,err]=FitSoJourn(1000,beta0,X,Y,PAR,Incidence,OnsetRatesF,DwellTimesF,AgeArrayLower,AgeArrayUpper, Mu_Tx, BaseMortality,RelativeMortNoTx,RelativeMortTx,0,randRuns);
            toc
            
            %% storing the results
            % diagnostic rates
            betaRuns(:,:,randRuns) = beta;
            % initial points
            initial_point_array(:,randRuns)=initial_point;
            % error
            run_error(:,randRuns) = err;
        end
        
        %%
        trial_detail=zeros(5,size(beta,2)+2,trial_run);
        trial_detail(:,1,:) = initial_point_array;
        %trial_detail(:,2,:) = run_error;
        trial_detail(:,3:end,:) = betaRuns;
        for i = 1:trial_run
            trial_details(:,:,i) = transpose(trial_detail(:,:,i));
        end

        %% selecting diagnostic rate matrix with the least error
        [~,idx] = min(run_error);
        DiagnosticRate = betaRuns(:,:,idx)';
        
        % plotting error
        figure
        plot(initial_point_array(1,:),run_error,'r+');
        xlabel('Staring point for algorithm');
        ylabel('Error');
        title(' Error for each trial ');

        %% collecting output in tractable format
        for i_l = 1:trial_run
            for j = 2:5
                stack((j-2)*14 + 1:(j-2)*14 + 14,i_l) = trial_details(:,j,i_l);
            end
        end
       
       % storing the results
       RES{13 + upper_start} = stack;
        
    end

% end of diagnostic rates calculations

%% Collecting results for sending it to 'process models' file

    %save results to be mapped to associaitons file
    RES{1}=PrevS;
    RES{2}=PrevC;
    RES{3}=Incidence;
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

function [PrevA,PrevS,PrevC,diagr,beta,err]=FitSoJourn(MaxFunEvals,beta0,X,Y,PAR,Incidence,OnsetRatesF,DwellTimesF,AgeArrayLower,AgeArrayUpper, Mu_Tx, BaseMortality,RelativeMortNoTx,RelativeMortTx,InsituOnsetRate,randRuns)
    
    % lower and the upper bounds
    LB = ones(5,size(AgeArrayLower,2))*0.01;
    UB = ones(5,size(AgeArrayLower,2))*4;
    
    % following values are made equal to zero because
    LB(:,1:3)   =0; % first three age groups are considered to have no cases
    UB(:,1:3)   =0;
    LB(1,:)     =0; % Pre-clinical 1 doesn't have any cases
    UB(1,:)     =0;
    
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
    store_beta = beta0;
    for age = 1 : size(AgeArrayLower,2)
        individualX = X(1,age);
        for state_c = 1:5
            [betafit_out,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@(x,individualX )...
            PopSim(x,individualX,PAR,Incidence,OnsetRatesF,DwellTimesF,AgeArrayLower,AgeArrayUpper,Mu_Tx, BaseMortality,RelativeMortNoTx,RelativeMortTx,0, age,store_beta,Y,state_c),...
            beta0(1,age),individualX,Y(state_c,age),LB(state_c,age),UB(state_c,age),OPTIONS);
        
            % storing values of diagnostic rates in different stages for a
            % specific age
            store_beta(state_c,age) = betafit_out;
        end
    end
    beta=store_beta;

    % beta = [0,0,0,8.04350879264037,2.06419265914802,1.05243573362542,1.32322102488666,0.865267058630671,0.373372280464398,0.189700814075861,0.105110504312803,0.0734006218779887];

    %plot the results
    [PrevA,diagr,PrevS,PrevC,err]=PopSim(beta,X,PAR,Incidence,OnsetRatesF,DwellTimesF,AgeArrayLower,AgeArrayUpper,Mu_Tx, BaseMortality, RelativeMortNoTx, RelativeMortTx,1,10^5,store_beta,Y,state_c);

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

function [diagr]=ExpSourjX(b,X, StageDist)
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

function [Z,diagr,PrevS,PrevC,err]=PopSim(fitPAR,X,PAR,Incidence,OnsetRatesF,DwellTimesF,AgeArrayLower,AgeArrayUpper, Mu_Tx, BaseMortality,RelativeMortNoTx,RelativeMortTx, PlotModel, age, store_beta,Y,state_c)
%

    Z=[];
    PrevS=[];
    PrevC=[];

    AgeArrayMid=(AgeArrayUpper+AgeArrayLower)/2;

    StageDist=PAR{6};

    %% calculating diagnostic rates for each age ( 1 to 101 )

    % PREVIOUS CALCULATIONS : In this case, we are fitting values to incidence 
    % in each age and NOT fitting to incidence in each age and each state
    %     if age < 13
    %         if length(fitPAR)<2
    %             current_beta = fitPAR;
    %             fitPAR = store_beta;
    %             fitPAR(1,age) = current_beta;
    %         end
    %     end  
    %     [diagrT]=ExpSourj(fitPAR,StageDist, AgeArrayLower, AgeArrayUpper);
    %     diagr = 1./diagrT

    % NEW CALCULATIONS: In this case we are fitting values to incidence in each
    % age and stage of cancer.
    diagr = zeros(5,101);
    for ex_p = 1:length(AgeArrayLower)
        for c_s = 1:5
            diagr(c_s,AgeArrayLower(1,ex_p):AgeArrayUpper(1,ex_p)) = store_beta(c_s,ex_p);
        end
    end
    if age ~= 10^5
        if age <= 6
            if mod(age,2) ~= 0
                diagr(state_c,(5*(age-1) + 1):(5*(age-1) + 1) + 3) = fitPAR;
            elseif mod(age,2) == 0
                diagr(state_c,(5*(age-1)):(5*(age-1)) + 5) = fitPAR;
            end
        elseif age <= 11
            diagr(state_c,(5*(2*age-8)):(5*(2*age-8)) + 9) = fitPAR;
        else
            diagr(state_c,80:101) = fitPAR;
        end
    elseif age == 10^5
        diagrT = fitPAR;
        for i=1:12
            if i <= 6
                if mod(i,2) ~= 0
                    for j=1:5
                        diagr(j,(5*(i-1) + 1):(5*(i-1) + 1) + 3) = diagrT(j,i);
                    end
                elseif mod(i,2) == 0
                    for j=1:5
                        diagr(j,(5*(i-1)):(5*(i-1)) + 5) = diagrT(j,i);
                    end
                end
            elseif i <= 11
                for j=1:5
                    diagr(j,(5*(2*i-8)):(5*(2*i-8)) + 9) = diagrT(j,i);
                end
            else
                for j=1:5
                    diagr(j,80:101) = diagrT(j,i);
                end
            end
        end
    end

    diagr = diagr';

    % diagr(:,1) =smooth(diagr(:,1) );
    % diagr(:,2) =smooth(diagr(:,2) );
    % diagr(:,3) =smooth(diagr(:,3) );
    % diagr(:,4) =smooth(diagr(:,4) );
    % diagr(:,5) =smooth(diagr(:,5) );
    % diagr(2:4) = StageDist(1:3)./diagrT(2:4);
    
    %%
    
    %diagr=1./diagrT;
    if StageDist(1) == 0
        diagr(:,1) = 0;
        diagrT(:,1) = 0;
    end
    if StageDist(5) == 0
        diagr(:,5) = 0;
        diagrT(:,5) = 0;
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

            IncidenceE(:,2) = IncidenceE(:,2) + dt*(diagr(:,2)).*POP(:,3);
            %POP(:,8)=POP(:,8)+dt*(-mu ).*POP(:,8) - dt*BaseMortality(2,:)'.*POP(:,8)+ dt*(diagr(:,2)).*POP(:,3);
            POP(:,8)=POP(:,8)+dt*(-mu  .* (1+ RelativeMortTx(2,:)') ).*POP(:,8) + dt*(diagr(:,2)).*POP(:,3);
            CFR(:,t)= CFR(:,t)+ dt*(mu ).*POP(:,8)+ dt*BaseMortality(2,:)'.*POP(:,8);

            IncidenceE(:,3) = IncidenceE(:,3) + dt*(diagr(:,3)).*POP(:,4);
            %POP(:,9)=POP(:,9)+dt*(-mu ).*POP(:,9) - dt*BaseMortality(3,:)'.*POP(:,9)+ dt*(diagr(:,3)).*POP(:,4);
            POP(:,9)=POP(:,9)+dt*(-mu.* (1+ RelativeMortTx(3,:)') ).*POP(:,9) + dt*(diagr(:,3)).*POP(:,4);
            CFR(:,t)= CFR(:,t)+dt*(mu ).*POP(:,9) + dt*BaseMortality(3,:)'.*POP(:,9);

            IncidenceE(:,4) = IncidenceE(:,4) + dt*(diagr(:,4)).*POP(:,5);
            %POP(:,10)=POP(:,10)+dt*(-mu ).*POP(:,10) - dt*BaseMortality(4,:)'.*POP(:,10)+ dt*(diagr(:,4)).*POP(:,5);
            POP(:,10)=POP(:,10)+dt*(-mu .* (1+ RelativeMortTx(4,:)')).*POP(:,10) + dt*(diagr(:,4)).*POP(:,5);
            CFR(:,t)= CFR(:,t)+dt*(mu ).*POP(:,10) + dt*BaseMortality(4,:)'.*POP(:,10);

            IncidenceE(:,5) = IncidenceE(:,5) + dt*(diagr(:,5)).*POP(:,6);
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

        inc=IncidenceE(indA,:);

        p_a=POP(indA,:);
        p_c=POP(indA,7:11);

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
        prevAI(1,indA)=(total_inc./sum(p_a(:,:)'));
        %      for k=indA(1,1):indA(end)
        %          prevAI(1,k)=sum(inc(k-indA(1,1)+1,:))/sum(p_a(k-indA(1,1)+1,:))*100;
        %      end
        %%%%%%
        PrevA(J)=sum(p_c(:))/sum(p_a(:));
        
        
        if age ~= 10^5
            incidenceCalc(J) = sum(inc(:,state_c))/sum(p_a(:));
        elseif age == 10^5
            for state_c = 1:5
                incidenceCalc(state_c,J) = sum(inc(:,state_c))/sum(p_a(:));
            end
        end

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

%         prevInd=xlsread('prevInd.xlsx','prevAI');
%         prevInd(PAR{13},1:101)=prevAI*100;
%         xlswrite('prevInd.xlsx',prevInd,'prevAI');

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


%         PrevA=Z(1:length(prev));%cervical cancer has incidence data not prevalence. Validating against incidence
%         figure
%         hold on
%         plot(AgeArrayMid,prev*100,'k*')
%         plot(AgeArrayMid,PrevA,'r-')
%         title('Prevalence of Diagnosed cases by age')


        %% incidence fit plot 
        if age ~= 10^5
            
            % calculating actual incidence for each state of cancer
            Incidence_state = StageDist'*Incidence;
            
            figure
            hold on
            plot(AgeArrayMid,Incidence_state(state_c,:)*100,'k*')
            plot(AgeArrayMid,incidenceCalc*100,'r-')
            title('Incidence of Diagnosed cases by age')
            
        elseif age == 10^5
            
            % calculating actual incidence for each state of cancer
            Incidence_state = StageDist'*Incidence;
            
            for state_c = 1:5
                figure
                hold on
                plot(AgeArrayMid,Incidence_state(state_c,:)*100,'k*')
                plot(AgeArrayMid,incidenceCalc(state_c,:)*100,'r-')
                formatSpec = 'Incidence of Diagnosed cases by age for cancer-state %d';
                str = sprintf(formatSpec,state_c);
                title(str)
            end
        end
        
        % calculating error
        err1 = sum(sum(abs(Incidence_state - incidenceCalc)*10^5));

        %%
%         figure
%         hold on
%         plot(0:100,OnsetRatesF )
% 
%         title('Onset Rate')
% 
%         figure
%         hold on
%         plot(0:100,CFR(:,t))
%         plot(0:100,CFR(:,t-100))
%         plot(0:100,CFR(:,t-200))
%         plot(0:100,CFR(:,t-300))
%         title('CFR')



%         DiagrT_weighted(:,1)=diagrT(:,1)'.*probAgeA;
%         DiagrT_weighted(:,2)=diagrT(:,2)'.*probAgeA;
%         DiagrT_weighted(:,3)=diagrT(:,3)'.*probAgeA;
%         DiagrT_weighted(:,4)=diagrT(:,4)'.*probAgeA;
%         DiagrT_weighted(:,5)=diagrT(:,5)'.*probAgeA;

%         figure
%         hold on
%         X=0:100;
%         plot(X,diagrT(:,1),'b')
%         plot(X,diagrT(:,2),'r')
%         plot(X,diagrT(:,3),'g')
%         plot(X,diagrT(:,4),'c')
%         plot(X,diagrT(:,5),'y')
%         legend('in-situ/0', 'local/I', 'regional/II', 'distant/III', 'IV');
%         title('Time to diagnosis')
% 
%         figure
%         hold on
%         X=0:100;
%         plot(X,DwellTimesF(1,:),'b')
%         plot(X,DwellTimesF(2,:),'r')
%         plot(X,DwellTimesF(3,:),'g')
%         plot(X,DwellTimesF(4,:),'c')
%         plot(X,DwellTimesF(5,:),'y')
%         legend('in-situ/0', 'local/I', 'regional/II', 'distant/III', 'IV');
% 
%         title('Dwell Time')

        %% plotting and calculating stage diagnosis fit
        figure
        hold on
        plot(1:5,DiagsA(1:5,t)','r-')
        plot(1:5,StageDist(1:5),'k*')
        title('Stage at diagnosis')
        
        % calculating error
        err2 = sum(abs(DiagsA(1:4,t)'-StageDist(1:4))*10^1);
        
        %%

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
        
        %% Returning the error values
        err = sum([err1 err2]);
%% Not sure whether following plot should be included or not

        
%         [~,~,TEMPMARKOV1]=xlsread('Markov Results.xlsx','Total Pop');
%         [~,~,TEMPMARKOV2]=xlsread('Markov Results.xlsx','prev & prevA');
%         [~,~,TEMPMARKOV3]=xlsread('Markov Results.xlsx','Onset Rate');
%         [~,~,TEMPMARKOV4]=xlsread('Markov Results.xlsx','CFR');
%         [~,~,TEMPMARKOV5]=xlsread('Markov Results.xlsx','Time to diagnosis');
%         [~,~,TEMPMARKOV6]=xlsread('Markov Results.xlsx','Dwell Time');
%         [~,~,TEMPMARKOV7]=xlsread('Markov Results.xlsx','Stage at diagnosis');
%         [~,~,TEMPMARKOV8]=xlsread('Markov Results.xlsx','DiagsA');


        %CountryName = PAR{14};
        %CountryIndex=PAR{13};
        %[AA,BB,CC,DD,EE,FF,GG,HH] = MikeWriteTemplateMarkov(TEMPMARKOV1,TEMPMARKOV2,TEMPMARKOV3,TEMPMARKOV4,TEMPMARKOV5,TEMPMARKOV6,TEMPMARKOV7,TEMPMARKOV8,PAR{13},POP,OnsetRatesF,CFR,diagr,DwellTimesF,DiagsA,StageDist,AgeArrayMid,prev,PrevA,CountryName)


% 
%         xlswrite('Markov Results.xlsx',AA,'Total Pop')
%         xlswrite('Markov Results.xlsx',BB,'prev & prevA')
%         xlswrite('Markov Results.xlsx',CC,'Onset Rate')
%         xlswrite('Markov Results.xlsx',DD,'CFR')
%         xlswrite('Markov Results.xlsx',EE,'Time to diagnosis')
%         xlswrite('Markov Results.xlsx',FF,'Dwell Time')
%         xlswrite('Markov Results.xlsx',GG,'Stage at diagnosis')
%         %   xlswrite('Markov Results.xlsx',HH,'DiagsA')

%            'Markov Results Written'

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

%ESTIMATE Previous Stage
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
        probability_hu(1,a) = max(0,(IC - sumKloop)/ denominator);
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
 
    OnsetRates = max(0,-log(1- probability_hu));
    Rates=max(0,probability_hu);
    
    
    
   % onset_ref=[0	0	0	0.000329083	0.000768594	0.00132553	0.001327089	0.002014344	0.001415079	0.001059426	0.001065192	0.000921512];
 onset_ref=[0,0,0,0.000319774967341117,0.000176135100179266,0.000472206331897201,0.000562072912696757,0.000823544048344283,0.000118330351268421,5.20318713942696e-05,0.000122336825117098,0.00106178919423848]

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


% apply stage at diagnosis to incidence to get persons by stage
% Use dwell times to estimate onset rate in Insitu


