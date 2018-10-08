function ProcessModelsMike

[~,~,TEMPLATE]=xlsread('BreastCancerTemplate.xls','BreastCancerTemplate');
[~,~,INTERVENTIONS]=xlsread('BreastCancerTemplate.xls','Interventions');
[~,~,DWeights]=xlsread('BreastCancerTemplate.xls','DW');

%[~,~,TEMPMORT]=xlsread('C:\Users\Mike\Documents\Umass\Spring 2015\Cancer Research\excel\Template for OutputData Markov.xls');
%[~,~,TEMPMARKOV]=xlsread('C:\Users\Mike\Documents\Umass\Spring 2015\Cancer Research\excel\Template for OutputData Mort.xls');


%% Intervention name
for k = 1:4
    INTVN{k,1}=INTERVENTIONS{2*(k+2),2};
    INTVN{k,2}=INTERVENTIONS{2*(k+2),5};
    DW{k,1}=DWeights{k,2};
    DW{k,2}=DWeights{k,3};
end

%pal care
DW{5,2}=DWeights{7,3};
DW{6,2}=DWeights{8,3};

%pal care
DW{5,1}=DW{4,1};
DW{6,1}=DW{4,1};

%% reading data

% population data
[~,~,RAW_D]=xlsread('res4.xlsx');
% cancer data
[~,~,FORMAT]=xlsread('DataInci.xls');

%% defining required matrices

P = (cell2mat(FORMAT(4:13,2:17)))';
S =(cell2mat(FORMAT(32:35,2:17)))';
M = (cell2mat(FORMAT(18:25,2:17)))';

prevalenceActualMatrix=P / 1000;
StageDistMatrixT = S;
StageDistMatrix = StageDistMatrixT;
cancerMortalityArrayMatrix = M / 1000;
AgeArrayL=[1 15 40 45 50  55	60	65 70 75];
AgeArrayU=[14 39 44  49 54	59	64	69 74 100];
for i=1:17
    Countries{i}=cell2mat(FORMAT(3,i+1));
end

%% commented stuff

% %DEBUG TOOL
% prevalenceActualMatrix;
% %DEBUG TOOL
% cancerMortalityArrayMatrix;

% Agedist=prevalenceActualMatrix(1,:)/sum(prevalenceActualMatrix(1,:));
% prevalenceActualMatrix(1,:)=0.000297*Agedist;
%DEBUG TOOL
% StageDistMatrix=zeros(15,4);
% StageDistMatrix(:,1)=0;
% StageDistMatrix(:,2)=StageDistMatrixT(:,1);
% StageDistMatrix(:,3)=StageDistMatrixT(:,2)+StageDistMatrixT(:,3);
% StageDistMatrix(:,4)=StageDistMatrixT(:,4);

%% commented stuff
%DwellTimes
% DwellTimes =[
% [0.83	0.83	0.83	1.01	1.17	1.47	2.00	2.71	3.68	4.98	6.76	9.76];
% [1.46	1.46	1.46	1.78	2.07	2.60	3.53	4.78	6.49	8.79	11.94	17.24];
% [3.14	3.14	3.14	3.83	4.44	5.59	7.59	10.29	13.97	18.91	25.67	37.08];
% [0.14	0.14	0.14	0.17	0.19	0.25	0.33	0.45	0.61	0.83	1.13	1.63]];

%%%  % %%Assuming Local = Stage I+II
%%%  % StageDistMatrix = [%[.864 .086 0.050] %US
%%%  %    [.236 .58 .184];  [.48 .42 .10];[.23 .50 .27]]; %USA, Jordan, Mexico, Ghana ;
%%%
%%%  %%Assuming Regional  = Stage II+III  %%%%% *****  %%%%%
%%%  StageDistMatrix = [%[.21 .48 .26 .04 0] %[0 .49 .46 0.050] %US
%%%     [0 .1 .72 .18 0];  [0 .11 .79 .10 0];[0 .03 .70 .27 0]]; %Jordan, Mexico, Ghana ;
%%%
% % %  prevalenceActualMatrix = [%[ 0.0000010 	 0.0000610 	 0.0008500 	 0.0047240 	 0.0118950 	 0.0217650 	 0.0296630 	 0.0344300 ]%US	from SEER works better, this is at a different age group than the rest
% % % 
% % %                        %[9.87588E-07	4.47848E-07	0.0004747	0.01166894	0.024693709	0.026197707	0.034809053	0.01533131]%US	
% % %                        [2.50589E-07	6.96893E-07	0.000161713	0.004127239	0.006161374	0.004302963	0.002898198	0.001819891	]%Jordan
% % %                        [7.11707E-06	1.79226E-06	0.000151617	0.004328815	0.009580667	0.008402266	0.007294113	0.004731743]% Mexico
% % %                        
% % %                        ]%;

%[ 1 20 30 40 50 60 70 80]%US %US 2011	from SEER works better, this is at a different age group than the rest
   % [1    5   15  30	45	60	70	80]%US 2004 from BRCdata.xls 
    %Jordan
%     [1    5   15  30	45	60	70	80]% Mexico
%     [1    5   15  30	45	60	70	80]% Ghana
%     ];

    %[19 29 39 49 59 69 79 100]
    
%     [14    39   44  49	54	59	64	69 74 100]
%     [4    14  29 	44	59	69	79	100]
%     [4    14  29 	44	59	69	79	100]];

% prevalenceActualMatrix (:,2:6)= 0;    %%%%   ******  %%%%%
%  for   i = 1: size( prevalenceActualMatrix,1)
%  prevalenceActualMatrix (i,:) = smooth(prevalenceActualMatrix(i,:));
%  end

% transpose(prevalenceActualMatrix)

%These are probabilities, converting into rates in MarkovChain
% % % cancerMortalityArrayMatrix = [%[0	0.112401522	0.008219873	0.009297556	0.01782396	0.031465983	0.035359524	0.135330108]%US
% % %                               [ 0	0	0.037194168	0.042929048	0.105885908	0.208111406	0.366232009	1]%Jordan
% % %                             ];

%% extracting country specific data

% for I=2 for AFRE, 
%      =12 for SEARB, 
%      =16 for Peru
%      =17 for Caribbean
I = 17;
    Countries{I}
    ind=find(strcmp(RAW_D(:,1),Countries{I}));
    DATA=RAW_D(ind,6);
    pA=cell2mat(DATA(163:243))';
    Mx=cell2mat(DATA(406:486))';
    Births=cell2mat(DATA(487:487));
    xa=101;
    pA(end+1:101)=0;
    Mx(end+1:101)=0;
    
%%  defining parameters for the Model-1 and Model-2 calculations  

    PAR{1}=pA/sum(pA); %distribtion of each age
    PAR{2}=Mx; %mortality by age
    PAR{3}=Births/sum(pA); %percentage of births compare to total population
    PAR{4}=prevalenceActualMatrix(I,:); %prevalence of the country by age (14 years range)
    PAR{5}=cancerMortalityArrayMatrix(I,:); %mortality by age (14 years range)
    PAR{6}=StageDistMatrix(I,:); %distribution of stages (insitu, local, regional, distant)
    PAR{7}=AgeArrayL(1,:); %lower range of ages
    PAR{8}=AgeArrayU(1,:); %upper range of ages
    PAR{9}=I; %country number
    PAR{10}=Countries(I); %country name

%% Parameters are now fed to the function which calculates the onset rate and diagnostic rates
    [RES]=MarkovChainCountry2_edit(PAR);

%% Writing the output
    [C_TEMPLATE]=WriteTemplate(RES,TEMPLATE);
    formatSpec = '%s-ProcessModOut-BreastCancer-%s.xlsx';
    str = sprintf(formatSpec,datetime('today','Format','yyyy-MM-dd'),Countries{I});
    xlswrite(str,C_TEMPLATE)
        
%% producing file for spectrum input

    % reading the association file
    [~,~,data1] = xlsread('BreastCancer_Asso');

    % reading the process model output
    [~,~,data2] = xlsread(str);

    % writing file for spectrum input
    data2(310:667,:) = data1(310:end,:);
    xlswrite(str,data2);
    'written'
    
%% Function for writing the results in required format
    
function [lTEMPLATE]=WriteTemplate(RES,TEMPLATE)

    lTEMPLATE=TEMPLATE;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Data from Markov and Cohort model

    PrevS=RES{1};
    PrevC=RES{2};
    Incidence=RES{3};
    OnsetRates=RES{4};
    OnsetRates(end)=OnsetRates(end-1);
    %DiagnosticRate=1./RES{5};
    
    %%
    DiagnosticRate=RES{5};

    %%
    % some values in diagnosticRate matrix becomes very large and then
    % gives error in spectrum simulations, following adjustment is to
    % avoid those problems
    for row = 1:3
        for col = 2:4
            if DiagnosticRate(row,col) >= 99
                DiagnosticRate(row,col) = 0;
            end
        end
    end
    %DiagnosticRate(:,1) = 0;

    DwellRate=1./RES{6};
    for J=1:12
        for K=1:4
            if  DwellRate(J,K)>10000
                DwellRate(J,K)=0;
            end
        end
    end

    MortRate= RES{7};%Relative Mortality without treatment
    MortRateDF= RES{8};% Disease free mortality
    Impact = RES{9};%Relative Mortality with treatment
    TxCoverage = RES{10};% Current coverage of treatment

    PrevS(1:3,:)= 0;
    PrevC(1:3,:) =0;
    OnsetRates(1:3)= 0;
    %Impact(:,1:6)=0;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Matrix with index to cell position in TEMPLATE;

    onset=1;
    %Pre-clinical records
    PC1_Prev0=2;PC1_CFR=3;PC1_Screen=4;PC1_Progr=5;
    PC2_Prev0=6;PC2_CFR=7;PC2_Screen=8;PC2_Progr=9;
    PC3_Prev0=10;PC3_CFR=11;PC3_Screen=12;PC3_Progr=13;
    PC4_Prev0=14;PC4_CFR=15;PC4_Screen=16;

    %Clinical records
    C1_Prev0=17;C1_CFR=18;C1_Progr=19;
    C2_Prev0=20;C2_CFR=21;C2_Progr=22;
    C3_Prev0=23;C3_CFR=24;C3_Progr=25;
    C4_Prev0=26;C4_CFR=27;

    %Treatment records
    I1_CFR=28;
    I2_CFR=29;
    I3_CFR=30;
    I4_CFR=31;

    %Palliative + Treatment records
    I5_CFR_Function=32;
    I6_CFR_Function=33;

    %Palliative records
    I7_Function=34;
    I8_Function=35;


    %Screening
    I1_Screen=36;
    I2_Screen=37;
    I3_Screen=38;

    %Diability weight records
    D_C1=39;
    D_C2=40;
    D_C3=41;
    D_C4=42;

    % Treatment coverage records
    Tx_C1 = 43;
    Tx_C2 = 44;
    Tx_C3 = 45;
    Tx_C4 = 46;

    Tx_C5 = 47;
    Tx_C6 = 48;
    Tx_C7 = 49;
    Tx_C8 = 50;
    %1=row 1, 2=row end, 3=row col

    %disease free
    R{onset,1}=9;R{onset,2}=20;R{onset,3}=7;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %pre clinical
    %stage 0
    R{PC1_Prev0,1}=141;R{PC1_Prev0,2}=152;R{PC1_Prev0,3}=4;
    R{PC1_CFR,1}=141;R{PC1_CFR,2}=152;R{PC1_CFR,3}=6;
    R{PC1_Screen,1}=141;R{PC1_Screen,2}=152;R{PC1_Screen,3}=9;
    R{PC1_Progr,1}=141;R{PC1_Progr,2}=152;R{PC1_Progr,3}=12;

    %stage 2
    R{PC2_Prev0,1}=163;R{PC2_Prev0,2}=174;R{PC2_Prev0,3}=4;
    R{PC2_CFR,1}=163;R{PC2_CFR,2}=174;R{PC2_CFR,3}=6;
    R{PC2_Screen,1}=163;R{PC2_Screen,2}=174;R{PC2_Screen,3}=9;
    R{PC2_Progr,1}=163;R{PC2_Progr,2}=174;R{PC2_Progr,3}=12;

    %stage 3
    R{PC3_Prev0,1}=185;R{PC3_Prev0,2}=196;R{PC3_Prev0,3}=4;
    R{PC3_CFR,1}=185;R{PC3_CFR,2}=196;R{PC3_CFR,3}=6;
    R{PC3_Screen,1}=185;R{PC3_Screen,2}=196;R{PC3_Screen,3}=9;
    R{PC3_Progr,1}=185;R{PC3_Progr,2}=196;R{PC3_Progr,3}=12;

    %stage 4
    R{PC4_Prev0,1}=207;R{PC4_Prev0,2}=218;R{PC4_Prev0,3}=4;
    R{PC4_CFR,1}=207;R{PC4_CFR,2}=218;R{PC4_CFR,3}=6;
    R{PC4_Screen,1}=207;R{PC4_Screen,2}=218;R{PC4_Screen,3}=9;
    %R{PC4_Progr,1}=207;R{PC4_Progr,2}=218;R{PC4_Progr,3}=12;  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %clinical 0
    R{C1_Prev0,1}=229;R{C1_Prev0,2}=240;R{C1_Prev0,3}=4;
    R{C1_CFR,1}=229;R{C1_CFR,2}=240;R{C1_CFR,3}=6;
    R{C1_Progr,1}=229;R{C1_Progr,2}=240;R{C1_Progr,3}=9;

    %clinical 1
    R{C2_Prev0,1}=251;R{C2_Prev0,2}=262;R{C2_Prev0,3}=4;
    R{C2_CFR,1}=251;R{C2_CFR,2}=262;R{C2_CFR,3}=6;
    R{C2_Progr,1}=251;R{C2_Progr,2}=262;R{C2_Progr,3}=9;

    %clinical 2
    R{C3_Prev0,1}=273;R{C3_Prev0,2}=284;R{C3_Prev0,3}=4;
    R{C3_CFR,1}=273;R{C3_CFR,2}=284;R{C3_CFR,3}=6;
    R{C3_Progr,1}=273;R{C3_Progr,2}=284;R{C3_Progr,3}=9;

    %clinical 3
    R{C4_Prev0,1}=295;R{C4_Prev0,2}=306;R{C4_Prev0,3}=4;
    R{C4_CFR,1}=295;R{C4_CFR,2}=306;R{C4_CFR,3}=6;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Interventions 7-CFR, 14 Function-DW
    R{I1_CFR,1}=339;R{I1_CFR,2}=339;R{I1_CFR,3}=7;R{I1_CFR,4}=14;
    R{I2_CFR,1}=361;R{I2_CFR,2}=361;R{I2_CFR,3}=7;R{I2_CFR,4}=14;
    R{I3_CFR,1}=383;R{I3_CFR,2}=383;R{I3_CFR,3}=7;R{I3_CFR,4}=14;
    R{I4_CFR,1}=405;R{I4_CFR,2}=405;R{I4_CFR,3}=7;R{I4_CFR,4}=14;

    R{I5_CFR_Function,1}=427;R{I5_CFR_Function,2}=427;R{I5_CFR_Function,3}=7;R{I5_CFR_Function,4}=14;
    R{I6_CFR_Function,1}=449;R{I6_CFR_Function,2}=449;R{I6_CFR_Function,3}=7;R{I6_CFR_Function,4}=14;

    R{I7_Function,1}=471;R{I7_Function,2}=471;R{I7_Function,4}=7;
    R{I8_Function,1}=493;R{I8_Function,2}=493;R{I8_Function,4}=7;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Screening Interventions
    %bi-annual screening
    R{I1_Screen,1}=515;R{I1_Screen,2}=515;R{I1_Screen,3}=7;R{I1_Screen,4}=14;R{I1_Screen,5}=21;R{I1_Screen,6}=28;
    R{I2_Screen,1}=537;R{I2_Screen,2}=537;R{I2_Screen,3}=7;R{I2_Screen,4}=14;R{I2_Screen,5}=21;R{I2_Screen,6}=28;
    %clinical examination
    R{I3_Screen,1}=559;R{I3_Screen,2}=559;R{I3_Screen,3}=7;R{I3_Screen,4}=14;R{I3_Screen,5}=21;R{I3_Screen,6}=28;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Disability weights
    %using row indices only
    R{D_C1,1}=137;R{D_C1,2}=225;
    R{D_C2,1}=159;R{D_C2,2}=247;
    R{D_C3,1}=181;R{D_C3,2}=269;
    R{D_C4,1}=203;R{D_C4,2}=291;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Current and scale-up coverage
    R{Tx_C1,1}=339;R{Tx_C1,2}=339;R{Tx_C1,3}=339;R{Tx_C1,4}=8;R{Tx_C1,5}=9;R{Tx_C1,6}=15;R{Tx_C1,7}=16;
    R{Tx_C2,1}=361;R{Tx_C2,2}=361;R{Tx_C2,3}=361;R{Tx_C2,4}=8;R{Tx_C2,5}=9;R{Tx_C2,6}=15;R{Tx_C2,7}=16;
    R{Tx_C3,1}=383;R{Tx_C3,2}=383;R{Tx_C3,3}=383;R{Tx_C3,4}=8;R{Tx_C3,5}=9;R{Tx_C3,6}=15;R{Tx_C3,7}=16;
    R{Tx_C4,1}=405;R{Tx_C4,2}=405;R{Tx_C4,3}=405;R{Tx_C4,4}=8;R{Tx_C4,5}=9;R{Tx_C4,6}=15;R{Tx_C4,7}=16;

    R{Tx_C5,1}=427;R{Tx_C5,2}=427;R{Tx_C5,3}=427;R{Tx_C5,4}=8;R{Tx_C5,5}=9;R{Tx_C5,6}=15;R{Tx_C5,7}=16;
    R{Tx_C6,1}=449;R{Tx_C6,2}=449;R{Tx_C6,3}=449;R{Tx_C6,4}=8;R{Tx_C6,5}=9;R{Tx_C6,6}=15;R{Tx_C6,7}=16;
    R{Tx_C7,1}=471;R{Tx_C7,2}=471;R{Tx_C7,3}=471;R{Tx_C7,4}=8;R{Tx_C7,5}=9;R{Tx_C7,6}=15;R{Tx_C7,7}=16;
    R{Tx_C8,1}=493;R{Tx_C8,2}=493;R{Tx_C8,3}=493;R{Tx_C8,4}=8;R{Tx_C8,5}=9;R{Tx_C8,6}=15;R{Tx_C8,7}=16;


    for I_f=1:12

    %onset to first stage
    lTEMPLATE{ R{onset,1}+I_f-1 ,R{onset,3}}=OnsetRates(I_f,1);    

    %progression
    lTEMPLATE{ R{PC1_Prev0,1}+I_f-1 ,R{PC1_Prev0,3} }=PrevS(I_f,1);
    lTEMPLATE{ R{PC2_Prev0,1}+I_f-1 ,R{PC2_Prev0,3} }=PrevS(I_f,2);
    lTEMPLATE{ R{PC3_Prev0,1}+I_f-1 ,R{PC3_Prev0,3} }=PrevS(I_f,3);
    lTEMPLATE{ R{PC4_Prev0,1}+I_f-1 ,R{PC4_Prev0,3} }=PrevS(I_f,4);

    %progression
    lTEMPLATE{ R{PC1_Progr,1}+I_f-1 ,R{PC1_Progr,3} }=DwellRate(I_f,1);
    lTEMPLATE{ R{PC2_Progr,1}+I_f-1 ,R{PC2_Progr,3} }=DwellRate(I_f,2);
    lTEMPLATE{ R{PC3_Progr,1}+I_f-1 ,R{PC3_Progr,3} }=DwellRate(I_f,3);
    %lTEMPLATE{ R{PC4_Progr,1}+I-1 ,R{PC4_Progr,3} }=DwellRate(I,4);

    %screening process
    lTEMPLATE{ R{PC1_Screen,1}+I_f-1 ,R{PC1_Screen,3} }=DiagnosticRate(I_f,1);
    lTEMPLATE{ R{PC2_Screen,1}+I_f-1 ,R{PC2_Screen,3} }=DiagnosticRate(I_f,2);
    lTEMPLATE{ R{PC3_Screen,1}+I_f-1 ,R{PC3_Screen,3} }=DiagnosticRate(I_f,3);
    lTEMPLATE{ R{PC4_Screen,1}+I_f-1 ,R{PC4_Screen,3} }=DiagnosticRate(I_f,4);

    %Cancer Mortality pre-Clinical
    lTEMPLATE{ R{PC1_CFR,1}+I_f-1 ,R{PC1_CFR,3} }=0;%MortRateDF(I,1);
    lTEMPLATE{ R{PC2_CFR,1}+I_f-1 ,R{PC2_CFR,3} }=0;%MortRateDF(I,1);
    lTEMPLATE{ R{PC3_CFR,1}+I_f-1 ,R{PC3_CFR,3} }=0;%MortRateDF(I,1);
    lTEMPLATE{ R{PC4_CFR,1}+I_f-1 ,R{PC4_CFR,3} }=MortRate(4,I_f);

    %Clinical
    %Base Cancer Mortality (after considering coverage of TX among those
    %diagnosed
    lTEMPLATE{ R{C1_CFR,1}+I_f-1 ,R{C1_CFR,3} }=MortRate(1,I_f);
    lTEMPLATE{ R{C2_CFR,1}+I_f-1 ,R{C2_CFR,3} }=MortRate(2,I_f);
    lTEMPLATE{ R{C3_CFR,1}+I_f-1 ,R{C3_CFR,3} }=MortRate(3,I_f);
    lTEMPLATE{ R{C4_CFR,1}+I_f-1 ,R{C4_CFR,3} }=MortRate(4,I_f);

    %prev 0
    lTEMPLATE{ R{C1_Prev0,1}+I_f-1 ,R{C1_Prev0,3} }=PrevC(I_f,1);
    lTEMPLATE{ R{C2_Prev0,1}+I_f-1 ,R{C2_Prev0,3} }=PrevC(I_f,2);
    lTEMPLATE{ R{C3_Prev0,1}+I_f-1 ,R{C3_Prev0,3} }=PrevC(I_f,3);
    lTEMPLATE{ R{C4_Prev0,1}+I_f-1 ,R{C4_Prev0,3} }=PrevC(I_f,4);

    %progression
    lTEMPLATE{ R{C1_Progr,1}+I_f-1 ,R{C1_Progr,3} }=0;%DwellRate(I,1);
    lTEMPLATE{ R{C2_Progr,1}+I_f-1 ,R{C2_Progr,3} }=0;%DwellRate(I,2);
    lTEMPLATE{ R{C3_Progr,1}+I_f-1 ,R{C3_Progr,3} }=0;%DwellRate(I,3);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Treatment impact on MORT
    lTEMPLATE{ R{I1_CFR,1}+I_f-1 ,R{I1_CFR,3} }=Impact(1,I_f);%INTVN{1,2}*100;
    lTEMPLATE{ R{I2_CFR,1}+I_f-1 ,R{I2_CFR,3} }=Impact(2,I_f);%INTVN{2,2}*100;
    lTEMPLATE{ R{I3_CFR,1}+I_f-1 ,R{I3_CFR,3} }=Impact(3,I_f);%INTVN{3,2}*100;
    lTEMPLATE{ R{I4_CFR,1}+I_f-1 ,R{I4_CFR,3} }=Impact(4,I_f);%INTVN{4,2}*100;

    % Treatment+palliative care impact on Mort
    lTEMPLATE{ R{I5_CFR_Function,1}+I_f-1 ,R{I5_CFR_Function,3} }=Impact(4,I_f);%INTVN{4,2}*100;
    lTEMPLATE{ R{I6_CFR_Function,1}+I_f-1 ,R{I6_CFR_Function,3} }=Impact(4,I_f);%INTVN{4,2}*100;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Baseline treatment coverage for treatment
    lTEMPLATE{ R{Tx_C1,1}+I_f-1 ,R{Tx_C1,4} }  = TxCoverage(1,I_f)*100;
    lTEMPLATE{ R{Tx_C2,1}+I_f-1 ,R{Tx_C2,4} }  = TxCoverage(2,I_f)*100;   
    lTEMPLATE{ R{Tx_C3,1}+I_f-1 ,R{Tx_C3,4} }  = TxCoverage(3,I_f)*100;
    lTEMPLATE{ R{Tx_C4,1}+I_f-1 ,R{Tx_C4,4} }  = TxCoverage(4,I_f)*100;

    lTEMPLATE{ R{Tx_C5,1}+I_f-1 ,R{Tx_C5,4} }  = TxCoverage(4,I_f)*100;
    lTEMPLATE{ R{Tx_C6,1}+I_f-1 ,R{Tx_C6,4} }  = TxCoverage(4,I_f)*100;
    lTEMPLATE{ R{Tx_C7,1}+I_f-1 ,R{Tx_C7,4} }  = TxCoverage(4,I_f)*100;
    lTEMPLATE{ R{Tx_C8,1}+I_f-1 ,R{Tx_C8,4} }  = TxCoverage(4,I_f)*100;


    %Default intervention treatment coverage= current + 80% of those not
    lTEMPLATE{ R{Tx_C1,1}+I_f-1 ,R{Tx_C1,5} }  = ((1-TxCoverage(1,I_f))*.8+TxCoverage(1,I_f)) *100;
    lTEMPLATE{ R{Tx_C2,1}+I_f-1 ,R{Tx_C2,5} }  = ((1-TxCoverage(2,I_f))*.8+TxCoverage(2,I_f)) *100;   
    lTEMPLATE{ R{Tx_C3,1}+I_f-1 ,R{Tx_C3,5} }  = ((1-TxCoverage(3,I_f))*.8+TxCoverage(3,I_f)) *100;
    lTEMPLATE{ R{Tx_C4,1}+I_f-1 ,R{Tx_C4,5} }  = ((1-TxCoverage(4,I_f))*.8+TxCoverage(4,I_f)) *100;

    lTEMPLATE{ R{Tx_C5,1}+I_f-1 ,R{Tx_C5,5} }  = ((1-TxCoverage(4,I_f))*.8+TxCoverage(4,I_f)) *100;
    lTEMPLATE{ R{Tx_C6,1}+I_f-1 ,R{Tx_C6,5} }  = ((1-TxCoverage(4,I_f))*.8+TxCoverage(4,I_f)) *100;
    lTEMPLATE{ R{Tx_C7,1}+I_f-1 ,R{Tx_C7,5} }  = ((1-TxCoverage(4,I_f))*.8+TxCoverage(4,I_f)) *100;
    lTEMPLATE{ R{Tx_C8,1}+I_f-1 ,R{Tx_C8,5} }  = ((1-TxCoverage(4,I_f))*.8+TxCoverage(4,I_f)) *100;

    % Baseline coverage for disability weight = same as Tx coverage
    lTEMPLATE{ R{Tx_C1,1}+I_f-1 ,R{Tx_C1,6} }  = TxCoverage(1,I_f)*100;
    lTEMPLATE{ R{Tx_C2,1}+I_f-1 ,R{Tx_C2,6} }  = TxCoverage(2,I_f)*100;   
    lTEMPLATE{ R{Tx_C3,1}+I_f-1 ,R{Tx_C3,6} }  = TxCoverage(3,I_f)*100;
    lTEMPLATE{ R{Tx_C4,1}+I_f-1 ,R{Tx_C4,6} }  = TxCoverage(4,I_f)*100;

    lTEMPLATE{ R{Tx_C5,1}+I_f-1 ,R{Tx_C5,6} }  = TxCoverage(4,I_f)*100;
    lTEMPLATE{ R{Tx_C6,1}+I_f-1 ,R{Tx_C6,6} }  = TxCoverage(4,I_f)*100;

    %Default intervention coverage for disability weight= current + 80% of those not
    lTEMPLATE{ R{Tx_C1,1}+I_f-1 ,R{Tx_C1,7} }  = ((1-TxCoverage(1,I_f))*.8+TxCoverage(1,I_f)) *100;
    lTEMPLATE{ R{Tx_C2,1}+I_f-1 ,R{Tx_C2,7} }  = ((1-TxCoverage(2,I_f))*.8+TxCoverage(2,I_f)) *100;   
    lTEMPLATE{ R{Tx_C3,1}+I_f-1 ,R{Tx_C3,7} }  = ((1-TxCoverage(3,I_f))*.8+TxCoverage(3,I_f)) *100;
    lTEMPLATE{ R{Tx_C4,1}+I_f-1 ,R{Tx_C4,7} }  = ((1-TxCoverage(4,I_f))*.8+TxCoverage(4,I_f)) *100;
    lTEMPLATE{ R{Tx_C5,1}+I_f-1 ,R{Tx_C5,7} }  = ((1-TxCoverage(4,I_f))*.8+TxCoverage(4,I_f)) *100;
    lTEMPLATE{ R{Tx_C6,1}+I_f-1 ,R{Tx_C6,7} }  = ((1-TxCoverage(4,I_f))*.8+TxCoverage(4,I_f)) *100;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Treatment
    %impact on disability weights, see def of DW in start of function
    lTEMPLATE{ R{I1_CFR,1}+I_f-1 ,R{I1_CFR,4} }=-(DW{1,1}-DW{1,2})/DW{1,1}*100;
    lTEMPLATE{ R{I2_CFR,1}+I_f-1 ,R{I2_CFR,4} }=-(DW{2,1}-DW{2,2})/DW{2,1}*100;
    lTEMPLATE{ R{I3_CFR,1}+I_f-1 ,R{I3_CFR,4} }=-(DW{3,1}-DW{3,2})/DW{3,1}*100;
    lTEMPLATE{ R{I4_CFR,1}+I_f-1 ,R{I4_CFR,4} }=-(DW{4,1}-DW{4,2})/DW{4,1}*100;

    %Treatment+ Palliative Care
    lTEMPLATE{ R{I5_CFR_Function,1}+I_f-1 ,R{I5_CFR_Function,4} }=-(DW{4,1}-DW{5,2})/DW{4,1}*100;
    lTEMPLATE{ R{I6_CFR_Function,1}+I_f-1 ,R{I6_CFR_Function,4} }=-(DW{4,1}-DW{6,2})/DW{4,1}*100;

    %Palliative Care
    lTEMPLATE{ R{I7_Function,1}+I_f-1 ,R{I7_Function,4} }=-(DW{4,1}-DW{5,2})/DW{4,1}*100;
    lTEMPLATE{ R{I8_Function,1}+I_f-1 ,R{I8_Function,4} }=-(DW{4,1}-DW{6,2})/DW{4,1}*100;


    %Bi annual screening of 50-69, age match not perfect
    lTEMPLATE{ R{I1_Screen,1}+I_f-1 ,R{I1_Screen,3} }=100*0/2*0.9;
    lTEMPLATE{ R{I1_Screen,1}+I_f-1 ,R{I1_Screen,4} }=100*0/2*0.9;
    lTEMPLATE{ R{I1_Screen,1}+I_f-1 ,R{I1_Screen,5} }=100*0/2*0.9;
    lTEMPLATE{ R{I1_Screen,1}+I_f-1 ,R{I1_Screen,6} }=100*0/2*0.9;
    if(ismember(I_f,[9,10]))
    lTEMPLATE{ R{I1_Screen,1}+I_f-1 ,R{I1_Screen,3} }=100*1/2*0.9;
    lTEMPLATE{ R{I1_Screen,1}+I_f-1 ,R{I1_Screen,4} }=100*1/2*0.9;
    lTEMPLATE{ R{I1_Screen,1}+I_f-1 ,R{I1_Screen,5} }=100*1/2*0.9;
    lTEMPLATE{ R{I1_Screen,1}+I_f-1 ,R{I1_Screen,6} }=100*1/2*0.9;
    end

    %Bi annual screening of 40-69, age match not perfect
    lTEMPLATE{ R{I2_Screen,1}+I_f-1 ,R{I2_Screen,3} }=100*0/2*0.9;
    lTEMPLATE{ R{I2_Screen,1}+I_f-1 ,R{I2_Screen,4} }=100*0/2*0.9;
    lTEMPLATE{ R{I2_Screen,1}+I_f-1 ,R{I2_Screen,5} }=100*0/2*0.9;
    lTEMPLATE{ R{I2_Screen,1}+I_f-1 ,R{I2_Screen,6} }=100*0/2*0.9;
    if(ismember(I_f,[8,9,10]))
    lTEMPLATE{ R{I2_Screen,1}+I_f-1 ,R{I2_Screen,3} }=100*1/2*0.9;
    lTEMPLATE{ R{I2_Screen,1}+I_f-1 ,R{I2_Screen,4} }=100*1/2*0.9;
    lTEMPLATE{ R{I2_Screen,1}+I_f-1 ,R{I2_Screen,5} }=100*1/2*0.9;
    lTEMPLATE{ R{I2_Screen,1}+I_f-1 ,R{I2_Screen,6} }=100*1/2*0.9;
    end

    %Screening-by clinical exam
    lTEMPLATE{ R{I3_Screen,1}+I_f-1 ,R{I3_Screen,3} }=100*0/2*0.35;
    lTEMPLATE{ R{I3_Screen,1}+I_f-1 ,R{I3_Screen,4} }=100*0/2*0.35;
    lTEMPLATE{ R{I3_Screen,1}+I_f-1 ,R{I3_Screen,5} }=100*0/2*0.35;
    lTEMPLATE{ R{I3_Screen,1}+I_f-1 ,R{I3_Screen,6} }=100*0/2*0.35;
    if(ismember(I_f,[8,9,10]))%assume it takes longer than bi-anual
    lTEMPLATE{ R{I3_Screen,1}+I_f-1 ,R{I3_Screen,3} }=100*1/4*0.35;
    lTEMPLATE{ R{I3_Screen,1}+I_f-1 ,R{I3_Screen,4} }=100*1/4*0.35;
    lTEMPLATE{ R{I3_Screen,1}+I_f-1 ,R{I3_Screen,5} }=100*1/4*0.35;
    lTEMPLATE{ R{I3_Screen,1}+I_f-1 ,R{I3_Screen,6} }=100*1/4*0.35;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %Disability weights
    %preclinical
    lTEMPLATE{ R{D_C1,1} , 3 }=DW{1,1};
    lTEMPLATE{ R{D_C2,1} , 3 }=DW{2,1};
    lTEMPLATE{ R{D_C3,1} , 3 }=DW{3,1};
    lTEMPLATE{ R{D_C4,1} , 3 }=DW{4,1};
    %clinical
    lTEMPLATE{ R{D_C1,2} , 3 }=DW{1,1};
    lTEMPLATE{ R{D_C2,2} , 3 }=DW{2,1};
    lTEMPLATE{ R{D_C3,2} , 3 }=DW{3,1};
    lTEMPLATE{ R{D_C4,2} , 3 }=DW{4,1};

end

end

end




