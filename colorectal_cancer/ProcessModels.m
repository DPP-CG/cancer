function ProcessModelsMike

    [~,~,TEMPLATE]=xlsread('ColorectalTemplate.xls','ColorectalTemplate');
    [~,~,INTERVENTIONS]=xlsread('ColorectalTemplate.xls','Interventions');
    [~,~,DWeights]=xlsread('ColorectalTemplate.xls','DW');

    %Intervention name

    INTVN{1,1}=INTERVENTIONS{6,2};
    INTVN{2,1}=INTERVENTIONS{8,2};
    INTVN{3,1}=INTERVENTIONS{10,2};
    INTVN{4,1}=INTERVENTIONS{12,2};

    %MORT
    INTVN{1,2}=INTERVENTIONS{6,5};
    INTVN{2,2}=INTERVENTIONS{8,5};
    INTVN{3,2}=INTERVENTIONS{10,5};
    INTVN{4,2}=INTERVENTIONS{12,5};

    %DW
    DW{1,1}=DWeights{1,2};
    DW{2,1}=DWeights{2,2};
    DW{3,1}=DWeights{3,2};
    DW{4,1}=DWeights{4,2};

    %DW
    DW{1,2}=DWeights{1,3};
    DW{2,2}=DWeights{2,3};
    DW{3,2}=DWeights{3,3};
    DW{4,2}=DWeights{4,3};


    %pal care
    DW{5,2}=DWeights{7,3};
    DW{6,2}=DWeights{8,3};

    %pal care
    DW{5,1}=DW{4,1};
    DW{6,1}=DW{4,1};

    [~,~,RAW_D]=xlsread('res4.xlsx');

    %%%%%%%%%%%%%%%%%%%%
    %%% MIKE'S INPUT %%%
    %%%%%%%%%%%%%%%%%%%%
    
%% Following loop is to do the calculations separately for Male and Female
sex = ["Female","Male"];
formatF = 'DataInci-%s.xlsx';
for s_loop = 1:2
    
    % cancer data for sex under consideration
    f_name = sprintf(formatF,sex(s_loop));
    [~,~,FORMAT]=xlsread(f_name);

    %%
    P = (cell2mat(FORMAT(4:13,2:22)))'; %incidence
    S =(cell2mat(FORMAT(32:35,2:22)))'; %distribution of stages
    M = (cell2mat(FORMAT(18:27,2:22)))'; %mortality
    cancerMortalityArrayMatrix = M / 1000;
    prevalenceActualMatrix=P / 1000; %prevalence = incidence/1000
    StageDistMatrix = S;
    for i=1:21 % i = 21 because we just have data columns for 16 region
        Countries{i}=cell2mat(FORMAT(3,i+1));   
    end
    AgeArrayL=[1 15 40 45 50  55	60	65 70 75];
    AgeArrayU=[14 39 44  49 54	59	64	69 74 100];
    
    %%

    % I = 12 for SEARB region and, I = 2 for AFRE
     % I = 16 for peru 
     % I = 17 for bolivia
     % I = 18 for eucador
     % I = 19 for colombia
     % I = 20 for Brazil 
     
    I = 16;
    
    %%
        Countries{I};
        ind=find(strcmp(RAW_D(:,1),Countries{I})); %finding a country index string at the res4.xlx file
        DATA=RAW_D(ind,6);
        pA=cell2mat(DATA(163:243))'; % population by age
        Mx=cell2mat(DATA(406:486))'; % natural mortality by age
        Births=cell2mat(DATA(487:487)); %New birth population proportion
        xa=101;
        pA(end+1:101)=0; %adding ages till 101
        Mx(end+1:101)=0; %adding ages till 101

    %%    
        PAR{1}=pA/sum(pA); %distribtion of each age
        PAR{2}=Mx; %mortality by age
        PAR{3}=Births/sum(pA); %percentage of births compare to total population
        PAR{4}=prevalenceActualMatrix(I,:); %prevalence of the country by age (14 years range)
        PAR{5}=cancerMortalityArrayMatrix(I,:); %mortality by age (14 years range)
        PAR{6}=StageDistMatrix(I,:); %distribution of stages (insitu, local, regional, distant)
    %   PAR{8}=AgeArrayU(I,:);
        PAR{7}=AgeArrayL(1,:); %lower range of ages
        PAR{8}=AgeArrayU(1,:); %upper range of ages
        PAR{9}=I; %country number
        PAR{10}=Countries(I); %country name
        
    %%
        [RES]=MarkovChainCountry2(PAR); %storing markov chain result 
        [C_TEMPLATE]=WriteTemplate(RES,TEMPLATE); %writing result from markov chain to template
        formatSpec = '%s-ProcessModOut--OLD--ColorectalCancer-%s.xlsx';
        str = sprintf(formatSpec,datetime('today','Format','yyyy-MM-dd'),Countries{I});
        xlswrite('ProcModOut_old',C_TEMPLATE)
        
        %% producing file for spectrum input

        % reading the association file
        [~,~,data_asso] = xlsread('ColorectalCancer_Asso',s_loop);

        % reading the process model output
        [~,~,data2] = xlsread('ProcModOut_old');

        % writing file for spectrum input
        data2(310:710,:) = data_asso(310:end,:);
        xlswrite(str,data2,sex(s_loop));
        
        formatS = ' ---- results are written for %s -----';
        which_result = sprintf(formatS,sex(s_loop));
        disp(which_result);
        
end

        %Write Template  
        function [lTEMPLATE]=WriteTemplate(RES,TEMPLATE)

        lTEMPLATE=TEMPLATE;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Data from Markov and Cohort model

        PrevS=RES{1};
        PrevC=RES{2};
        Incidence=RES{3};
        OnsetRates=RES{4};
        OnsetRates(end)=OnsetRates(end-1);
        PolyponsetRateAge = RES{5};
        DiagnosticRate=1./RES{6};
        DwellRate=1./RES{7};

        for J=1:12
            for K=1:4
                if  DwellRate(J,K)>10000
                    DwellRate(J,K)=0;
                end
            end
        end

        MortRate= RES{8};%Relative Mortality without treatment
        MortRateDF= RES{9};% Disease free mortality
        Impact = RES{10};%Relative Mortality with treatment
        TxCoverage = RES{11};% Current coverage of treatment
        PrevPolys = RES{12}; %Prevalance of Polyps

        PrevS(1:3,:)= 0;    %Forcing Prevalance of Preclinical Cancer to zero for first 3 agegroups
        PrevC(1:3,:) =0;    %Forcing Prevalance of Cancer to zero for first 3 agegroups
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

        %Pre-cancer Stages
        %Polyp
        PCPolypLT5mm_Prev0=51;PCPolypLT5mm_CFR=52;PCPolypLT5mm_Progr=53;
        PCPolyp6to10mm_Prev0=54;PCPolyp6to10mm_CFR=55;PCPolyp6to10mm_Progr=56;
        PCPolypGT10mm_Prev0=57;PCPolypGT10mm_CFR=58;PCPolypGT10mm_Progr=59;

        %disease free to PolypGT10mmm
        DF_PC1 =60;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %disease free
        R{onset,1}=9;R{onset,2}=20;R{onset,3}=7;
        R{DF_PC1,1}=9;R{DF_PC1,2}=20;R{DF_PC1,3}=10;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %pre-cancer
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %PolypLT5mm
        R{PCPolypLT5mm_Prev0,1}=75;R{PCPolypLT5mm_Prev0,2}=86;R{PCPolypLT5mm_Prev0,3}=4;
        R{PCPolypLT5mm_CFR,1}=75;R{PCPolypLT5mm_CFR,2}=86;R{PCPolypLT5mm_CFR,3}=6;
        R{PCPolypLT5mm_Progr,1}=75;R{PCPolypLT5mm_Progr,2}=86;R{PCPolypLT5mm_Progr,3}=9;

        %Polyp6to10mm
        R{PCPolyp6to10mm_Prev0,1}=97;R{PCPolyp6to10mm_Prev0,2}=108;R{PCPolyp6to10mm_Prev0,3}=4;
        R{PCPolyp6to10mm_CFR,1}=97;R{PCPolyp6to10mm_CFR,2}=108;R{PCPolyp6to10mm_CFR,3}=6;
        R{PCPolyp6to10mm_Progr,1}=97;R{PCPolyp6to10mm_Progr,2}=108;R{PCPolyp6to10mm_Progr,3}=9;

        %PolypGT10mm
        R{PCPolypGT10mm_Prev0,1}=119;R{PCPolypGT10mm_Prev0,2}=130;R{PCPolypGT10mm_Prev0,3}=4;
        R{PCPolypGT10mm_CFR,1}=119;R{PCPolypGT10mm_CFR,2}=130;R{PCPolypGT10mm_CFR,3}=6;
        R{PCPolypGT10mm_Progr,1}=119;R{PCPolypGT10mm_Progr,2}=130;R{PCPolypGT10mm_Progr,3}=9;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %Cancer stages
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %pre clinical
        %stage 0
        R{PC1_Prev0,1}=141;R{PC1_Prev0,2}=152;R{PC1_Prev0,3}=4; %Baseline Prev
        R{PC1_CFR,1}=141;R{PC1_CFR,2}=152;R{PC1_CFR,3}=6; %CFR rate
        R{PC1_Screen,1}=141;R{PC1_Screen,2}=152;R{PC1_Screen,3}=9; %Screened rate
        R{PC1_Progr,1}=141;R{PC1_Progr,2}=152;R{PC1_Progr,3}=12; %Pre-cancer 1 rate

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
        R{I1_CFR,1}=405;R{I1_CFR,2}=405;R{I1_CFR,3}=7;R{I1_CFR,4}=14;
        R{I2_CFR,1}=427;R{I2_CFR,2}=427;R{I2_CFR,3}=7;R{I2_CFR,4}=14;
        R{I3_CFR,1}=449;R{I3_CFR,2}=449;R{I3_CFR,3}=7;R{I3_CFR,4}=14;
        R{I4_CFR,1}=471;R{I4_CFR,2}=471;R{I4_CFR,3}=7;R{I4_CFR,4}=14;

        R{I5_CFR_Function,1}=493;R{I5_CFR_Function,2}=493;R{I5_CFR_Function,3}=7;R{I5_CFR_Function,4}=14;
        R{I6_CFR_Function,1}=515;R{I6_CFR_Function,2}=515;R{I6_CFR_Function,3}=7;R{I6_CFR_Function,4}=14;

        R{I7_Function,1}=537;R{I7_Function,2}=537;R{I7_Function,4}=7;
        R{I8_Function,1}=559;R{I8_Function,2}=559;R{I8_Function,4}=7;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Screening Interventions
        %bi-annual screening
        R{I1_Screen,1}=651;R{I1_Screen,2}=651;R{I1_Screen,3}=7;R{I1_Screen,4}=14;R{I1_Screen,5}=21;R{I1_Screen,6}=28;
        %R{I2_Screen,1}=605;R{I2_Screen,2}=605;R{I2_Screen,3}=7;R{I2_Screen,4}=14;R{I2_Screen,5}=21;R{I2_Screen,6}=28;
        %clinical examination
        %R{I3_Screen,1}=628;R{I3_Screen,2}=628;R{I3_Screen,3}=7;R{I3_Screen,4}=14;R{I3_Screen,5}=21;R{I3_Screen,6}=28;
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
        R{Tx_C1,1}=405;R{Tx_C1,2}=405;R{Tx_C1,3}=405;R{Tx_C1,4}=8;R{Tx_C1,5}=9;R{Tx_C1,6}=15;R{Tx_C1,7}=16;
        R{Tx_C2,1}=427;R{Tx_C2,2}=427;R{Tx_C2,3}=427;R{Tx_C2,4}=8;R{Tx_C2,5}=9;R{Tx_C2,6}=15;R{Tx_C2,7}=16;
        R{Tx_C3,1}=449;R{Tx_C3,2}=449;R{Tx_C3,3}=449;R{Tx_C3,4}=8;R{Tx_C3,5}=9;R{Tx_C3,6}=15;R{Tx_C3,7}=16;
        R{Tx_C4,1}=471;R{Tx_C4,2}=471;R{Tx_C4,3}=471;R{Tx_C4,4}=8;R{Tx_C4,5}=9;R{Tx_C4,6}=15;R{Tx_C4,7}=16;

        R{Tx_C5,1}=493;R{Tx_C5,2}=493;R{Tx_C5,3}=493;R{Tx_C5,4}=8;R{Tx_C5,5}=9;R{Tx_C5,6}=15;R{Tx_C5,7}=16;
        R{Tx_C6,1}=515;R{Tx_C6,2}=515;R{Tx_C6,3}=515;R{Tx_C6,4}=8;R{Tx_C6,5}=9;R{Tx_C6,6}=15;R{Tx_C6,7}=16;
        R{Tx_C7,1}=537;R{Tx_C7,2}=537;R{Tx_C7,3}=537;R{Tx_C7,4}=8;R{Tx_C7,5}=9;R{Tx_C7,6}=15;R{Tx_C7,7}=16;
        R{Tx_C8,1}=559;R{Tx_C8,2}=559;R{Tx_C8,3}=559;R{Tx_C8,4}=8;R{Tx_C8,5}=9;R{Tx_C8,6}=15;R{Tx_C8,7}=16;


        for I=1:12

        %onset to first stage
        lTEMPLATE{ R{onset,1}+I-1 ,R{onset,3}}=PolyponsetRateAge(I,1); 

        %healthy to pre-clinical
        lTEMPLATE{ R{DF_PC1,1}+I-1 ,R{DF_PC1,3}}=PolyponsetRateAge(I,5);    


        %progression
        lTEMPLATE{ R{PC1_Prev0,1}+I-1 ,R{PC1_Prev0,3} }=PrevS(I,1);
        lTEMPLATE{ R{PC2_Prev0,1}+I-1 ,R{PC2_Prev0,3} }=PrevS(I,2);
        lTEMPLATE{ R{PC3_Prev0,1}+I-1 ,R{PC3_Prev0,3} }=PrevS(I,3);
        lTEMPLATE{ R{PC4_Prev0,1}+I-1 ,R{PC4_Prev0,3} }=PrevS(I,4);

        %progression
        lTEMPLATE{ R{PC1_Progr,1}+I-1 ,R{PC1_Progr,3} }=DwellRate(I,1);
        lTEMPLATE{ R{PC2_Progr,1}+I-1 ,R{PC2_Progr,3} }=DwellRate(I,2);
        lTEMPLATE{ R{PC3_Progr,1}+I-1 ,R{PC3_Progr,3} }=DwellRate(I,3);
        %lTEMPLATE{ R{PC4_Progr,1}+I-1 ,R{PC4_Progr,3} }=DwellRate(I,4);

        %screening process
        lTEMPLATE{ R{PC1_Screen,1}+I-1 ,R{PC1_Screen,3} }=DiagnosticRate(I,1);
        lTEMPLATE{ R{PC2_Screen,1}+I-1 ,R{PC2_Screen,3} }=DiagnosticRate(I,2);
        lTEMPLATE{ R{PC3_Screen,1}+I-1 ,R{PC3_Screen,3} }=DiagnosticRate(I,3);
        lTEMPLATE{ R{PC4_Screen,1}+I-1 ,R{PC4_Screen,3} }=DiagnosticRate(I,4);

        %Cancer Mortality pre-Clinical
        lTEMPLATE{ R{PC1_CFR,1}+I-1 ,R{PC1_CFR,3} }=0;%MortRateDF(I,1);
        lTEMPLATE{ R{PC2_CFR,1}+I-1 ,R{PC2_CFR,3} }=0;%MortRateDF(I,1);
        lTEMPLATE{ R{PC3_CFR,1}+I-1 ,R{PC3_CFR,3} }=0;%MortRateDF(I,1);
        lTEMPLATE{ R{PC4_CFR,1}+I-1 ,R{PC4_CFR,3} }=MortRate(4,I);

        %Clinical
        %Base Cancer Mortality (after considering coverage of TX among those
        %diagnosed
        lTEMPLATE{ R{C1_CFR,1}+I-1 ,R{C1_CFR,3} }=MortRate(1,I);
        lTEMPLATE{ R{C2_CFR,1}+I-1 ,R{C2_CFR,3} }=MortRate(2,I);
        lTEMPLATE{ R{C3_CFR,1}+I-1 ,R{C3_CFR,3} }=MortRate(3,I);
        lTEMPLATE{ R{C4_CFR,1}+I-1 ,R{C4_CFR,3} }=MortRate(4,I);

        %prev 0
        lTEMPLATE{ R{C1_Prev0,1}+I-1 ,R{C1_Prev0,3} }=PrevC(I,1);
        lTEMPLATE{ R{C2_Prev0,1}+I-1 ,R{C2_Prev0,3} }=PrevC(I,2);
        lTEMPLATE{ R{C3_Prev0,1}+I-1 ,R{C3_Prev0,3} }=PrevC(I,3);
        lTEMPLATE{ R{C4_Prev0,1}+I-1 ,R{C4_Prev0,3} }=PrevC(I,4);

        %progression
        lTEMPLATE{ R{C1_Progr,1}+I-1 ,R{C1_Progr,3} }=0;%DwellRate(I,1);
        lTEMPLATE{ R{C2_Progr,1}+I-1 ,R{C2_Progr,3} }=0;%DwellRate(I,2);
        lTEMPLATE{ R{C3_Progr,1}+I-1 ,R{C3_Progr,3} }=0;%DwellRate(I,3);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Treatment impact on MORT
        lTEMPLATE{ R{I1_CFR,1}+I-1 ,R{I1_CFR,3} }=Impact(1,I);%INTVN{1,2}*100;
        lTEMPLATE{ R{I2_CFR,1}+I-1 ,R{I2_CFR,3} }=Impact(2,I);%INTVN{2,2}*100;
        lTEMPLATE{ R{I3_CFR,1}+I-1 ,R{I3_CFR,3} }=Impact(3,I);%INTVN{3,2}*100;
        lTEMPLATE{ R{I4_CFR,1}+I-1 ,R{I4_CFR,3} }=Impact(4,I);%INTVN{4,2}*100;

        % Treatment+palliative care impact on Mort
        lTEMPLATE{ R{I5_CFR_Function,1}+I-1 ,R{I5_CFR_Function,3} }=Impact(4,I);%INTVN{4,2}*100;
        lTEMPLATE{ R{I6_CFR_Function,1}+I-1 ,R{I6_CFR_Function,3} }=Impact(4,I);%INTVN{4,2}*100;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Baseline treatment coverage for treatment
        lTEMPLATE{ R{Tx_C1,1}+I-1 ,R{Tx_C1,4} }  = TxCoverage(1,I)*100;
        lTEMPLATE{ R{Tx_C2,1}+I-1 ,R{Tx_C2,4} }  = TxCoverage(2,I)*100;   
        lTEMPLATE{ R{Tx_C3,1}+I-1 ,R{Tx_C3,4} }  = TxCoverage(3,I)*100;
        lTEMPLATE{ R{Tx_C4,1}+I-1 ,R{Tx_C4,4} }  = TxCoverage(4,I)*100;

        lTEMPLATE{ R{Tx_C5,1}+I-1 ,R{Tx_C5,4} }  = TxCoverage(4,I)*100;
        lTEMPLATE{ R{Tx_C6,1}+I-1 ,R{Tx_C6,4} }  = TxCoverage(4,I)*100;
        lTEMPLATE{ R{Tx_C7,1}+I-1 ,R{Tx_C7,4} }  = TxCoverage(4,I)*100;
        lTEMPLATE{ R{Tx_C8,1}+I-1 ,R{Tx_C8,4} }  = TxCoverage(4,I)*100;


        %Default intervention treatment coverage= current + 80% of those not
        lTEMPLATE{ R{Tx_C1,1}+I-1 ,R{Tx_C1,5} }  = ((1-TxCoverage(1,I))*.8+TxCoverage(1,I)) *100;
        lTEMPLATE{ R{Tx_C2,1}+I-1 ,R{Tx_C2,5} }  = ((1-TxCoverage(2,I))*.8+TxCoverage(2,I)) *100;   
        lTEMPLATE{ R{Tx_C3,1}+I-1 ,R{Tx_C3,5} }  = ((1-TxCoverage(3,I))*.8+TxCoverage(3,I)) *100;
        lTEMPLATE{ R{Tx_C4,1}+I-1 ,R{Tx_C4,5} }  = ((1-TxCoverage(4,I))*.8+TxCoverage(4,I)) *100;

        lTEMPLATE{ R{Tx_C5,1}+I-1 ,R{Tx_C5,5} }  = ((1-TxCoverage(4,I))*.8+TxCoverage(4,I)) *100;
        lTEMPLATE{ R{Tx_C6,1}+I-1 ,R{Tx_C6,5} }  = ((1-TxCoverage(4,I))*.8+TxCoverage(4,I)) *100;
        lTEMPLATE{ R{Tx_C7,1}+I-1 ,R{Tx_C7,5} }  = ((1-TxCoverage(4,I))*.8+TxCoverage(4,I)) *100;
        lTEMPLATE{ R{Tx_C8,1}+I-1 ,R{Tx_C8,5} }  = ((1-TxCoverage(4,I))*.8+TxCoverage(4,I)) *100;

        % Baseline coverage for disability weight = same as Tx coverage
        lTEMPLATE{ R{Tx_C1,1}+I-1 ,R{Tx_C1,6} }  = TxCoverage(1,I)*100;
        lTEMPLATE{ R{Tx_C2,1}+I-1 ,R{Tx_C2,6} }  = TxCoverage(2,I)*100;   
        lTEMPLATE{ R{Tx_C3,1}+I-1 ,R{Tx_C3,6} }  = TxCoverage(3,I)*100;
        lTEMPLATE{ R{Tx_C4,1}+I-1 ,R{Tx_C4,6} }  = TxCoverage(4,I)*100;

        lTEMPLATE{ R{Tx_C5,1}+I-1 ,R{Tx_C5,6} }  = TxCoverage(4,I)*100;
        lTEMPLATE{ R{Tx_C6,1}+I-1 ,R{Tx_C6,6} }  = TxCoverage(4,I)*100;

        %Default intervention coverage for disability weight= current + 80% of those not
        lTEMPLATE{ R{Tx_C1,1}+I-1 ,R{Tx_C1,7} }  = ((1-TxCoverage(1,I))*.8+TxCoverage(1,I)) *100;
        lTEMPLATE{ R{Tx_C2,1}+I-1 ,R{Tx_C2,7} }  = ((1-TxCoverage(2,I))*.8+TxCoverage(2,I)) *100;   
        lTEMPLATE{ R{Tx_C3,1}+I-1 ,R{Tx_C3,7} }  = ((1-TxCoverage(3,I))*.8+TxCoverage(3,I)) *100;
        lTEMPLATE{ R{Tx_C4,1}+I-1 ,R{Tx_C4,7} }  = ((1-TxCoverage(4,I))*.8+TxCoverage(4,I)) *100;
        lTEMPLATE{ R{Tx_C5,1}+I-1 ,R{Tx_C5,7} }  = ((1-TxCoverage(4,I))*.8+TxCoverage(4,I)) *100;
        lTEMPLATE{ R{Tx_C6,1}+I-1 ,R{Tx_C6,7} }  = ((1-TxCoverage(4,I))*.8+TxCoverage(4,I)) *100;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Treatment
        %impact on disability weights, see def of DW in start of function
        lTEMPLATE{ R{I1_CFR,1}+I-1 ,R{I1_CFR,4} }=-(DW{1,1}-DW{1,2})/DW{1,1}*100;
        lTEMPLATE{ R{I2_CFR,1}+I-1 ,R{I2_CFR,4} }=-(DW{2,1}-DW{2,2})/DW{2,1}*100;
        lTEMPLATE{ R{I3_CFR,1}+I-1 ,R{I3_CFR,4} }=-(DW{3,1}-DW{3,2})/DW{3,1}*100;
        lTEMPLATE{ R{I4_CFR,1}+I-1 ,R{I4_CFR,4} }=-(DW{4,1}-DW{4,2})/DW{4,1}*100;

        %Treatment+ Palliative Care
        lTEMPLATE{ R{I5_CFR_Function,1}+I-1 ,R{I5_CFR_Function,4} }=-(DW{4,1}-DW{5,2})/DW{4,1}*100;
        lTEMPLATE{ R{I6_CFR_Function,1}+I-1 ,R{I6_CFR_Function,4} }=-(DW{4,1}-DW{6,2})/DW{4,1}*100;

        %Palliative Care
        lTEMPLATE{ R{I7_Function,1}+I-1 ,R{I7_Function,4} }=-(DW{4,1}-DW{5,2})/DW{4,1}*100;
        lTEMPLATE{ R{I8_Function,1}+I-1 ,R{I8_Function,4} }=-(DW{4,1}-DW{6,2})/DW{4,1}*100;


        % %Bi annual screening of 50-69, age match not perfect
        % lTEMPLATE{ R{I1_Screen,1}+I-1 ,R{I1_Screen,3} }=100*0/2*0.9;
        % lTEMPLATE{ R{I1_Screen,1}+I-1 ,R{I1_Screen,4} }=100*0/2*0.9;
        % lTEMPLATE{ R{I1_Screen,1}+I-1 ,R{I1_Screen,5} }=100*0/2*0.9;
        % lTEMPLATE{ R{I1_Screen,1}+I-1 ,R{I1_Screen,6} }=100*0/2*0.9;
        % if(ismember(I,[9,10]))
        % lTEMPLATE{ R{I1_Screen,1}+I-1 ,R{I1_Screen,3} }=100*1/2*0.9;
        % lTEMPLATE{ R{I1_Screen,1}+I-1 ,R{I1_Screen,4} }=100*1/2*0.9;
        % lTEMPLATE{ R{I1_Screen,1}+I-1 ,R{I1_Screen,5} }=100*1/2*0.9;
        % lTEMPLATE{ R{I1_Screen,1}+I-1 ,R{I1_Screen,6} }=100*1/2*0.9;
        % end
        % 
        % %Bi annual screening of 40-69, age match not perfect
        % lTEMPLATE{ R{I2_Screen,1}+I-1 ,R{I2_Screen,3} }=100*0/2*0.9;
        % lTEMPLATE{ R{I2_Screen,1}+I-1 ,R{I2_Screen,4} }=100*0/2*0.9;
        % lTEMPLATE{ R{I2_Screen,1}+I-1 ,R{I2_Screen,5} }=100*0/2*0.9;
        % lTEMPLATE{ R{I2_Screen,1}+I-1 ,R{I2_Screen,6} }=100*0/2*0.9;
        % if(ismember(I,[8,9,10]))
        % lTEMPLATE{ R{I2_Screen,1}+I-1 ,R{I2_Screen,3} }=100*1/2*0.9;
        % lTEMPLATE{ R{I2_Screen,1}+I-1 ,R{I2_Screen,4} }=100*1/2*0.9;
        % lTEMPLATE{ R{I2_Screen,1}+I-1 ,R{I2_Screen,5} }=100*1/2*0.9;
        % lTEMPLATE{ R{I2_Screen,1}+I-1 ,R{I2_Screen,6} }=100*1/2*0.9;
        % end
        % 
        % %Screening-by clinical exam
        % lTEMPLATE{ R{I3_Screen,1}+I-1 ,R{I3_Screen,3} }=100*0/2*0.35;
        % lTEMPLATE{ R{I3_Screen,1}+I-1 ,R{I3_Screen,4} }=100*0/2*0.35;
        % lTEMPLATE{ R{I3_Screen,1}+I-1 ,R{I3_Screen,5} }=100*0/2*0.35;
        % lTEMPLATE{ R{I3_Screen,1}+I-1 ,R{I3_Screen,6} }=100*0/2*0.35;
        % if(ismember(I,[8,9,10]))%assume it takes longer than bi-anual
        % lTEMPLATE{ R{I3_Screen,1}+I-1 ,R{I3_Screen,3} }=100*1/4*0.35;
        % lTEMPLATE{ R{I3_Screen,1}+I-1 ,R{I3_Screen,4} }=100*1/4*0.35;
        % lTEMPLATE{ R{I3_Screen,1}+I-1 ,R{I3_Screen,5} }=100*1/4*0.35;
        % lTEMPLATE{ R{I3_Screen,1}+I-1 ,R{I3_Screen,6} }=100*1/4*0.35;
        % end

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


        %pre-cancer stages
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %onset to first stage
        %lTEMPLATE{ R{onset,1}+I-1 ,R{onset,3}}=OnsetRates(I,1);

        %Baseline Prevalence
        lTEMPLATE{ R{PCPolypLT5mm_Prev0,1}+I-1 ,R{PCPolypLT5mm_Prev0,3} }=PrevPolys(I,1);
        lTEMPLATE{ R{PCPolyp6to10mm_Prev0,1}+I-1 ,R{PCPolyp6to10mm_Prev0,3} }=PrevPolys(I,2);
        lTEMPLATE{ R{PCPolypGT10mm_Prev0,1}+I-1 ,R{PCPolypGT10mm_Prev0,3} }=PrevPolys(I,3);

        %CFR rate
        lTEMPLATE{ R{PCPolypLT5mm_CFR,1}+I-1 ,R{PCPolypLT5mm_CFR,3} }=0;
        lTEMPLATE{ R{PCPolyp6to10mm_CFR,1}+I-1 ,R{PCPolyp6to10mm_CFR,3} }=0;
        lTEMPLATE{ R{PCPolypGT10mm_CFR,1}+I-1 ,R{PCPolypGT10mm_CFR,3} }=0;

        %Progression
        lTEMPLATE{ R{PCPolypLT5mm_Progr,1}+I-1 ,R{PCPolypLT5mm_Progr,3} }=0.021;
        lTEMPLATE{ R{PCPolyp6to10mm_Progr,1}+I-1 ,R{PCPolyp6to10mm_Progr,3} }=0.057;
        lTEMPLATE{ R{PCPolypGT10mm_Progr,1}+I-1 ,R{PCPolypGT10mm_Progr,3} }=0.063;


        end

        end

end
