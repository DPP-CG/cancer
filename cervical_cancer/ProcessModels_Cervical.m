function ProcessModels_Cervical


[~,~,TEMPLATEW]=xlsread('CervicalCancerTemplate.xls','CervicalCancerTemplate');
[~,~,TEMPLATEM]=xlsread('CervicalCancerMen.xlsx');
[~,~,INTERVENTIONS]=xlsread('CervicalCancerTemplate.xls','Interventions');
[~,~,DWeights]=xlsread('CervicalCancerTemplate.xls','DW');

%[~,~,TEMPLATEW]=xlsread('CervicalCancerTemplate-2016-02-03.xls','CervicalCancerTemplate');
%[~,~,TEMPLATEM]=xlsread('CervicalCancerMen-2016-02-03.xlsx');
%[~,~,INTERVENTIONS]=xlsread('CervicalCancerTemplate-2016-02-03.xls','Interventions');
%[~,~,DWeights]=xlsread('CervicalCancerTemplate-2016-02-03.xls','DW');

%% Intervention name
for k = 1:4
    INTVN{k,1}=INTERVENTIONS{2*(k+2),2};
    INTVN{k,2}=INTERVENTIONS{2*(k+2),5};
    DW{k,1}=DWeights{k,2};
    DW{k,2}=DWeights{k,3};
end

%pal care
DW{6,2}=DWeights{7,3};
DW{7,2}=DWeights{8,3};

%pal care
DW{6,1}=DW{4,1};
DW{7,1}=DW{4,1};

%% reading data

% population data
[~,~,RAW_D]=xlsread('res4.xlsx');
% cancer data
[~,Countries] = xlsread('DataInci.xlsx','B3:Q3');
[StageDistMatrix] = xlsread('DataInci.xlsx','B32:Q36'); 
[prevalenceActualMatrix] = xlsread('DataInci.xlsx','B4:Q13');
[distHPV] = xlsread('DataInci.xlsx','B42:Q43');
[distHigh] = xlsread('DataInci.xlsx','B48:Q49');
[cancerMortalityArrayMatrix] = xlsread('DataInci.xlsx','B18:Q27');

%% defining required matrices

AgeHPVprev =[
[16	25	35	45	55	65]
[24	34	44	54	64	70]];
StageDistMatrix=StageDistMatrix';
prevalenceActualMatrix = prevalenceActualMatrix'/1000;               
distHPV=distHPV';
distHigh=distHigh';
AgeArrayL=[1 15 40 45 50  55	60	65 70 75];
AgeArrayU=[14 39 44  49 54	59	64	69 74 100];
cancerMortalityArrayMatrix = cancerMortalityArrayMatrix'/1000;

%% diagnostic rates and onset rates calculations starts here
% for I=2 for AFRE, 
%      =12 for SEARB, 
%      =16 for Peru
% for I=1:length(Countries)
for I = 2
    
    %%
    
    % reading country specific data and defining different parameters
    Countries{I}
    ind=find(strcmp(RAW_D(:,1),Countries{I}));
    DATA=RAW_D(ind,6);
    pA=cell2mat(DATA(163:243))'; %Population of female by age
    pB=cell2mat(DATA(82:162))'; %%Population of male by age
    MxW=cell2mat(DATA(406:486))'; %Natrual mortality rate by age for female
    MxM=cell2mat(DATA(325:405))'; %Natrual mortality rate by age for male
    Births=cell2mat(DATA(487:487)); %Birth population
    xa=101; pA(end+1:101)=0; pB(end+1:101)=0; MxW(end+1:101)=0; MxM(end+1:101)=0;

%% Defining country specific input parameters for MarkovChainCountry function

    PAR{1}=pA/sum(pA);
    PAR{2}=MxW;
    PAR{3}=Births/sum(pA);
    PAR{4}=prevalenceActualMatrix(I,:);
    PAR{5}=cancerMortalityArrayMatrix(I,:);
    PAR{6}= StageDistMatrix(I,:);
    PAR{7}=AgeArrayL(1,:);
    PAR{8}=AgeArrayU(1,:);
    PAR{9}=PAR{4};
    PAR{10} =distHPV(I,:);
    PAR{11} = AgeHPVprev;
    PAR{12} = distHigh(I,:);
    PAR{13}=I;
    PAR{14}=Countries(I);
    PAR{15}=pB/sum(pB);
    PAR{16}=MxM;
    
%% Parameters are now fed to the function which calculates the onset rate and diagnostic rates
    [RES,RES2]=MarkovChainCountry2(PAR);
    
%% Writing the output 

    [C_TEMPLATE]=WriteTemplateW(RES,TEMPLATEW);
    xlswrite('Temp Output1.xlsx',C_TEMPLATE);

    % some changes
    Output2=zeros(749,52);
    Output2(:)=NaN;
    Output2=num2cell(Output2);
    [~,~,Output1]=xlsread('Temp Output1.xlsx');
    Output2(1:133,1:52)=Output1(1:133,1:52);
    Output2(134:155,1:52)=Output1(706:727,1:52);
    Output2(156:243,1:52)=Output1(882:969,1:52);
    Output2(244:309,1:52)=Output1(640:705,1:52);
    Output2(310:330,1:52)=Output1(970:990,1:52);
    Output2(331,1:52)=Output1(969,1:52);
    Output2(332:353,1:52)=Output1(728:749,1:52);
    Output2(354:749,1:52)=Output1(134:529,1:52);
    
%% we now have process model output, making some chnages for producing the input file for Spectrum     
    
%% debug

%     res2_1 = [0,0,0,0;0,0,0,0;0,0,0,0;0.00379668826936542,0.00569503240404813,0.0328001679989529,0.000719548389275260;0.0769876455836706,0.115481468375506,0.653519177660332,0.0697749335292599;0.0763405304385020,0.114510795657753,0.637068475048247,0.165617127010724;0.0680137000861607,0.102020550129241,0.548498342780633,0.222838028043453;0.0690794264788688,0.103619139718303,0.484166502684547,0.234140471247613;0.0671661297701167,0.100749194655175,0.425138026592261,0.232857170121473;0.0509809131422551,0.0764713697133827,0.313871894907452,0.200205746444364;0.0436923946812650,0.0655385920218975,0.250779356700323,0.163866037386049;0.0397404475714879,0.0596106713572319,0.209534978170043,0.137261663345098];
%     res2_2 = [0,2.22044604925031e-14;0,2.22044604925031e-14;0,2.22044604925031e-14;0.697662841063462,2.22044604925031e-14;0.662790086709170,0.994284061144148;0.292560414016658,0.611144545382092;0.0374924145663416,0.00106331011631988;0.0457354918367405,0.0845633958956593;0.00348667938480050,4.02650002800210e-14;0.00544995338214347,0.00854942549700973;0.00354254835314299,0.00637834498550333;0.00428590593818790,0.00696691188822496];
%     res2_3 = [0,0,0;0,0,0;0,0,0;0,0,0;0.235765502019711,0.353648253029567,1.82408102052488;0.283960459461811,0.425940689192716,1.82408642690535;0.000750198966458661,0.00112529844968799,0.00358459709560631;0.0445948333506551,0.0668922500259827,0.206064115113882;1.48106640361751e-14,2.22159960542626e-14,6.96879089092701e-14;0.00263004438152286,0.00394506657228428,0.0108700728152445;0.00200926238895766,0.00301389358343649,0.00718568137973248;0.00237420937793574,0.00356131406690361,0.00722395452349790];
%     res2_4 = [0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0.357234012000000,0.522377791000000,0.119867991000000,0.000520165000000000,4.13429000000000e-08,6.01843000000000e-14,1.60468000000000e-21,3.86399000000000e-36;0,0,0,0,0.103761012000000,0.412439816000000,0.468759712000000,0.0150306330000000,8.82727000000000e-06,9.49505000000000e-11,1.87064000000000e-17,9.04736000000000e-31;0,0,0,0,0.00184802700000000,0.0363836160000000,0.534920730000000,0.420783550000000,0.00606247800000000,1.59979000000000e-06,7.73213000000000e-12,4.54406000000000e-23;0,0,0,0,2.43644000000000e-07,3.54440000000000e-05,0.0127840410000000,0.549054445000000,0.431901524000000,0.00622266100000000,1.64206000000000e-06,3.89315000000000e-15;0,0,0,0,5.73209000000000e-13,6.16153000000000e-10,5.45201000000000e-06,0.0127844490000000,0.549071944000000,0.431915289000000,0.00622285900000000,5.95208000000000e-09;0,0,0,0,2.48495000000000e-20,1.97370000000000e-16,4.28442000000000e-11,5.48523000000000e-06,0.0128623460000000,0.552417496000000,0.434546993000000,0.000167680000000000;0,0,0,0,3.01117000000000e-29,1.76721000000000e-24,9.41114000000000e-18,6.57844000000000e-11,8.42221000000000e-06,0.0197492840000000,0.848200665000000,0.132041629000000;0,0,0,0,4.34976000000000e-46,5.12746000000000e-40,3.31794000000000e-31,9.35658000000000e-22,4.83267000000000e-14,4.57171000000000e-08,0.000792124000000000,0.999207830000000];
%     res2_5 = [0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0.538135176000000,0.368009879000000,0.0925833930000000,0.00127138900000000,1.63305000000000e-07,3.84188000000000e-13,1.65542000000000e-20,1.30646000000000e-29,1.60203000000000e-46;0,0,0,0.237150835000000,0.440846338000000,0.301477798000000,0.0205055660000000,1.94618000000000e-05,3.38310000000000e-10,1.07713000000000e-16,6.28125000000000e-25,1.54704000000000e-40;0,0,0,0.0143939950000000,0.132530273000000,0.448904914000000,0.394970506000000,0.00919639100000000,3.92186000000000e-06,3.06330000000000e-11,4.38236000000000e-18,1.31153000000000e-31;0,0,0,1.15825000000000e-05,0.000787995000000000,0.0197220380000000,0.425700768000000,0.541171740000000,0.0126005020000000,5.37357000000000e-06,4.19720000000000e-11,5.06753000000000e-22;0,0,0,1.26404000000000e-10,6.35437000000000e-08,1.17514000000000e-05,0.00622278600000000,0.431910189000000,0.549065461000000,0.0127842980000000,5.45195000000000e-06,2.65556000000000e-14;0,0,0,2.49035000000000e-17,9.25042000000000e-14,1.26406000000000e-10,1.64212000000000e-06,0.00622288300000000,0.431916927000000,0.549074026000000,0.0127844970000000,2.51220000000000e-08;0,0,0,9.09856000000000e-26,2.49725000000000e-21,2.52149000000000e-17,8.03595000000000e-12,1.66265000000000e-06,0.00630069200000000,0.437317466000000,0.555939457000000,0.000440722000000000;0,0,0,1.69731000000000e-41,9.35694000000000e-36,1.89763000000000e-30,7.34862000000000e-23,6.13390000000000e-15,9.37757000000000e-09,0.000262583000000000,0.134667680000000,0.865069728000000];
%     RES2 = [{res2_1},{res2_2},{res2_3},{res2_4},{res2_5}];

    %% Female sheet
    
    % reading the association file
    if (I == 16 || I == 12)
        female_sheet = 3;
    else
        female_sheet = 1;
    end
        
    [~,~,data_asso] = xlsread('CervicalCancer_Asso',female_sheet);
    
    % making Output2's format exact as asso_file
    Output2(90:770,:) = Output2(68:748,:);
    Output2(772:951,1:80) = data_asso(772:951,1:80);
    Output2(68:88,1:10) = data_asso(68:88,1:10);
    
    % writing the file
    formatSpec = '%s-ProcessModOut-CervicalCancer-%s.xlsx';
    str = sprintf(formatSpec,datetime('today','Format','yyyy-MM-dd'),Countries{I});
    xlswrite(str,Output2,'Female');
    
    %% Male sheet
    
    % writing output for men
    [D_TEMPLATE2] = WriteTemplateM(RES2,TEMPLATEM);
    xlswrite(str,D_TEMPLATE2,'Male');
    
    %% Transmission sheet
    
    % for I=2 for AFRE, 
%      =12 for SEARB, 
%      =16 for Peru
% for I=1:length(Countries)
    
    if I == 16
        
        % reading transmission sheet from association file 
        [~,~,data_asso] = xlsread('CervicalCancer_Asso',5,'A1:M44');
        
        % creat a matrix of same size as data asso
        Transmission = num2cell(zeros(size(data_asso)));
        RES2{2} = num2cell(RES2{2}); RES2{4} = num2cell(RES2{4});  RES2{5} = num2cell(RES2{5});
        Transmission(:,1) = data_asso(:,1);
        Transmission(1,1) = cellstr('PERU'); Transmission(1,2:end) = num2cell(NaN);
        % results from RES2{3,4,5} are for transmision sheet. Therefore
        % arranging them as template (i.e. data asso)
        % 1. trunover rate
        Transmission(4:15,2:3) = RES2{2};Transmission(4:15,4:end) = num2cell(NaN);
        Transmission(2:3,:) = data_asso(2:3,:);
        % 2. Mixing matrix female
        Transmission(19:30,2:end) = RES2{4};
        Transmission(16:18,:) = data_asso(16:18,:);
        % 3. MIxing matrix male
        Transmission(33:44,2:end) = RES2{5};
        Transmission(31:32,:) = data_asso(31:32,:);        
        
    elseif I == 12
        
        % reading transmission sheet from association file 
        [~,~,data_asso] = xlsread('CervicalCancer_Asso',5);
        
        % keeping transmission sheet as it is in association file
        Transmission = data_asso(46:end, 1:end);
        
    elseif I == 2
        
        % reading transmission sheet from association file 
        [~,~,data_asso] = xlsread('CervicalCancer_Asso',5);
        
        % keeping transmission sheet as it is in association file
        Transmission = data_asso(1:45, 1:end);
       
    end
    
    % wrting the sheet for transmission
    xlswrite(str,Transmission,'Transmission');
    
%     e = actxserver('Excel.Application'); % # open Activex server
%     ewb = e.Workbooks.Open('C:\Users\Vijeta Deshpande\Documents\Research_all files\Cancer\Vijeta\Cervical Cancer Codes-12th Dec 2017\2018-02-12-ProcessModOut-CervicalCancer-PERU.xlsx'); % # open file (enter full path!)
%     blanksheet = ewb.Worksheets.Item(1);
%     sheet3 = ewb.Worksheets.Item(3);
%     sheet2 = ewb.Worksheets.Item(2);
%     blanksheet.Delete;
%     ewb.Save;
%     ewb.Close;
%     Excel.ActiveWorkbook.Save;
%     Excel.Quit;
%     e.Quit;

     
    'Output for Spectrum written done'
    
end


function [lTEMPLATE]=WriteTemplateW(RES,TEMPLATE)

    lTEMPLATE=TEMPLATE;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Data from Markov and Cohort model

    PrevS=RES{1};%preclinical states prevalence
    PrevC=RES{2};%clinical states prevalence
    Incidence=RES{3};
    OnsetRates=RES{4}';
    OnsetRates(end)=OnsetRates(end-1);
    DiagnosticRate=RES{5};
    DwellRate=1./RES{6};
    MortRate= RES{7};%Relative Mortality without treatment
    MortRateDF= RES{8};% Disease free mortality
    Impact = RES{9};%Relative Mortality with treatment
    TxCoverage = RES{10};% Current coverage of treatment
    PrevG=RES{11};
    HPVonsetG=RES{12};
    PrevS(1:4,:)= 0;
    PrevC(1:4,:) =0;
    OnsetRates(1:4)= 0;
    %Impact(:,1:6)=0;
    Immunity=RES{13};

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Matrix with index to cell position in TEMPLATE;

    onset=1;
    %Pre-clinical records
    PC1_Prev0=2;PC1_CFR=3;PC1_Screen=4;PC1_Progr=5;
    PC2_Prev0=6;PC2_CFR=7;PC2_Screen=8;PC2_Progr=9;
    PC3_Prev0=10;PC3_CFR=11;PC3_Screen=12;PC3_Progr=13;
    PC4_Prev0=14;PC4_CFR=15;PC4_Screen=16;PC4_Progr=17;
    PC5_Prev0=18;PC5_CFR=19;PC5_Screen=20;

    %Clinical records
    C1_Prev0=21;C1_CFR=22;C1_Progr=23;
    C2_Prev0=24;C2_CFR=25;C2_Progr=26;
    C3_Prev0=27;C3_CFR=28;C3_Progr=29;
    C4_Prev0=30;C4_CFR=31;C4_Progr=32;
    C5_Prev0=33;C5_CFR=34;

    %Treatment records
    I1_CFR=35;
    I2_CFR=36;
    I3_CFR=37;
    I4_CFR=38;
    I5_CFR=39;

    %Palliative + Treatment records
    I6_CFR_Function=40;
    I7_CFR_Function=41;

    %Palliative records
    I8_Function=42;
    I9_Function=43;


    %Screening
    I1_Screen=44;
    I2_Screen=45;
    I3_Screen=46;

    %Disability weight records
    D_C1=47;
    D_C2=48;
    D_C3=49;
    D_C4=50;
    D_C5=51;

    % Treatment coverage records
    Tx_C1 = 52;
    Tx_C2 = 53;
    Tx_C3 = 54;
    Tx_C4 = 55;
    Tx_C5 = 56;

    Tx_C6 = 57;
    Tx_C7 = 58;
    Tx_C8 = 59;
    Tx_C9 = 60;
    Tx_C10 = 60;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Disease Free
    DSF1_incH=95;DSF1_incL=96;
    DSF1_inc1618=97;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Pre-cancer Stages
    %HPV 16/18
    PCH1_1618_Prev0=61;PCH1_1618_CFR=62;PCH1_1618_Rem=63;PCH1_1618_Progr=64;
    PCH2_1618_Prev0=65;PCH2_1618_CFR=66;PCH2_1618_Rem1=67;PCH2_1618_Rem2=68;PCH2_1618_Progr=69;
    PCH3_1618_Prev0=70;PCH3_1618_CFR=71;PCH3_1618_Rem1=72;PCH3_1618_Rem2=73;PCH3_1618_Progr=74;

    %HPV High Risk
    PCH1_H_Prev0=104;PCH1_H_CFR=105;PCH1_H_Rem=106;PCH1_H_Progr=107;
    PCH2_H_Prev0=108;PCH2_H_CFR=109;PCH2_H_Rem1=110;PCH2_H_Rem2=111;PCH2_H_Progr=112;
    PCH3_H_Prev0=113;PCH3_H_CFR=114;PCH3_H_Rem1=115;PCH3_H_Rem2=116;PCH3_H_Progr=117;

    %HPV Low Risk
    PCL1_Prev0=75;PCL1_CFR=76;PCL1_Rem=77;PCL1_Progr=78;
    PCL2_Prev0=79;PCL2_CFR=80;PCL2_Rem1=81;PCL2_Rem2=82;PCL2_Progr=83;
    PCL3_Prev0=84;PCL3_CFR=85;PCL3_Rem1=86;PCL3_Rem2=87;PCL3_Progr=88;

    %HPV Regression
    PCR1_Prev0=89;PCR1_CFR=90;PCR1_Rem=91;
    PCR2_Prev0=98;PCR2_CFR=99;PCR2_Rem=100;
    PCR3_Prev0=101;PCR3_CFR=102;PCR3_Rem=103;
    %Immunity
    PCM1_Prev0=92;PCM1_CFR=93;PCM1_Progr=94;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %1=row 1, 2=row end, 3=row col
    %pre-cancer stages
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %HPV 16/18
    R{PCH1_1618_Prev0,1}=75;R{PCH1_1618_Prev0,2}=86;R{PCH1_1618_Prev0,3}=4;
    R{PCH1_1618_CFR,1}=75;R{PCH1_1618_CFR,2}=86;R{PCH1_1618_CFR,3}=6;
    R{PCH1_1618_Rem,1}=75;R{PCH1_1618_Rem,2}=86;R{PCH1_1618_Rem,3}=9;
    R{PCH1_1618_Progr,1}=75;R{PCH1_1618_Progr,2}=86;R{PCH1_1618_Progr,3}=12;

    %CIN1 16/18
    R{PCH2_1618_Prev0,1}=97;R{PCH2_1618_Prev0,2}=108;R{PCH2_1618_Prev0,3}=4;
    R{PCH2_1618_CFR,1}=97;R{PCH2_1618_CFR,2}=108;R{PCH2_1618_CFR,3}=6;
    R{PCH2_1618_Rem1,1}=97;R{PCH2_1618_Rem1,2}=108;R{PCH2_1618_Rem1,3}=9;
    R{PCH2_1618_Rem2,1}=97;R{PCH2_1618_Rem2,2}=108;R{PCH2_1618_Rem2,3}=12;
    R{PCH2_1618_Progr,1}=97;R{PCH2_1618_Progr,2}=108;R{PCH2_1618_Progr,3}=15;

    %CIN2/3 16/18
    R{PCH3_1618_Prev0,1}=119;R{PCH3_1618_Prev0,2}=130;R{PCH3_1618_Prev0,3}=4;
    R{PCH3_1618_CFR,1}=119;R{PCH3_1618_CFR,2}=130;R{PCH3_1618_CFR,3}=6;
    R{PCH3_1618_Rem1,1}=119;R{PCH3_1618_Rem1,2}=130;R{PCH3_1618_Rem1,3}=9;
    R{PCH3_1618_Rem2,1}=119;R{PCH3_1618_Rem2,2}=130;R{PCH3_1618_Rem2,3}=12;
    R{PCH3_1618_Progr,1}=119;R{PCH3_1618_Progr,2}=130;R{PCH3_1618_Progr,3}=15;

    %HPV High Risk
    R{PCH1_H_Prev0,1}=889;R{PCH1_H_Prev0,2}=900;R{PCH1_H_Prev0,3}=4;
    R{PCH1_H_CFR,1}=889;R{PCH1_H_CFR,2}=900;R{PCH1_H_CFR,3}=6;
    R{PCH1_H_Rem,1}=889;R{PCH1_H_Rem,2}=900;R{PCH1_H_Rem,3}=9;
    R{PCH1_H_Progr,1}=889;R{PCH1_H_Progr,2}=900;R{PCH1_H_Progr,3}=12;

    %CIN1 High Risk
    R{PCH2_H_Prev0,1}=911;R{PCH2_H_Prev0,2}=922;R{PCH2_H_Prev0,3}=4;
    R{PCH2_H_CFR,1}=911;R{PCH2_H_CFR,2}=922;R{PCH2_H_CFR,3}=6;
    R{PCH2_H_Rem1,1}=911;R{PCH2_H_Rem1,2}=922;R{PCH2_H_Rem1,3}=9;
    R{PCH2_H_Rem2,1}=911;R{PCH2_H_Rem2,2}=922;R{PCH2_H_Rem2,3}=12;
    R{PCH2_H_Progr,1}=911;R{PCH2_H_Progr,2}=922;R{PCH2_H_Progr,3}=15;

    %CIN2/3 High Risk
    R{PCH3_H_Prev0,1}=933;R{PCH3_H_Prev0,2}=944;R{PCH3_H_Prev0,3}=4;
    R{PCH3_H_CFR,1}=933;R{PCH3_H_CFR,2}=944;R{PCH3_H_CFR,3}=6;
    R{PCH3_H_Rem1,1}=933;R{PCH3_H_Rem1,2}=944;R{PCH3_H_Rem1,3}=9;
    R{PCH3_H_Rem2,1}=933;R{PCH3_H_Rem2,2}=944;R{PCH3_H_Rem2,3}=12;
    R{PCH3_H_Progr,1}=933;R{PCH3_H_Progr,2}=944;R{PCH3_H_Progr,3}=15;

    %HPV Low Risk
    R{PCL1_Prev0,1}=647;R{PCL1_Prev0,2}=658;R{PCL1_Prev0,3}=4;
    R{PCL1_CFR,1}=647;R{PCL1_CFR,2}=658;R{PCL1_CFR,3}=6;
    R{PCL1_Rem,1}=647;R{PCL1_Rem,2}=658;R{PCL1_Rem,3}=9;
    R{PCL1_Progr,1}=647;R{PCL1_Progr,2}=658;R{PCL1_Progr,3}=12;

    %CIN1 Low Risk
    R{PCL2_Prev0,1}=669;R{PCL2_Prev0,2}=680;R{PCL2_Prev0,3}=4;
    R{PCL2_CFR,1}=669;R{PCL2_CFR,2}=680;R{PCL2_CFR,3}=6;
    R{PCL2_Rem1,1}=669;R{PCL2_Rem1,2}=680;R{PCL2_Rem1,3}=9;
    R{PCL2_Rem2,1}=669;R{PCL2_Rem2,2}=680;R{PCL2_Rem2,3}=12;
    R{PCL2_Progr,1}=669;R{PCL2_Progr,2}=680;R{PCL2_Progr,3}=15;

    %CIN2/3 Low Risk
    R{PCL3_Prev0,1}=691;R{PCL3_Prev0,2}=702;R{PCL3_Prev0,3}=4;
    R{PCL3_CFR,1}=691;R{PCL3_CFR,2}=702;R{PCL3_CFR,3}=6;
    R{PCL3_Rem1,1}=691;R{PCL3_Rem1,2}=702;R{PCL3_Rem1,3}=9;
    R{PCL3_Rem2,1}=691;R{PCL3_Rem2,2}=702;R{PCL3_Rem2,3}=12;
    R{PCL3_Progr,1}=691;R{PCL3_Progr,2}=702;R{PCL3_Progr,3}=15;

    %HPV Regression
    R{PCR1_Prev0,1}=713;R{PCR1_Prev0,2}=724;R{PCR1_Prev0,3}=4;
    R{PCR1_CFR,1}=713;R{PCR1_CFR,2}=724;R{PCR1_CFR,3}=6;
    R{PCR1_Rem,1}=713;R{PCR1_Rem,2}=724;R{PCR1_Rem,3}=9;
    R{PCR2_Prev0,1}=955;R{PCR2_Prev0,2}=966;R{PCR2_Prev0,3}=4;
    R{PCR2_CFR,1}=955;R{PCR2_CFR,2}=966;R{PCR2_CFR,3}=6;
    R{PCR2_Rem,1}=955;R{PCR2_Rem,2}=966;R{PCR2_Rem,3}=9;
    R{PCR3_Prev0,1}=977;R{PCR3_Prev0,2}=988;R{PCR3_Prev0,3}=4;
    R{PCR3_CFR,1}=977;R{PCR3_CFR,2}=988;R{PCR3_CFR,3}=6;
    R{PCR3_Rem,1}=977;R{PCR3_Rem,2}=988;R{PCR3_Rem,3}=9;




    %Immunity
    R{PCM1_Prev0,1}=735;R{PCM1_Prev0,2}=746;R{PCM1_Prev0,3}=4;
    R{PCM1_CFR,1}=735;R{PCM1_CFR,2}=746;R{PCM1_CFR,3}=6;
    R{PCM1_Progr,1}=735;R{PCM1_Progr,2}=746;R{PCM1_Progr,3}=9;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




    %cancer stages
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %disease free
    R{onset,1}=9;R{onset,2}=20;R{onset,3}=7;
    R{DSF1_inc1618,1}=9;R{DSF1_inc1618,2}=20;R{DSF1_inc1618,3}=7;
    R{DSF1_incH,1}=9;R{DSF1_incH,2}=20;R{DSF1_incH,3}=10;
    R{DSF1_incL,1}=9;R{DSF1_incL,2}=20;R{DSF1_incL,3}=13;
    %pre clinical
    %stage 0
    R{PC1_Prev0,1}=141;R{PC1_Prev0,2}=152;R{PC1_Prev0,3}=4;
    R{PC1_CFR,1}=141;R{PC1_CFR,2}=152;R{PC1_CFR,3}=6;
    R{PC1_Screen,1}=141;R{PC1_Screen,2}=152;R{PC1_Screen,3}=9;
    R{PC1_Progr,1}=141;R{PC1_Progr,2}=152;R{PC1_Progr,3}=12;

    %stage 1
    R{PC2_Prev0,1}=163;R{PC2_Prev0,2}=174;R{PC2_Prev0,3}=4;
    R{PC2_CFR,1}=163;R{PC2_CFR,2}=174;R{PC2_CFR,3}=6;
    R{PC2_Screen,1}=163;R{PC2_Screen,2}=174;R{PC2_Screen,3}=9;
    R{PC2_Progr,1}=163;R{PC2_Progr,2}=174;R{PC2_Progr,3}=12;

    %stage 2
    R{PC3_Prev0,1}=185;R{PC3_Prev0,2}=196;R{PC3_Prev0,3}=4;
    R{PC3_CFR,1}=185;R{PC3_CFR,2}=196;R{PC3_CFR,3}=6;
    R{PC3_Screen,1}=185;R{PC3_Screen,2}=196;R{PC3_Screen,3}=9;
    R{PC3_Progr,1}=185;R{PC3_Progr,2}=196;R{PC3_Progr,3}=12;

    %stage 3
    R{PC4_Prev0,1}=207;R{PC4_Prev0,2}=218;R{PC4_Prev0,3}=4;
    R{PC4_CFR,1}=207;R{PC4_CFR,2}=218;R{PC4_CFR,3}=6;
    R{PC4_Screen,1}=207;R{PC4_Screen,2}=218;R{PC4_Screen,3}=9;
    R{PC4_Progr,1}=207;R{PC4_Progr,2}=218;R{PC4_Progr,3}=12;  

    %stage 4
    R{PC5_Prev0,1}=229;R{PC5_Prev0,2}=240;R{PC5_Prev0,3}=4;
    R{PC5_CFR,1}=229;R{PC5_CFR,2}=240;R{PC5_CFR,3}=6;
    R{PC5_Screen,1}=229;R{PC5_Screen,2}=240;R{PC5_Screen,3}=9;
    %R{PC4_Progr,1}=207;R{PC4_Progr,2}=218;R{PC4_Progr,3}=12;  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %clinical 0
    R{C1_Prev0,1}=251;R{C1_Prev0,2}=262;R{C1_Prev0,3}=4;
    R{C1_CFR,1}=251;R{C1_CFR,2}=262;R{C1_CFR,3}=6;
    R{C1_Progr,1}=251;R{C1_Progr,2}=262;R{C1_Progr,3}=9;

    %clinical 1
    R{C2_Prev0,1}=273;R{C2_Prev0,2}=284;R{C2_Prev0,3}=4;
    R{C2_CFR,1}=273;R{C2_CFR,2}=284;R{C2_CFR,3}=6;
    R{C2_Progr,1}=273;R{C2_Progr,2}=284;R{C2_Progr,3}=9;

    %clinical 2
    R{C3_Prev0,1}=295;R{C3_Prev0,2}=306;R{C3_Prev0,3}=4;
    R{C3_CFR,1}=295;R{C3_CFR,2}=306;R{C3_CFR,3}=6;
    R{C3_Progr,1}=295;R{C3_Progr,2}=306;R{C3_Progr,3}=9;

    %clinical 3
    R{C4_Prev0,1}=317;R{C4_Prev0,2}=328;R{C4_Prev0,3}=4;
    R{C4_CFR,1}=317;R{C4_CFR,2}=328;R{C4_CFR,3}=6;
    R{C4_Progr,1}=317;R{C4_Progr,2}=328;R{C4_Progr,3}=9;

    %clinical 4
    R{C5_Prev0,1}=339;R{C5_Prev0,2}=350;R{C5_Prev0,3}=4;
    R{C5_CFR,1}=339;R{C5_CFR,2}=350;R{C5_CFR,3}=6;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Interventions 7-CFR, 14 Function-DW
    R{I1_CFR,1}=383;R{I1_CFR,2}=383;R{I1_CFR,3}=7;R{I1_CFR,4}=14;
    R{I2_CFR,1}=405;R{I2_CFR,2}=405;R{I2_CFR,3}=7;R{I2_CFR,4}=14;
    R{I3_CFR,1}=427;R{I3_CFR,2}=427;R{I3_CFR,3}=7;R{I3_CFR,4}=14;
    R{I4_CFR,1}=449;R{I4_CFR,2}=449;R{I4_CFR,3}=7;R{I4_CFR,4}=14;
    R{I5_CFR,1}=471;R{I5_CFR,2}=471;R{I5_CFR,3}=7;R{I5_CFR,4}=14;

    R{I6_CFR_Function,1}=493;R{I6_CFR_Function,2}=493;R{I6_CFR_Function,3}=7;R{I6_CFR_Function,4}=14;
    R{I7_CFR_Function,1}=515;R{I7_CFR_Function,2}=515;R{I7_CFR_Function,3}=7;R{I7_CFR_Function,4}=14;

    R{I8_Function,1}=537;R{I8_Function,2}=537;R{I8_Function,4}=7;
    R{I9_Function,1}=559;R{I9_Function,2}=559;R{I9_Function,4}=7;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Screening Interventions
    %bi-annual screening
    R{I1_Screen,1}=581;R{I1_Screen,2}=581;R{I1_Screen,3}=7;R{I1_Screen,4}=14;R{I1_Screen,5}=21;R{I1_Screen,6}=28;
    R{I2_Screen,1}=603;R{I2_Screen,2}=603;R{I2_Screen,3}=7;R{I2_Screen,4}=14;R{I2_Screen,5}=21;R{I2_Screen,6}=28;
    %clinical examination
    R{I3_Screen,1}=625;R{I3_Screen,2}=625;R{I3_Screen,3}=7;R{I3_Screen,4}=14;R{I3_Screen,5}=21;R{I3_Screen,6}=28;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Disability weights
    %using row indices only
    R{D_C1,1}=137;R{D_C1,2}=247;
    R{D_C2,1}=159;R{D_C2,2}=269;
    R{D_C3,1}=181;R{D_C3,2}=291;
    R{D_C4,1}=203;R{D_C4,2}=313;
    R{D_C5,1}=225;R{D_C5,2}=335;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Current and scale-up coverage
    R{Tx_C1,1}=383;R{Tx_C1,2}=383;R{Tx_C1,3}=383;R{Tx_C1,4}=8;R{Tx_C1,5}=9;R{Tx_C1,6}=15;R{Tx_C1,7}=16;
    R{Tx_C2,1}=405;R{Tx_C2,2}=405;R{Tx_C2,3}=405;R{Tx_C2,4}=8;R{Tx_C2,5}=9;R{Tx_C2,6}=15;R{Tx_C2,7}=16;
    R{Tx_C3,1}=427;R{Tx_C3,2}=427;R{Tx_C3,3}=427;R{Tx_C3,4}=8;R{Tx_C3,5}=9;R{Tx_C3,6}=15;R{Tx_C3,7}=16;
    R{Tx_C4,1}=449;R{Tx_C4,2}=449;R{Tx_C4,3}=449;R{Tx_C4,4}=8;R{Tx_C4,5}=9;R{Tx_C4,6}=15;R{Tx_C4,7}=16;
    R{Tx_C5,1}=471;R{Tx_C5,2}=471;R{Tx_C5,3}=471;R{Tx_C5,4}=8;R{Tx_C5,5}=9;R{Tx_C5,6}=15;R{Tx_C5,7}=16;


    R{Tx_C6,1}=493;R{Tx_C6,2}=493;R{Tx_C6,3}=493;R{Tx_C6,4}=8;R{Tx_C6,5}=9;R{Tx_C6,6}=15;R{Tx_C6,7}=16;
    R{Tx_C7,1}=515;R{Tx_C7,2}=515;R{Tx_C7,3}=515;R{Tx_C7,4}=8;R{Tx_C7,5}=9;R{Tx_C7,6}=15;R{Tx_C7,7}=16;
    R{Tx_C8,1}=537;R{Tx_C8,2}=537;R{Tx_C8,3}=537;R{Tx_C8,4}=8;R{Tx_C8,5}=9;R{Tx_C8,6}=15;R{Tx_C8,7}=16;
    R{Tx_C9,1}=559;R{Tx_C9,2}=559;R{Tx_C9,3}=559;R{Tx_C9,4}=8;R{Tx_C9,5}=9;R{Tx_C9,6}=15;R{Tx_C9,7}=16;


    for I_loop=1:12

        %onset to first stage
        lTEMPLATE{ R{onset,1}+I_loop-1 ,R{onset,3}}=OnsetRates(I_loop,1);    
        lTEMPLATE{ R{DSF1_inc1618,1}+I_loop-1 ,R{DSF1_inc1618,3} }=HPVonsetG(I_loop,1);
        lTEMPLATE{ R{DSF1_incH,1}+I_loop-1 ,R{DSF1_incH,3} }=HPVonsetG(I_loop,2);
        lTEMPLATE{ R{DSF1_incL,1}+I_loop-1 ,R{DSF1_incL,3} }=HPVonsetG(I_loop,3);


        %progression
        lTEMPLATE{ R{PC1_Prev0,1}+I_loop-1 ,R{PC1_Prev0,3} }=PrevS(I_loop,1);
        lTEMPLATE{ R{PC2_Prev0,1}+I_loop-1 ,R{PC2_Prev0,3} }=PrevS(I_loop,2);
        lTEMPLATE{ R{PC3_Prev0,1}+I_loop-1 ,R{PC3_Prev0,3} }=PrevS(I_loop,3);
        lTEMPLATE{ R{PC4_Prev0,1}+I_loop-1 ,R{PC4_Prev0,3} }=PrevS(I_loop,4);
        lTEMPLATE{ R{PC5_Prev0,1}+I_loop-1 ,R{PC5_Prev0,3} }=PrevS(I_loop,5);

        %progression
        lTEMPLATE{ R{PC1_Progr,1}+I_loop-1 ,R{PC1_Progr,3} }=DwellRate(I_loop,1);
        lTEMPLATE{ R{PC2_Progr,1}+I_loop-1 ,R{PC2_Progr,3} }=DwellRate(I_loop,2);
        lTEMPLATE{ R{PC3_Progr,1}+I_loop-1 ,R{PC3_Progr,3} }=DwellRate(I_loop,3);
        lTEMPLATE{ R{PC4_Progr,1}+I_loop-1 ,R{PC4_Progr,3} }=DwellRate(I_loop,4);

        %screening process
        lTEMPLATE{ R{PC1_Screen,1}+I_loop-1 ,R{PC1_Screen,3} }=DiagnosticRate(I_loop,1);
        lTEMPLATE{ R{PC2_Screen,1}+I_loop-1 ,R{PC2_Screen,3} }=DiagnosticRate(I_loop,2);
        lTEMPLATE{ R{PC3_Screen,1}+I_loop-1 ,R{PC3_Screen,3} }=DiagnosticRate(I_loop,3);
        lTEMPLATE{ R{PC4_Screen,1}+I_loop-1 ,R{PC4_Screen,3} }=DiagnosticRate(I_loop,4);
        lTEMPLATE{ R{PC5_Screen,1}+I_loop-1 ,R{PC5_Screen,3} }=DiagnosticRate(I_loop,5);

        %Cancer Mortality pre-Clinical
        lTEMPLATE{ R{PC1_CFR,1}+I_loop-1 ,R{PC1_CFR,3} }=0;%MortRateDF(I,1);
        lTEMPLATE{ R{PC2_CFR,1}+I_loop-1 ,R{PC2_CFR,3} }=0;%MortRateDF(I,1);
        lTEMPLATE{ R{PC3_CFR,1}+I_loop-1 ,R{PC3_CFR,3} }=0;%MortRateDF(I,1);
        lTEMPLATE{ R{PC4_CFR,1}+I_loop-1 ,R{PC4_CFR,3} }=0;%MortRate(4,I);
        lTEMPLATE{ R{PC5_CFR,1}+I_loop-1 ,R{PC5_CFR,3} }=MortRate(5,I_loop);

        %Clinical
        %Base Cancer Mortality (after considering coverage of TX among those
        %diagnosed
        lTEMPLATE{ R{C1_CFR,1}+I_loop-1 ,R{C1_CFR,3} }=MortRate(1,I_loop);
        lTEMPLATE{ R{C2_CFR,1}+I_loop-1 ,R{C2_CFR,3} }=MortRate(2,I_loop);
        lTEMPLATE{ R{C3_CFR,1}+I_loop-1 ,R{C3_CFR,3} }=MortRate(3,I_loop);
        lTEMPLATE{ R{C4_CFR,1}+I_loop-1 ,R{C4_CFR,3} }=MortRate(4,I_loop);
        lTEMPLATE{ R{C5_CFR,1}+I_loop-1 ,R{C5_CFR,3} }=MortRate(5,I_loop);

        %prev 0
        lTEMPLATE{ R{C1_Prev0,1}+I_loop-1 ,R{C1_Prev0,3} }=PrevC(I_loop,1);
        lTEMPLATE{ R{C2_Prev0,1}+I_loop-1 ,R{C2_Prev0,3} }=PrevC(I_loop,2);
        lTEMPLATE{ R{C3_Prev0,1}+I_loop-1 ,R{C3_Prev0,3} }=PrevC(I_loop,3);
        lTEMPLATE{ R{C4_Prev0,1}+I_loop-1 ,R{C4_Prev0,3} }=PrevC(I_loop,4);
        lTEMPLATE{ R{C5_Prev0,1}+I_loop-1 ,R{C5_Prev0,3} }=PrevC(I_loop,5);

        %progression
        lTEMPLATE{ R{C1_Progr,1}+I_loop-1 ,R{C1_Progr,3} }=0;%DwellRate(I,1);
        lTEMPLATE{ R{C2_Progr,1}+I_loop-1 ,R{C2_Progr,3} }=0;%DwellRate(I,2);
        lTEMPLATE{ R{C3_Progr,1}+I_loop-1 ,R{C3_Progr,3} }=0;%DwellRate(I,3);
        lTEMPLATE{ R{C4_Progr,1}+I_loop-1 ,R{C4_Progr,3} }=0;%DwellRate(I,3);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Treatment impact on MORT
        lTEMPLATE{ R{I1_CFR,1}+I_loop-1 ,R{I1_CFR,3} }=Impact(1,I_loop);%INTVN{1,2}*100;
        lTEMPLATE{ R{I2_CFR,1}+I_loop-1 ,R{I2_CFR,3} }=Impact(2,I_loop);%INTVN{2,2}*100;
        lTEMPLATE{ R{I3_CFR,1}+I_loop-1 ,R{I3_CFR,3} }=Impact(3,I_loop);%INTVN{3,2}*100;
        lTEMPLATE{ R{I4_CFR,1}+I_loop-1 ,R{I4_CFR,3} }=Impact(4,I_loop);%INTVN{4,2}*100;
        lTEMPLATE{ R{I5_CFR,1}+I_loop-1 ,R{I5_CFR,3} }=Impact(5,I_loop);%INTVN{4,2}*100;

        % Treatment+palliative care impact on Mort
        lTEMPLATE{ R{I6_CFR_Function,1}+I_loop-1 ,R{I6_CFR_Function,3} }=Impact(5,I_loop);%INTVN{4,2}*100;
        lTEMPLATE{ R{I7_CFR_Function,1}+I_loop-1 ,R{I7_CFR_Function,3} }=Impact(5,I_loop);%INTVN{4,2}*100;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Baseline treatment coverage for treatment
        lTEMPLATE{ R{Tx_C1,1}+I_loop-1 ,R{Tx_C1,4} }  = TxCoverage(1,I_loop)*100;
        lTEMPLATE{ R{Tx_C2,1}+I_loop-1 ,R{Tx_C2,4} }  = TxCoverage(2,I_loop)*100;   
        lTEMPLATE{ R{Tx_C3,1}+I_loop-1 ,R{Tx_C3,4} }  = TxCoverage(3,I_loop)*100;
        lTEMPLATE{ R{Tx_C4,1}+I_loop-1 ,R{Tx_C4,4} }  = TxCoverage(4,I_loop)*100;
        lTEMPLATE{ R{Tx_C5,1}+I_loop-1 ,R{Tx_C5,4} }  = TxCoverage(5,I_loop)*100;

        lTEMPLATE{ R{Tx_C6,1}+I_loop-1 ,R{Tx_C6,4} }  = TxCoverage(5,I_loop)*100;
        lTEMPLATE{ R{Tx_C7,1}+I_loop-1 ,R{Tx_C7,4} }  = TxCoverage(5,I_loop)*100;
        lTEMPLATE{ R{Tx_C8,1}+I_loop-1 ,R{Tx_C8,4} }  = TxCoverage(5,I_loop)*100;
        lTEMPLATE{ R{Tx_C9,1}+I_loop-1 ,R{Tx_C9,4} }  = TxCoverage(5,I_loop)*100;


        %Default intervention treatment coverage= current + 80% of those not
        lTEMPLATE{ R{Tx_C1,1}+I_loop-1 ,R{Tx_C1,5} }  = ((1-TxCoverage(1,I_loop))*.8+TxCoverage(1,I_loop)) *100;
        lTEMPLATE{ R{Tx_C2,1}+I_loop-1 ,R{Tx_C2,5} }  = ((1-TxCoverage(2,I_loop))*.8+TxCoverage(2,I_loop)) *100;   
        lTEMPLATE{ R{Tx_C3,1}+I_loop-1 ,R{Tx_C3,5} }  = ((1-TxCoverage(3,I_loop))*.8+TxCoverage(3,I_loop)) *100;
        lTEMPLATE{ R{Tx_C4,1}+I_loop-1 ,R{Tx_C4,5} }  = ((1-TxCoverage(4,I_loop))*.8+TxCoverage(4,I_loop)) *100;
        lTEMPLATE{ R{Tx_C5,1}+I_loop-1 ,R{Tx_C5,5} }  = ((1-TxCoverage(5,I_loop))*.8+TxCoverage(5,I_loop)) *100;

        lTEMPLATE{ R{Tx_C6,1}+I_loop-1 ,R{Tx_C6,5} }  = ((1-TxCoverage(5,I_loop))*.8+TxCoverage(5,I_loop)) *100;
        lTEMPLATE{ R{Tx_C7,1}+I_loop-1 ,R{Tx_C7,5} }  = ((1-TxCoverage(5,I_loop))*.8+TxCoverage(5,I_loop)) *100;
        lTEMPLATE{ R{Tx_C8,1}+I_loop-1 ,R{Tx_C8,5} }  = ((1-TxCoverage(5,I_loop))*.8+TxCoverage(5,I_loop)) *100;
        lTEMPLATE{ R{Tx_C9,1}+I_loop-1 ,R{Tx_C9,5} }  = ((1-TxCoverage(5,I_loop))*.8+TxCoverage(5,I_loop)) *100;

        % Baseline coverage for disability weight = same as Tx coverage
        lTEMPLATE{ R{Tx_C1,1}+I_loop-1 ,R{Tx_C1,6} }  = TxCoverage(1,I_loop)*100;
        lTEMPLATE{ R{Tx_C2,1}+I_loop-1 ,R{Tx_C2,6} }  = TxCoverage(2,I_loop)*100;   
        lTEMPLATE{ R{Tx_C3,1}+I_loop-1 ,R{Tx_C3,6} }  = TxCoverage(3,I_loop)*100;
        lTEMPLATE{ R{Tx_C4,1}+I_loop-1 ,R{Tx_C4,6} }  = TxCoverage(4,I_loop)*100;
        lTEMPLATE{ R{Tx_C5,1}+I_loop-1 ,R{Tx_C5,6} }  = TxCoverage(5,I_loop)*100;

        lTEMPLATE{ R{Tx_C6,1}+I_loop-1 ,R{Tx_C6,6} }  = TxCoverage(5,I_loop)*100;
        lTEMPLATE{ R{Tx_C7,1}+I_loop-1 ,R{Tx_C7,6} }  = TxCoverage(5,I_loop)*100;

        %Default intervention coverage for disability weight= current + 80% of those not
        lTEMPLATE{ R{Tx_C1,1}+I_loop-1 ,R{Tx_C1,7} }  = ((1-TxCoverage(1,I_loop))*.8+TxCoverage(1,I_loop)) *100;
        lTEMPLATE{ R{Tx_C2,1}+I_loop-1 ,R{Tx_C2,7} }  = ((1-TxCoverage(2,I_loop))*.8+TxCoverage(2,I_loop)) *100;   
        lTEMPLATE{ R{Tx_C3,1}+I_loop-1 ,R{Tx_C3,7} }  = ((1-TxCoverage(3,I_loop))*.8+TxCoverage(3,I_loop)) *100;
        lTEMPLATE{ R{Tx_C4,1}+I_loop-1 ,R{Tx_C4,7} }  = ((1-TxCoverage(4,I_loop))*.8+TxCoverage(4,I_loop)) *100;
        lTEMPLATE{ R{Tx_C5,1}+I_loop-1 ,R{Tx_C5,7} }  = ((1-TxCoverage(5,I_loop))*.8+TxCoverage(5,I_loop)) *100;

        lTEMPLATE{ R{Tx_C6,1}+I_loop-1 ,R{Tx_C6,7} }  = ((1-TxCoverage(5,I_loop))*.8+TxCoverage(5,I_loop)) *100;
        lTEMPLATE{ R{Tx_C7,1}+I_loop-1 ,R{Tx_C7,7} }  = ((1-TxCoverage(5,I_loop))*.8+TxCoverage(5,I_loop)) *100;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Treatment
        %impact on disability weights, see def of DW in start of function
        lTEMPLATE{ R{I1_CFR,1}+I_loop-1 ,R{I1_CFR,4} }=-(DW{1,1}-DW{1,2})/DW{1,1}*100;
        lTEMPLATE{ R{I2_CFR,1}+I_loop-1 ,R{I2_CFR,4} }=-(DW{2,1}-DW{2,2})/DW{2,1}*100;
        lTEMPLATE{ R{I3_CFR,1}+I_loop-1 ,R{I3_CFR,4} }=-(DW{3,1}-DW{3,2})/DW{3,1}*100;
        lTEMPLATE{ R{I4_CFR,1}+I_loop-1 ,R{I4_CFR,4} }=-(DW{4,1}-DW{4,2})/DW{4,1}*100;
        lTEMPLATE{ R{I5_CFR,1}+I_loop-1 ,R{I5_CFR,4} }=-(DW{5,1}-DW{5,2})/DW{5,1}*100;

        %Treatment+ Palliative Care
        lTEMPLATE{ R{I6_CFR_Function,1}+I_loop-1 ,R{I6_CFR_Function,4} }=-(DW{5,1}-DW{6,2})/DW{5,1}*100;
        lTEMPLATE{ R{I7_CFR_Function,1}+I_loop-1 ,R{I7_CFR_Function,4} }=-(DW{5,1}-DW{7,2})/DW{5,1}*100;

        %Palliative Care
        lTEMPLATE{ R{I8_Function,1}+I_loop-1 ,R{I8_Function,4} }=-(DW{5,1}-DW{6,2})/DW{5,1}*100;
        lTEMPLATE{ R{I9_Function,1}+I_loop-1 ,R{I9_Function,4} }=-(DW{5,1}-DW{7,2})/DW{5,1}*100;


        %Bi annual screening of 50-69, age match not perfect
        lTEMPLATE{ R{I1_Screen,1}+I_loop-1 ,R{I1_Screen,3} }=100*0/2;
        lTEMPLATE{ R{I1_Screen,1}+I_loop-1 ,R{I1_Screen,4} }=100*0/2;
        lTEMPLATE{ R{I1_Screen,1}+I_loop-1 ,R{I1_Screen,5} }=100*0/2;
        lTEMPLATE{ R{I1_Screen,1}+I_loop-1 ,R{I1_Screen,6} }=100*0/2;
        %lTEMPLATE{ R{I1_Screen,1}+I-1 ,R{I1_Screen,7} }=100*0/2;
        if(ismember(I_loop,[9,10]))
            lTEMPLATE{ R{I1_Screen,1}+I_loop-1 ,R{I1_Screen,3} }=100*1/2;
            lTEMPLATE{ R{I1_Screen,1}+I_loop-1 ,R{I1_Screen,4} }=100*1/2;
            lTEMPLATE{ R{I1_Screen,1}+I_loop-1 ,R{I1_Screen,5} }=100*1/2;
            lTEMPLATE{ R{I1_Screen,1}+I_loop-1 ,R{I1_Screen,6} }=100*1/2;
            %lTEMPLATE{ R{I1_Screen,1}+I-1 ,R{I1_Screen,7} }=100*1/2;
        end

        %Bi annual screening of 40-69, age match not perfect
        lTEMPLATE{ R{I2_Screen,1}+I_loop-1 ,R{I2_Screen,3} }=100*0/2;
        lTEMPLATE{ R{I2_Screen,1}+I_loop-1 ,R{I2_Screen,4} }=100*0/2;
        lTEMPLATE{ R{I2_Screen,1}+I_loop-1 ,R{I2_Screen,5} }=100*0/2;
        lTEMPLATE{ R{I2_Screen,1}+I_loop-1 ,R{I2_Screen,6} }=100*0/2;
        %lTEMPLATE{ R{I2_Screen,1}+I-1 ,R{I2_Screen,7} }=100*0/2;
        if(ismember(I_loop,[8,9,10]))
            lTEMPLATE{ R{I2_Screen,1}+I_loop-1 ,R{I2_Screen,3} }=100*1/2;
            lTEMPLATE{ R{I2_Screen,1}+I_loop-1 ,R{I2_Screen,4} }=100*1/2;
            lTEMPLATE{ R{I2_Screen,1}+I_loop-1 ,R{I2_Screen,5} }=100*1/2;
            lTEMPLATE{ R{I2_Screen,1}+I_loop-1 ,R{I2_Screen,6} }=100*1/2;
            %lTEMPLATE{ R{I2_Screen,1}+I-1 ,R{I2_Screen,7} }=100*1/2;
        end

        %Screening-by clinical exam
        lTEMPLATE{ R{I3_Screen,1}+I_loop-1 ,R{I3_Screen,3} }=100*0/2;
        lTEMPLATE{ R{I3_Screen,1}+I_loop-1 ,R{I3_Screen,4} }=100*0/2;
        lTEMPLATE{ R{I3_Screen,1}+I_loop-1 ,R{I3_Screen,5} }=100*0/2;
        lTEMPLATE{ R{I3_Screen,1}+I_loop-1 ,R{I3_Screen,6} }=100*0/2;
        %lTEMPLATE{ R{I3_Screen,1}+I-1 ,R{I3_Screen,7} }=100*0/2;
        if(ismember(I_loop,[8,9,10]))%assume it takes longer than bi-anual
            lTEMPLATE{ R{I3_Screen,1}+I_loop-1 ,R{I3_Screen,3} }=100*1/4;
            lTEMPLATE{ R{I3_Screen,1}+I_loop-1 ,R{I3_Screen,4} }=100*1/4;
            lTEMPLATE{ R{I3_Screen,1}+I_loop-1 ,R{I3_Screen,5} }=100*1/4;
            lTEMPLATE{ R{I3_Screen,1}+I_loop-1 ,R{I3_Screen,6} }=100*1/4;
            %lTEMPLATE{ R{I3_Screen,1}+I-1 ,R{I3_Screen,7} }=100*1/4;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %Disability weights
        %preclinical
        lTEMPLATE{ R{D_C1,1} , 3 }=DW{1,1};
        lTEMPLATE{ R{D_C2,1} , 3 }=DW{2,1};
        lTEMPLATE{ R{D_C3,1} , 3 }=DW{3,1};
        lTEMPLATE{ R{D_C4,1} , 3 }=DW{4,1};
        lTEMPLATE{ R{D_C5,1} , 3 }=DW{5,1};
        %clinical
        lTEMPLATE{ R{D_C1,2} , 3 }=DW{1,1};
        lTEMPLATE{ R{D_C2,2} , 3 }=DW{2,1};
        lTEMPLATE{ R{D_C3,2} , 3 }=DW{3,1};
        lTEMPLATE{ R{D_C4,2} , 3 }=DW{4,1};
        lTEMPLATE{ R{D_C5,2} , 3 }=DW{5,1};

        %pre-cancer stages
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %onset to first stage
        %lTEMPLATE{ R{onset,1}+I-1 ,R{onset,3}}=OnsetRates(I,1);

        %Baseline Prevalence
        lTEMPLATE{ R{PCH1_1618_Prev0,1}+I_loop-1 ,R{PCH1_1618_Prev0,3} }=PrevG(I_loop,1);
        lTEMPLATE{ R{PCH2_1618_Prev0,1}+I_loop-1 ,R{PCH2_1618_Prev0,3} }=PrevG(I_loop,2);
        lTEMPLATE{ R{PCH3_1618_Prev0,1}+I_loop-1 ,R{PCH3_1618_Prev0,3} }=PrevG(I_loop,3);
        lTEMPLATE{ R{PCL1_Prev0,1}+I_loop-1 ,R{PCL1_Prev0,3} }=PrevG(I_loop,4);
        lTEMPLATE{ R{PCL2_Prev0,1}+I_loop-1 ,R{PCL2_Prev0,3} }=PrevG(I_loop,5);
        lTEMPLATE{ R{PCL3_Prev0,1}+I_loop-1 ,R{PCL3_Prev0,3} }=PrevG(I_loop,6);
        lTEMPLATE{ R{PCR1_Prev0,1}+I_loop-1 ,R{PCR1_Prev0,3} }=PrevG(I_loop,7);
        lTEMPLATE{ R{PCM1_Prev0,1}+I_loop-1 ,R{PCM1_Prev0,3} }=PrevG(I_loop,8);
        lTEMPLATE{ R{PCR2_Prev0,1}+I_loop-1 ,R{PCR2_Prev0,3} }=PrevG(I_loop,13);
        lTEMPLATE{ R{PCH1_H_Prev0,1}+I_loop-1 ,R{PCH1_H_Prev0,3} }=PrevG(I_loop,10);
        lTEMPLATE{ R{PCH2_H_Prev0,1}+I_loop-1 ,R{PCH2_H_Prev0,3} }=PrevG(I_loop,11);
        lTEMPLATE{ R{PCH3_H_Prev0,1}+I_loop-1 ,R{PCH3_H_Prev0,3} }=PrevG(I_loop,12);
        lTEMPLATE{ R{PCR3_Prev0,1}+I_loop-1 ,R{PCR3_Prev0,3} }=PrevG(I_loop,9);
        %CFR rate
        lTEMPLATE{ R{PCH1_1618_CFR,1}+I_loop-1 ,R{PCH1_1618_CFR,3} }=0;
        lTEMPLATE{ R{PCH2_1618_CFR,1}+I_loop-1 ,R{PCH2_1618_CFR,3} }=0;
        lTEMPLATE{ R{PCH3_1618_CFR,1}+I_loop-1 ,R{PCH3_1618_CFR,3} }=0;
        lTEMPLATE{ R{PCH1_H_CFR,1}+I_loop-1 ,R{PCH1_H_CFR,3} }=0;
        lTEMPLATE{ R{PCH2_H_CFR,1}+I_loop-1 ,R{PCH2_H_CFR,3} }=0;
        lTEMPLATE{ R{PCH3_H_CFR,1}+I_loop-1 ,R{PCH3_H_CFR,3} }=0;
        lTEMPLATE{ R{PCL1_CFR,1}+I_loop-1 ,R{PCL1_CFR,3} }=0;
        lTEMPLATE{ R{PCL2_CFR,1}+I_loop-1 ,R{PCL2_CFR,3} }=0;
        lTEMPLATE{ R{PCL3_CFR,1}+I_loop-1 ,R{PCL3_CFR,3} }=0;
        lTEMPLATE{ R{PCR1_CFR,1}+I_loop-1 ,R{PCR1_CFR,3} }=0;
        lTEMPLATE{ R{PCR2_CFR,1}+I_loop-1 ,R{PCR2_CFR,3} }=0;
        lTEMPLATE{ R{PCR3_CFR,1}+I_loop-1 ,R{PCR3_CFR,3} }=0;
        lTEMPLATE{ R{PCM1_CFR,1}+I_loop-1 ,R{PCM1_CFR,3} }=0;

        %Progression
        lTEMPLATE{ R{PCH1_1618_Progr,1}+I_loop-1 ,R{PCH1_1618_Progr,3} }=0.0931;
        lTEMPLATE{ R{PCH2_1618_Progr,1}+I_loop-1 ,R{PCH2_1618_Progr,3} }=0.2107;
        lTEMPLATE{ R{PCH1_H_Progr,1}+I_loop-1 ,R{PCH1_H_Progr,3} }=0.0931;
        lTEMPLATE{ R{PCH2_H_Progr,1}+I_loop-1 ,R{PCH2_H_Progr,3} }=0.2107;
        lTEMPLATE{ R{PCL1_Progr,1}+I_loop-1 ,R{PCL1_Progr,3} }=0.0568;
        lTEMPLATE{ R{PCL2_Progr,1}+I_loop-1 ,R{PCL2_Progr,3} }=0.0921;
        lTEMPLATE{ R{PCM1_Progr,1}+I_loop-1 ,R{PCM1_Progr,3} }=1/Immunity;

        if I_loop<7
            lTEMPLATE{ R{PCH3_1618_Progr,1}+I_loop-1 ,R{PCH3_1618_Progr,3} }=0.0292;
            lTEMPLATE{ R{PCH3_H_Progr,1}+I_loop-1 ,R{PCH3_H_Progr,3} }=0.0292;
            lTEMPLATE{ R{PCL3_Progr,1}+I_loop-1 ,R{PCL3_Progr,3} }=0.0070;
        end
        if (I_loop<8 && I_loop>=7)
            lTEMPLATE{ R{PCH3_1618_Progr,1}+I_loop-1 ,R{PCH3_1618_Progr,3} }=0.0506;
            lTEMPLATE{ R{PCH3_H_Progr,1}+I_loop-1 ,R{PCH3_H_Progr,3} }=0.0506;
            lTEMPLATE{ R{PCL3_Progr,1}+I_loop-1 ,R{PCL3_Progr,3} }=0.0140;
        end
        if (I_loop<9 && I_loop>=8)
            lTEMPLATE{ R{PCH3_1618_Progr,1}+I_loop-1 ,R{PCH3_1618_Progr,3} }=0.1344;
            lTEMPLATE{ R{PCH3_H_Progr,1}+I_loop-1 ,R{PCH3_H_Progr,3} }=0.1344;
            lTEMPLATE{ R{PCL3_Progr,1}+I_loop-1 ,R{PCL3_Progr,3} }=0.0221;
        end
        if I_loop>=9
            lTEMPLATE{ R{PCH3_1618_Progr,1}+I_loop-1 ,R{PCH3_1618_Progr,3} }=0.1952; 
            lTEMPLATE{ R{PCH3_H_Progr,1}+I_loop-1 ,R{PCH3_H_Progr,3} }=0.1952; 
            lTEMPLATE{ R{PCL3_Progr,1}+I_loop-1 ,R{PCL3_Progr,3} }=0.0445;
        end

        %Regression
        lTEMPLATE{ R{PCH1_1618_Rem,1}+I_loop-1 ,R{PCH1_1618_Rem,3} }=0.2693;
        lTEMPLATE{ R{PCH2_1618_Rem1,1}+I_loop-1 ,R{PCH2_1618_Rem1,3} }=0.1188;
        lTEMPLATE{ R{PCH2_1618_Rem2,1}+I_loop-1 ,R{PCH2_1618_Rem2,3} }=0.1188;
        lTEMPLATE{ R{PCH3_1618_Rem1,1}+I_loop-1 ,R{PCH3_1618_Rem1,3} }=0.01715;
        lTEMPLATE{ R{PCH3_1618_Rem2,1}+I_loop-1 ,R{PCH3_1618_Rem2,3} }=0.01715;

        lTEMPLATE{ R{PCH1_H_Rem,1}+I_loop-1 ,R{PCH1_H_Rem,3} }=0.2693;
        lTEMPLATE{ R{PCH2_H_Rem1,1}+I_loop-1 ,R{PCH2_H_Rem1,3} }=0.1188;
        lTEMPLATE{ R{PCH2_H_Rem2,1}+I_loop-1 ,R{PCH2_H_Rem2,3} }=0.1188;
        lTEMPLATE{ R{PCH3_H_Rem1,1}+I_loop-1 ,R{PCH3_H_Rem1,3} }=0.01715;
        lTEMPLATE{ R{PCH3_H_Rem2,1}+I_loop-1 ,R{PCH3_H_Rem2,3} }=0.01715;



        lTEMPLATE{ R{PCL1_Rem,1}+I_loop-1 ,R{PCL1_Rem,3} }=0.2693;
        lTEMPLATE{ R{PCL2_Rem1,1}+I_loop-1 ,R{PCL2_Rem1,3} }=0.1059;
        lTEMPLATE{ R{PCL2_Rem2,1}+I_loop-1 ,R{PCL2_Rem2,3} }=0.1059;
        lTEMPLATE{ R{PCL3_Rem1,1}+I_loop-1 ,R{PCL3_Rem1,3} }=0.0704;
        lTEMPLATE{ R{PCL3_Rem2,1}+I_loop-1 ,R{PCL3_Rem2,3} }=0.0704;
        lTEMPLATE{ R{PCR1_Rem,1}+I_loop-1 ,R{PCR1_Rem,3} }=0.0363;
        lTEMPLATE{ R{PCR2_Rem,1}+I_loop-1 ,R{PCR2_Rem,3} }=0.0363;
        lTEMPLATE{ R{PCR3_Rem,1}+I_loop-1 ,R{PCR3_Rem,3} }=0.0363;
    end

end

function [lTEMPLATE2]=WriteTemplateM(RES2,TEMPLATEM)
    
    lTEMPLATE2=TEMPLATEM;
    PrevM=RES2{1};
    T=RES2{2};
    OnsetM=RES2{3};

    OnsetR_1618=1;OnsetR_HR=2;OnsetR_LR=3;
    Prev_1618=4;Prev_HR=5;Prev_LR=6;Prev_Im=7;
    Reg_1618=8;Reg_HR=9;Reg_LR=10;Reg_Im=11;

    %Onset rate of HPV
    R{OnsetR_1618,1}=9;R{OnsetR_1618,2}=20;R{OnsetR_1618,3}=7;
    R{OnsetR_HR,1}=9;R{OnsetR_HR,2}=20;R{OnsetR_HR,3}=10;
    R{OnsetR_LR,1}=9;R{OnsetR_LR,2}=20;R{OnsetR_LR,3}=13;

    %Regression rate
    R{Reg_1618,1}=75;R{Reg_1618,2}=86;R{Reg_1618,3}=9;
    R{Reg_HR,1}=97;R{Reg_HR,2}=108;R{Reg_HR,3}=9;
    R{Reg_LR,1}=119;R{Reg_LR,2}=130;R{Reg_LR,3}=9;
    R{Reg_Im,1}=141;R{Reg_Im,2}=152;R{Reg_Im,3}=9;

    %Prevalence
    P{Prev_1618,1}=75;P{Prev_1618,2}=86;P{Prev_1618,3}=4;
    P{Prev_HR,1}=97;P{Prev_HR,2}=108;P{Prev_HR,3}=4;
    P{Prev_LR,1}=119;P{Prev_LR,2}=130;P{Prev_LR,3}=4;
    P{Prev_Im,1}=141;P{Prev_Im,2}=152;P{Prev_Im,3}=4;

    for I_loop=1:12
        lTEMPLATE2{ R{OnsetR_1618,1}+I_loop-1 ,R{OnsetR_1618,3}}=OnsetM(I_loop,1);  
        lTEMPLATE2{ R{OnsetR_HR,1}+I_loop-1 ,R{OnsetR_HR,3}}=OnsetM(I_loop,2);  
        lTEMPLATE2{ R{OnsetR_LR,1}+I_loop-1 ,R{OnsetR_LR,3}}=OnsetM(I_loop,3);  

        lTEMPLATE2{ R{Reg_1618,1}+I_loop-1 ,R{Reg_1618,3}}=0.0363; 
        lTEMPLATE2{ R{Reg_HR,1}+I_loop-1 ,R{Reg_HR,3}}=0.0363; 
        lTEMPLATE2{ R{Reg_LR,1}+I_loop-1 ,R{Reg_LR,3}}=0.0363; 
        lTEMPLATE2{ R{Reg_Im,1}+I_loop-1 ,R{Reg_Im,3}}=0.1; 

        lTEMPLATE2{ P{Prev_1618,1}+I_loop-1 ,P{Prev_1618,3}}=PrevM(I_loop,1); 
        lTEMPLATE2{ P{Prev_HR,1}+I_loop-1 ,P{Prev_HR,3}}=PrevM(I_loop,2); 
        lTEMPLATE2{ P{Prev_LR,1}+I_loop-1 ,P{Prev_LR,3}}=PrevM(I_loop,3); 
        lTEMPLATE2{ P{Prev_Im,1}+I_loop-1 ,P{Prev_Im,3}}=PrevM(I_loop,4); 

    end
end

end

% 
% function ProcessModels_Cervical
% 
% [~,~,TEMPLATEW]=xlsread('CervicalCancerTemplate.xls','CervicalCancerTemplate');
% [~,~,TEMPLATEM]=xlsread('CervicalCancerMen.xlsx');
% [~,~,INTERVENTIONS]=xlsread('CervicalCancerTemplate.xls','Interventions');
% [~,~,DWeights]=xlsread('CervicalCancerTemplate.xls','DW');
% 
% %[~,~,TEMPLATEW]=xlsread('CervicalCancerTemplate-2016-02-03.xls','CervicalCancerTemplate');
% %[~,~,TEMPLATEM]=xlsread('CervicalCancerMen-2016-02-03.xlsx');
% %[~,~,INTERVENTIONS]=xlsread('CervicalCancerTemplate-2016-02-03.xls','Interventions');
% %[~,~,DWeights]=xlsread('CervicalCancerTemplate-2016-02-03.xls','DW');
% 
% %Intervention name
% INTVN{1,1}=INTERVENTIONS{6,2};
% INTVN{2,1}=INTERVENTIONS{8,2};
% INTVN{3,1}=INTERVENTIONS{10,2};
% INTVN{4,1}=INTERVENTIONS{12,2};
% INTVN{5,1}=INTERVENTIONS{14,2};
% 
% %MORT
% INTVN{1,2}=INTERVENTIONS{6,5};
% INTVN{2,2}=INTERVENTIONS{8,5};
% INTVN{3,2}=INTERVENTIONS{10,5};
% INTVN{4,2}=INTERVENTIONS{12,5};
% INTVN{5,2}=INTERVENTIONS{14,5};
% 
% %DW
% DW{1,1}=DWeights{1,2};
% DW{2,1}=DWeights{2,2};
% DW{3,1}=DWeights{3,2};
% DW{4,1}=DWeights{4,2};
% DW{5,1}=DWeights{5,2};
% 
% %DW
% DW{1,2}=DWeights{1,3};
% DW{2,2}=DWeights{2,3};
% DW{3,2}=DWeights{3,3};
% DW{4,2}=DWeights{4,3};
% DW{5,2}=DWeights{5,3};
% 
% 
% %pal care
% DW{6,2}=DWeights{7,3};
% DW{7,2}=DWeights{8,3};
% 
% %pal care
% DW{6,1}=DW{4,1};
% DW{7,1}=DW{4,1};
% 
% [~,~,RAW_D]=xlsread('res4.xlsx');
% 
% AgeHPVprev =[
% [16	25	35	45	55	65]
% [24	34	44	54	64	70]];
% 
% [~,Countries] = xlsread('DataInci.xlsx','B3:P3');
% 
% [StageDistMatrix] = xlsread('DataInci.xlsx','B32:P36');  
% StageDistMatrix=StageDistMatrix';
% [prevalenceActualMatrix] = xlsread('DataInci.xlsx','B4:P13');
%  prevalenceActualMatrix = prevalenceActualMatrix'/1000;               
% [distHPV] = xlsread('DataInci.xlsx','B42:P43'); 
% distHPV=distHPV';
% [distHigh] = xlsread('DataInci.xlsx','B48:P49'); 
% distHigh=distHigh';
% AgeArrayL=[1 15 40 45 50  55	60	65 70 75];
% AgeArrayU=[14 39 44  49 54	59	64	69 74 100];
% 
% [cancerMortalityArrayMatrix] = xlsread('DataInci.xlsx','B18:P27'); 
% cancerMortalityArrayMatrix = cancerMortalityArrayMatrix'/1000;
% % for I=2 for AFRE, 
% %      =12 for SEARB, 
% %      =16 for Peru
% for I=2
%     Countries{I}
%     ind=find(strcmp(RAW_D(:,1),Countries{I}));
%     DATA=RAW_D(ind,6);
%     
%     pA=cell2mat(DATA(163:243))'; %Population of female by age
%     pB=cell2mat(DATA(82:162))'; %%Population of male by age
%     MxW=cell2mat(DATA(406:486))'; %Natrual mortality rate by age for female
%     MxM=cell2mat(DATA(325:405))'; %Natrual mortality rate by age for male
%     Births=cell2mat(DATA(487:487)); %Birth population
%     xa=101;
%     pA(end+1:101)=0;
%     pB(end+1:101)=0;
%     MxW(end+1:101)=0;
%     MxM(end+1:101)=0;
%     PAR{1}=pA/sum(pA);
%     PAR{2}=MxW;
%     PAR{3}=Births/sum(pA);
%     PAR{4}=prevalenceActualMatrix(I,:);
%     PAR{5}=cancerMortalityArrayMatrix(I,:);
%     PAR{6}= StageDistMatrix(I,:);
%     PAR{7}=AgeArrayL(1,:);
%     PAR{8}=AgeArrayU(1,:);
%     PAR{9}=PAR{4};
%     PAR{10} =distHPV(I,:);
%     PAR{11} = AgeHPVprev;
%     PAR{12} = distHigh(I,:);
%     PAR{13}=I;
%     PAR{14}=Countries(I);
%     PAR{15}=pB/sum(pB);
%     PAR{16}=MxM;
%     [RES,RES2]=MarkovChainCountry2(PAR);
%     
%      [C_TEMPLATE]=WriteTemplateW(RES,TEMPLATEW);
%      xlswrite('Temp Output1.xlsx',C_TEMPLATE);
%      
%      Output2=zeros(749,52);
%      Output2(:)=NaN;
%      Output2=num2cell(Output2);
%      [~,~,Output1]=xlsread('Temp Output1.xlsx');
%      Output2(1:133,1:52)=Output1(1:133,1:52);
%      Output2(134:155,1:52)=Output1(706:727,1:52);
%      Output2(156:243,1:52)=Output1(882:969,1:52);
%      Output2(244:309,1:52)=Output1(640:705,1:52);
%      Output2(310:330,1:52)=Output1(970:990,1:52);
%      Output2(331,1:52)=Output1(969,1:52);
%      Output2(332:353,1:52)=Output1(728:749,1:52);
%      Output2(354:705,1:52)=Output1(134:485,1:52);
%      xlswrite('Output for Spectrum1.xlsx',Output2);
%     
%      [D_TEMPLATE2]=WriteTemplateM(RES2,TEMPLATEM);
%      xlswrite('Output for men1.xlsx',D_TEMPLATE2);
%      
%      
%      
%      
%      
%      
%     'Output for Spectrum written done'
% %      [C_TEMPLATE]=WriteTemplate(RES,TEMPLATE);
% %      xlswrite('Output for Spectrum.xlsx',C_TEMPLATE)
% %     'written'
% %     %pause
%     
% end
% 
% 
% function [lTEMPLATE]=WriteTemplateW(RES,TEMPLATE)
% 
% lTEMPLATE=TEMPLATE;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Data from Markov and Cohort model
% 
% PrevS=RES{1};
% PrevC=RES{2};
% Incidence=RES{3};
% OnsetRates=RES{4}';
% OnsetRates(end)=OnsetRates(end-1);
% DiagnosticRate=1./RES{5}
% DwellRate=1./RES{6}
% MortRate= RES{7};%Relative Mortality without treatment
% MortRateDF= RES{8};% Disease free mortality
% Impact = RES{9};%Relative Mortality with treatment
% TxCoverage = RES{10};% Current coverage of treatment
% PrevG=RES{11};
% HPVonsetG=RES{12};
% PrevS(1:4,:)= 0;
% PrevC(1:4,:) =0;
% OnsetRates(1:4)= 0;
% %Impact(:,1:6)=0;
% Immunity=RES{13};
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %Matrix with index to cell position in TEMPLATE;
% 
% onset=1;
% %Pre-clinical records
% PC1_Prev0=2;PC1_CFR=3;PC1_Screen=4;PC1_Progr=5;
% PC2_Prev0=6;PC2_CFR=7;PC2_Screen=8;PC2_Progr=9;
% PC3_Prev0=10;PC3_CFR=11;PC3_Screen=12;PC3_Progr=13;
% PC4_Prev0=14;PC4_CFR=15;PC4_Screen=16;PC4_Progr=17;
% PC5_Prev0=18;PC5_CFR=19;PC5_Screen=20;
% 
% %Clinical records
% C1_Prev0=21;C1_CFR=22;C1_Progr=23;
% C2_Prev0=24;C2_CFR=25;C2_Progr=26;
% C3_Prev0=27;C3_CFR=28;C3_Progr=29;
% C4_Prev0=30;C4_CFR=31;C4_Progr=32;
% C5_Prev0=33;C5_CFR=34;
% 
% %Treatment records
% I1_CFR=35;
% I2_CFR=36;
% I3_CFR=37;
% I4_CFR=38;
% I5_CFR=39;
% 
% %Palliative + Treatment records
% I6_CFR_Function=40;
% I7_CFR_Function=41;
% 
% %Palliative records
% I8_Function=42;
% I9_Function=43;
% 
% 
% %Screening
% I1_Screen=44;
% I2_Screen=45;
% I3_Screen=46;
% 
% %Disability weight records
% D_C1=47;
% D_C2=48;
% D_C3=49;
% D_C4=50;
% D_C5=51;
% 
% % Treatment coverage records
% Tx_C1 = 52;
% Tx_C2 = 53;
% Tx_C3 = 54;
% Tx_C4 = 55;
% Tx_C5 = 56;
% 
% Tx_C6 = 57;
% Tx_C7 = 58;
% Tx_C8 = 59;
% Tx_C9 = 60;
% Tx_C10 = 60;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Disease Free
% DSF1_incH=95;DSF1_incL=96;
% DSF1_inc1618=97;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Pre-cancer Stages
% %HPV 16/18
% PCH1_1618_Prev0=61;PCH1_1618_CFR=62;PCH1_1618_Rem=63;PCH1_1618_Progr=64;
% PCH2_1618_Prev0=65;PCH2_1618_CFR=66;PCH2_1618_Rem1=67;PCH2_1618_Rem2=68;PCH2_1618_Progr=69;
% PCH3_1618_Prev0=70;PCH3_1618_CFR=71;PCH3_1618_Rem1=72;PCH3_1618_Rem2=73;PCH3_1618_Progr=74;
% 
% %HPV High Risk
% PCH1_H_Prev0=104;PCH1_H_CFR=105;PCH1_H_Rem=106;PCH1_H_Progr=107;
% PCH2_H_Prev0=108;PCH2_H_CFR=109;PCH2_H_Rem1=110;PCH2_H_Rem2=111;PCH2_H_Progr=112;
% PCH3_H_Prev0=113;PCH3_H_CFR=114;PCH3_H_Rem1=115;PCH3_H_Rem2=116;PCH3_H_Progr=117;
% 
% %HPV Low Risk
% PCL1_Prev0=75;PCL1_CFR=76;PCL1_Rem=77;PCL1_Progr=78;
% PCL2_Prev0=79;PCL2_CFR=80;PCL2_Rem1=81;PCL2_Rem2=82;PCL2_Progr=83;
% PCL3_Prev0=84;PCL3_CFR=85;PCL3_Rem1=86;PCL3_Rem2=87;PCL3_Progr=88;
% 
% %HPV Regression
% PCR1_Prev0=89;PCR1_CFR=90;PCR1_Rem=91;
% PCR2_Prev0=98;PCR2_CFR=99;PCR2_Rem=100;
% PCR3_Prev0=101;PCR3_CFR=102;PCR3_Rem=103;
% %Immunity
% PCM1_Prev0=92;PCM1_CFR=93;PCM1_Progr=94;
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% %1=row 1, 2=row end, 3=row col
% %pre-cancer stages
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %HPV 16/18
% R{PCH1_1618_Prev0,1}=75;R{PCH1_1618_Prev0,2}=86;R{PCH1_1618_Prev0,3}=4;
% R{PCH1_1618_CFR,1}=75;R{PCH1_1618_CFR,2}=86;R{PCH1_1618_CFR,3}=6;
% R{PCH1_1618_Rem,1}=75;R{PCH1_1618_Rem,2}=86;R{PCH1_1618_Rem,3}=9;
% R{PCH1_1618_Progr,1}=75;R{PCH1_1618_Progr,2}=86;R{PCH1_1618_Progr,3}=12;
% 
% %CIN1 16/18
% R{PCH2_1618_Prev0,1}=97;R{PCH2_1618_Prev0,2}=108;R{PCH2_1618_Prev0,3}=4;
% R{PCH2_1618_CFR,1}=97;R{PCH2_1618_CFR,2}=108;R{PCH2_1618_CFR,3}=6;
% R{PCH2_1618_Rem1,1}=97;R{PCH2_1618_Rem1,2}=108;R{PCH2_1618_Rem1,3}=9;
% R{PCH2_1618_Rem2,1}=97;R{PCH2_1618_Rem2,2}=108;R{PCH2_1618_Rem2,3}=12;
% R{PCH2_1618_Progr,1}=97;R{PCH2_1618_Progr,2}=108;R{PCH2_1618_Progr,3}=15;
% 
% %CIN2/3 16/18
% R{PCH3_1618_Prev0,1}=119;R{PCH3_1618_Prev0,2}=130;R{PCH3_1618_Prev0,3}=4;
% R{PCH3_1618_CFR,1}=119;R{PCH3_1618_CFR,2}=130;R{PCH3_1618_CFR,3}=6;
% R{PCH3_1618_Rem1,1}=119;R{PCH3_1618_Rem1,2}=130;R{PCH3_1618_Rem1,3}=9;
% R{PCH3_1618_Rem2,1}=119;R{PCH3_1618_Rem2,2}=130;R{PCH3_1618_Rem2,3}=12;
% R{PCH3_1618_Progr,1}=119;R{PCH3_1618_Progr,2}=130;R{PCH3_1618_Progr,3}=15;
% 
% %HPV High Risk
% R{PCH1_H_Prev0,1}=889;R{PCH1_H_Prev0,2}=900;R{PCH1_H_Prev0,3}=4;
% R{PCH1_H_CFR,1}=889;R{PCH1_H_CFR,2}=900;R{PCH1_H_CFR,3}=6;
% R{PCH1_H_Rem,1}=889;R{PCH1_H_Rem,2}=900;R{PCH1_H_Rem,3}=9;
% R{PCH1_H_Progr,1}=889;R{PCH1_H_Progr,2}=900;R{PCH1_H_Progr,3}=12;
% 
% %CIN1 High Risk
% R{PCH2_H_Prev0,1}=911;R{PCH2_H_Prev0,2}=922;R{PCH2_H_Prev0,3}=4;
% R{PCH2_H_CFR,1}=911;R{PCH2_H_CFR,2}=922;R{PCH2_H_CFR,3}=6;
% R{PCH2_H_Rem1,1}=911;R{PCH2_H_Rem1,2}=922;R{PCH2_H_Rem1,3}=9;
% R{PCH2_H_Rem2,1}=911;R{PCH2_H_Rem2,2}=922;R{PCH2_H_Rem2,3}=12;
% R{PCH2_H_Progr,1}=911;R{PCH2_H_Progr,2}=922;R{PCH2_H_Progr,3}=15;
% 
% %CIN2/3 High Risk
% R{PCH3_H_Prev0,1}=933;R{PCH3_H_Prev0,2}=944;R{PCH3_H_Prev0,3}=4;
% R{PCH3_H_CFR,1}=933;R{PCH3_H_CFR,2}=944;R{PCH3_H_CFR,3}=6;
% R{PCH3_H_Rem1,1}=933;R{PCH3_H_Rem1,2}=944;R{PCH3_H_Rem1,3}=9;
% R{PCH3_H_Rem2,1}=933;R{PCH3_H_Rem2,2}=944;R{PCH3_H_Rem2,3}=12;
% R{PCH3_H_Progr,1}=933;R{PCH3_H_Progr,2}=944;R{PCH3_H_Progr,3}=15;
% 
% %HPV Low Risk
% R{PCL1_Prev0,1}=647;R{PCL1_Prev0,2}=658;R{PCL1_Prev0,3}=4;
% R{PCL1_CFR,1}=647;R{PCL1_CFR,2}=658;R{PCL1_CFR,3}=6;
% R{PCL1_Rem,1}=647;R{PCL1_Rem,2}=658;R{PCL1_Rem,3}=9;
% R{PCL1_Progr,1}=647;R{PCL1_Progr,2}=658;R{PCL1_Progr,3}=12;
% 
% %CIN1 Low Risk
% R{PCL2_Prev0,1}=669;R{PCL2_Prev0,2}=680;R{PCL2_Prev0,3}=4;
% R{PCL2_CFR,1}=669;R{PCL2_CFR,2}=680;R{PCL2_CFR,3}=6;
% R{PCL2_Rem1,1}=669;R{PCL2_Rem1,2}=680;R{PCL2_Rem1,3}=9;
% R{PCL2_Rem2,1}=669;R{PCL2_Rem2,2}=680;R{PCL2_Rem2,3}=12;
% R{PCL2_Progr,1}=669;R{PCL2_Progr,2}=680;R{PCL2_Progr,3}=15;
% 
% %CIN2/3 Low Risk
% R{PCL3_Prev0,1}=691;R{PCL3_Prev0,2}=702;R{PCL3_Prev0,3}=4;
% R{PCL3_CFR,1}=691;R{PCL3_CFR,2}=702;R{PCL3_CFR,3}=6;
% R{PCL3_Rem1,1}=691;R{PCL3_Rem1,2}=702;R{PCL3_Rem1,3}=9;
% R{PCL3_Rem2,1}=691;R{PCL3_Rem2,2}=702;R{PCL3_Rem2,3}=12;
% R{PCL3_Progr,1}=691;R{PCL3_Progr,2}=702;R{PCL3_Progr,3}=15;
% 
% %HPV Regression
% R{PCR1_Prev0,1}=713;R{PCR1_Prev0,2}=724;R{PCR1_Prev0,3}=4;
% R{PCR1_CFR,1}=713;R{PCR1_CFR,2}=724;R{PCR1_CFR,3}=6;
% R{PCR1_Rem,1}=713;R{PCR1_Rem,2}=724;R{PCR1_Rem,3}=9;
% R{PCR2_Prev0,1}=955;R{PCR2_Prev0,2}=966;R{PCR2_Prev0,3}=4;
% R{PCR2_CFR,1}=955;R{PCR2_CFR,2}=966;R{PCR2_CFR,3}=6;
% R{PCR2_Rem,1}=955;R{PCR2_Rem,2}=966;R{PCR2_Rem,3}=9;
% R{PCR3_Prev0,1}=977;R{PCR3_Prev0,2}=988;R{PCR3_Prev0,3}=4;
% R{PCR3_CFR,1}=977;R{PCR3_CFR,2}=988;R{PCR3_CFR,3}=6;
% R{PCR3_Rem,1}=977;R{PCR3_Rem,2}=988;R{PCR3_Rem,3}=9;
% 
% 
% 
% 
% %Immunity
% R{PCM1_Prev0,1}=735;R{PCM1_Prev0,2}=746;R{PCM1_Prev0,3}=4;
% R{PCM1_CFR,1}=735;R{PCM1_CFR,2}=746;R{PCM1_CFR,3}=6;
% R{PCM1_Progr,1}=735;R{PCM1_Progr,2}=746;R{PCM1_Progr,3}=9;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% 
% %cancer stages
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %disease free
% R{onset,1}=9;R{onset,2}=20;R{onset,3}=7;
% R{DSF1_inc1618,1}=9;R{DSF1_inc1618,2}=20;R{DSF1_inc1618,3}=7;
% R{DSF1_incH,1}=9;R{DSF1_incH,2}=20;R{DSF1_incH,3}=10;
% R{DSF1_incL,1}=9;R{DSF1_incL,2}=20;R{DSF1_incL,3}=13;
% %pre clinical
% %stage 0
% R{PC1_Prev0,1}=141;R{PC1_Prev0,2}=152;R{PC1_Prev0,3}=4;
% R{PC1_CFR,1}=141;R{PC1_CFR,2}=152;R{PC1_CFR,3}=6;
% R{PC1_Screen,1}=141;R{PC1_Screen,2}=152;R{PC1_Screen,3}=9;
% R{PC1_Progr,1}=141;R{PC1_Progr,2}=152;R{PC1_Progr,3}=12;
% 
% %stage 1
% R{PC2_Prev0,1}=163;R{PC2_Prev0,2}=174;R{PC2_Prev0,3}=4;
% R{PC2_CFR,1}=163;R{PC2_CFR,2}=174;R{PC2_CFR,3}=6;
% R{PC2_Screen,1}=163;R{PC2_Screen,2}=174;R{PC2_Screen,3}=9;
% R{PC2_Progr,1}=163;R{PC2_Progr,2}=174;R{PC2_Progr,3}=12;
% 
% %stage 2
% R{PC3_Prev0,1}=185;R{PC3_Prev0,2}=196;R{PC3_Prev0,3}=4;
% R{PC3_CFR,1}=185;R{PC3_CFR,2}=196;R{PC3_CFR,3}=6;
% R{PC3_Screen,1}=185;R{PC3_Screen,2}=196;R{PC3_Screen,3}=9;
% R{PC3_Progr,1}=185;R{PC3_Progr,2}=196;R{PC3_Progr,3}=12;
% 
% %stage 3
% R{PC4_Prev0,1}=207;R{PC4_Prev0,2}=218;R{PC4_Prev0,3}=4;
% R{PC4_CFR,1}=207;R{PC4_CFR,2}=218;R{PC4_CFR,3}=6;
% R{PC4_Screen,1}=207;R{PC4_Screen,2}=218;R{PC4_Screen,3}=9;
% R{PC4_Progr,1}=207;R{PC4_Progr,2}=218;R{PC4_Progr,3}=12;  
% 
% %stage 4
% R{PC5_Prev0,1}=229;R{PC5_Prev0,2}=240;R{PC5_Prev0,3}=4;
% R{PC5_CFR,1}=229;R{PC5_CFR,2}=240;R{PC5_CFR,3}=6;
% R{PC5_Screen,1}=229;R{PC5_Screen,2}=240;R{PC5_Screen,3}=9;
% %R{PC4_Progr,1}=207;R{PC4_Progr,2}=218;R{PC4_Progr,3}=12;  
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %clinical 0
% R{C1_Prev0,1}=251;R{C1_Prev0,2}=262;R{C1_Prev0,3}=4;
% R{C1_CFR,1}=251;R{C1_CFR,2}=262;R{C1_CFR,3}=6;
% R{C1_Progr,1}=251;R{C1_Progr,2}=262;R{C1_Progr,3}=9;
%     
% %clinical 1
% R{C2_Prev0,1}=273;R{C2_Prev0,2}=284;R{C2_Prev0,3}=4;
% R{C2_CFR,1}=273;R{C2_CFR,2}=284;R{C2_CFR,3}=6;
% R{C2_Progr,1}=273;R{C2_Progr,2}=284;R{C2_Progr,3}=9;
% 
% %clinical 2
% R{C3_Prev0,1}=295;R{C3_Prev0,2}=306;R{C3_Prev0,3}=4;
% R{C3_CFR,1}=295;R{C3_CFR,2}=306;R{C3_CFR,3}=6;
% R{C3_Progr,1}=295;R{C3_Progr,2}=306;R{C3_Progr,3}=9;
% 
% %clinical 3
% R{C4_Prev0,1}=317;R{C4_Prev0,2}=328;R{C4_Prev0,3}=4;
% R{C4_CFR,1}=317;R{C4_CFR,2}=328;R{C4_CFR,3}=6;
% R{C4_Progr,1}=317;R{C4_Progr,2}=328;R{C4_Progr,3}=9;
% 
% %clinical 4
% R{C5_Prev0,1}=339;R{C5_Prev0,2}=350;R{C5_Prev0,3}=4;
% R{C5_CFR,1}=339;R{C5_CFR,2}=350;R{C5_CFR,3}=6;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Interventions 7-CFR, 14 Function-DW
% R{I1_CFR,1}=383;R{I1_CFR,2}=383;R{I1_CFR,3}=7;R{I1_CFR,4}=14;
% R{I2_CFR,1}=405;R{I2_CFR,2}=405;R{I2_CFR,3}=7;R{I2_CFR,4}=14;
% R{I3_CFR,1}=427;R{I3_CFR,2}=427;R{I3_CFR,3}=7;R{I3_CFR,4}=14;
% R{I4_CFR,1}=449;R{I4_CFR,2}=449;R{I4_CFR,3}=7;R{I4_CFR,4}=14;
% R{I5_CFR,1}=471;R{I5_CFR,2}=471;R{I5_CFR,3}=7;R{I5_CFR,4}=14;
% 
% R{I6_CFR_Function,1}=493;R{I6_CFR_Function,2}=493;R{I6_CFR_Function,3}=7;R{I6_CFR_Function,4}=14;
% R{I7_CFR_Function,1}=515;R{I7_CFR_Function,2}=515;R{I7_CFR_Function,3}=7;R{I7_CFR_Function,4}=14;
% 
% R{I8_Function,1}=537;R{I8_Function,2}=537;R{I8_Function,4}=7;
% R{I9_Function,1}=559;R{I9_Function,2}=559;R{I9_Function,4}=7;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Screening Interventions
% %bi-annual screening
% R{I1_Screen,1}=581;R{I1_Screen,2}=581;R{I1_Screen,3}=7;R{I1_Screen,4}=14;R{I1_Screen,5}=21;R{I1_Screen,6}=28;
% R{I2_Screen,1}=603;R{I2_Screen,2}=603;R{I2_Screen,3}=7;R{I2_Screen,4}=14;R{I2_Screen,5}=21;R{I2_Screen,6}=28;
% %clinical examination
% R{I3_Screen,1}=625;R{I3_Screen,2}=625;R{I3_Screen,3}=7;R{I3_Screen,4}=14;R{I3_Screen,5}=21;R{I3_Screen,6}=28;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Disability weights
% %using row indices only
% R{D_C1,1}=137;R{D_C1,2}=247;
% R{D_C2,1}=159;R{D_C2,2}=269;
% R{D_C3,1}=181;R{D_C3,2}=291;
% R{D_C4,1}=203;R{D_C4,2}=313;
% R{D_C5,1}=225;R{D_C5,2}=335;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Current and scale-up coverage
% R{Tx_C1,1}=383;R{Tx_C1,2}=383;R{Tx_C1,3}=383;R{Tx_C1,4}=8;R{Tx_C1,5}=9;R{Tx_C1,6}=15;R{Tx_C1,7}=16;
% R{Tx_C2,1}=405;R{Tx_C2,2}=405;R{Tx_C2,3}=405;R{Tx_C2,4}=8;R{Tx_C2,5}=9;R{Tx_C2,6}=15;R{Tx_C2,7}=16;
% R{Tx_C3,1}=427;R{Tx_C3,2}=427;R{Tx_C3,3}=427;R{Tx_C3,4}=8;R{Tx_C3,5}=9;R{Tx_C3,6}=15;R{Tx_C3,7}=16;
% R{Tx_C4,1}=449;R{Tx_C4,2}=449;R{Tx_C4,3}=449;R{Tx_C4,4}=8;R{Tx_C4,5}=9;R{Tx_C4,6}=15;R{Tx_C4,7}=16;
% R{Tx_C5,1}=471;R{Tx_C5,2}=471;R{Tx_C5,3}=471;R{Tx_C5,4}=8;R{Tx_C5,5}=9;R{Tx_C5,6}=15;R{Tx_C5,7}=16;
% 
% 
% R{Tx_C6,1}=493;R{Tx_C6,2}=493;R{Tx_C6,3}=493;R{Tx_C6,4}=8;R{Tx_C6,5}=9;R{Tx_C6,6}=15;R{Tx_C6,7}=16;
% R{Tx_C7,1}=515;R{Tx_C7,2}=515;R{Tx_C7,3}=515;R{Tx_C7,4}=8;R{Tx_C7,5}=9;R{Tx_C7,6}=15;R{Tx_C7,7}=16;
% R{Tx_C8,1}=537;R{Tx_C8,2}=537;R{Tx_C8,3}=537;R{Tx_C8,4}=8;R{Tx_C8,5}=9;R{Tx_C8,6}=15;R{Tx_C8,7}=16;
% R{Tx_C9,1}=559;R{Tx_C9,2}=559;R{Tx_C9,3}=559;R{Tx_C9,4}=8;R{Tx_C9,5}=9;R{Tx_C9,6}=15;R{Tx_C9,7}=16;
% 
% 
% for I=1:12
%  
% %onset to first stage
% lTEMPLATE{ R{onset,1}+I-1 ,R{onset,3}}=OnsetRates(I,1);    
% lTEMPLATE{ R{DSF1_inc1618,1}+I-1 ,R{DSF1_inc1618,3} }=HPVonsetG(I,1);
% lTEMPLATE{ R{DSF1_incH,1}+I-1 ,R{DSF1_incH,3} }=HPVonsetG(I,2);
% lTEMPLATE{ R{DSF1_incL,1}+I-1 ,R{DSF1_incL,3} }=HPVonsetG(I,3);
% 
% 
% %progression
% lTEMPLATE{ R{PC1_Prev0,1}+I-1 ,R{PC1_Prev0,3} }=PrevS(I,1);
% lTEMPLATE{ R{PC2_Prev0,1}+I-1 ,R{PC2_Prev0,3} }=PrevS(I,2);
% lTEMPLATE{ R{PC3_Prev0,1}+I-1 ,R{PC3_Prev0,3} }=PrevS(I,3);
% lTEMPLATE{ R{PC4_Prev0,1}+I-1 ,R{PC4_Prev0,3} }=PrevS(I,4);
% lTEMPLATE{ R{PC5_Prev0,1}+I-1 ,R{PC5_Prev0,3} }=PrevS(I,5);
% 
% %progression
% lTEMPLATE{ R{PC1_Progr,1}+I-1 ,R{PC1_Progr,3} }=DwellRate(I,1);
% lTEMPLATE{ R{PC2_Progr,1}+I-1 ,R{PC2_Progr,3} }=DwellRate(I,2);
% lTEMPLATE{ R{PC3_Progr,1}+I-1 ,R{PC3_Progr,3} }=DwellRate(I,3);
% lTEMPLATE{ R{PC4_Progr,1}+I-1 ,R{PC4_Progr,3} }=DwellRate(I,4);
% 
% %screening process
% lTEMPLATE{ R{PC1_Screen,1}+I-1 ,R{PC1_Screen,3} }=DiagnosticRate(I,1);
% lTEMPLATE{ R{PC2_Screen,1}+I-1 ,R{PC2_Screen,3} }=DiagnosticRate(I,2);
% lTEMPLATE{ R{PC3_Screen,1}+I-1 ,R{PC3_Screen,3} }=DiagnosticRate(I,3);
% lTEMPLATE{ R{PC4_Screen,1}+I-1 ,R{PC4_Screen,3} }=DiagnosticRate(I,4);
% lTEMPLATE{ R{PC5_Screen,1}+I-1 ,R{PC5_Screen,3} }=DiagnosticRate(I,5);
% 
% %Cancer Mortality pre-Clinical
% lTEMPLATE{ R{PC1_CFR,1}+I-1 ,R{PC1_CFR,3} }=0;%MortRateDF(I,1);
% lTEMPLATE{ R{PC2_CFR,1}+I-1 ,R{PC2_CFR,3} }=0;%MortRateDF(I,1);
% lTEMPLATE{ R{PC3_CFR,1}+I-1 ,R{PC3_CFR,3} }=0;%MortRateDF(I,1);
% lTEMPLATE{ R{PC4_CFR,1}+I-1 ,R{PC4_CFR,3} }=0;%MortRate(4,I);
% lTEMPLATE{ R{PC5_CFR,1}+I-1 ,R{PC5_CFR,3} }=MortRate(5,I);
% 
% %Clinical
% %Base Cancer Mortality (after considering coverage of TX among those
% %diagnosed
% lTEMPLATE{ R{C1_CFR,1}+I-1 ,R{C1_CFR,3} }=MortRate(1,I);
% lTEMPLATE{ R{C2_CFR,1}+I-1 ,R{C2_CFR,3} }=MortRate(2,I);
% lTEMPLATE{ R{C3_CFR,1}+I-1 ,R{C3_CFR,3} }=MortRate(3,I);
% lTEMPLATE{ R{C4_CFR,1}+I-1 ,R{C4_CFR,3} }=MortRate(4,I);
% lTEMPLATE{ R{C5_CFR,1}+I-1 ,R{C5_CFR,3} }=MortRate(5,I);
% 
% %prev 0
% lTEMPLATE{ R{C1_Prev0,1}+I-1 ,R{C1_Prev0,3} }=PrevC(I,1);
% lTEMPLATE{ R{C2_Prev0,1}+I-1 ,R{C2_Prev0,3} }=PrevC(I,2);
% lTEMPLATE{ R{C3_Prev0,1}+I-1 ,R{C3_Prev0,3} }=PrevC(I,3);
% lTEMPLATE{ R{C4_Prev0,1}+I-1 ,R{C4_Prev0,3} }=PrevC(I,4);
% lTEMPLATE{ R{C5_Prev0,1}+I-1 ,R{C5_Prev0,3} }=PrevC(I,5);
% 
% %progression
% lTEMPLATE{ R{C1_Progr,1}+I-1 ,R{C1_Progr,3} }=0;%DwellRate(I,1);
% lTEMPLATE{ R{C2_Progr,1}+I-1 ,R{C2_Progr,3} }=0;%DwellRate(I,2);
% lTEMPLATE{ R{C3_Progr,1}+I-1 ,R{C3_Progr,3} }=0;%DwellRate(I,3);
% lTEMPLATE{ R{C4_Progr,1}+I-1 ,R{C4_Progr,3} }=0;%DwellRate(I,3);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Treatment impact on MORT
% lTEMPLATE{ R{I1_CFR,1}+I-1 ,R{I1_CFR,3} }=Impact(1,I);%INTVN{1,2}*100;
% lTEMPLATE{ R{I2_CFR,1}+I-1 ,R{I2_CFR,3} }=Impact(2,I);%INTVN{2,2}*100;
% lTEMPLATE{ R{I3_CFR,1}+I-1 ,R{I3_CFR,3} }=Impact(3,I);%INTVN{3,2}*100;
% lTEMPLATE{ R{I4_CFR,1}+I-1 ,R{I4_CFR,3} }=Impact(4,I);%INTVN{4,2}*100;
% lTEMPLATE{ R{I5_CFR,1}+I-1 ,R{I5_CFR,3} }=Impact(5,I);%INTVN{4,2}*100;
% 
% % Treatment+palliative care impact on Mort
% lTEMPLATE{ R{I6_CFR_Function,1}+I-1 ,R{I6_CFR_Function,3} }=Impact(5,I);%INTVN{4,2}*100;
% lTEMPLATE{ R{I7_CFR_Function,1}+I-1 ,R{I7_CFR_Function,3} }=Impact(5,I);%INTVN{4,2}*100;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Baseline treatment coverage for treatment
% lTEMPLATE{ R{Tx_C1,1}+I-1 ,R{Tx_C1,4} }  = TxCoverage(1,I)*100;
% lTEMPLATE{ R{Tx_C2,1}+I-1 ,R{Tx_C2,4} }  = TxCoverage(2,I)*100;   
% lTEMPLATE{ R{Tx_C3,1}+I-1 ,R{Tx_C3,4} }  = TxCoverage(3,I)*100;
% lTEMPLATE{ R{Tx_C4,1}+I-1 ,R{Tx_C4,4} }  = TxCoverage(4,I)*100;
% lTEMPLATE{ R{Tx_C5,1}+I-1 ,R{Tx_C5,4} }  = TxCoverage(5,I)*100;
% 
% lTEMPLATE{ R{Tx_C6,1}+I-1 ,R{Tx_C6,4} }  = TxCoverage(5,I)*100;
% lTEMPLATE{ R{Tx_C7,1}+I-1 ,R{Tx_C7,4} }  = TxCoverage(5,I)*100;
% lTEMPLATE{ R{Tx_C8,1}+I-1 ,R{Tx_C8,4} }  = TxCoverage(5,I)*100;
% lTEMPLATE{ R{Tx_C9,1}+I-1 ,R{Tx_C9,4} }  = TxCoverage(5,I)*100;
% 
% 
% %Default intervention treatment coverage= current + 80% of those not
% lTEMPLATE{ R{Tx_C1,1}+I-1 ,R{Tx_C1,5} }  = ((1-TxCoverage(1,I))*.8+TxCoverage(1,I)) *100;
% lTEMPLATE{ R{Tx_C2,1}+I-1 ,R{Tx_C2,5} }  = ((1-TxCoverage(2,I))*.8+TxCoverage(2,I)) *100;   
% lTEMPLATE{ R{Tx_C3,1}+I-1 ,R{Tx_C3,5} }  = ((1-TxCoverage(3,I))*.8+TxCoverage(3,I)) *100;
% lTEMPLATE{ R{Tx_C4,1}+I-1 ,R{Tx_C4,5} }  = ((1-TxCoverage(4,I))*.8+TxCoverage(4,I)) *100;
% lTEMPLATE{ R{Tx_C5,1}+I-1 ,R{Tx_C5,5} }  = ((1-TxCoverage(5,I))*.8+TxCoverage(5,I)) *100;
% 
% lTEMPLATE{ R{Tx_C6,1}+I-1 ,R{Tx_C6,5} }  = ((1-TxCoverage(5,I))*.8+TxCoverage(5,I)) *100;
% lTEMPLATE{ R{Tx_C7,1}+I-1 ,R{Tx_C7,5} }  = ((1-TxCoverage(5,I))*.8+TxCoverage(5,I)) *100;
% lTEMPLATE{ R{Tx_C8,1}+I-1 ,R{Tx_C8,5} }  = ((1-TxCoverage(5,I))*.8+TxCoverage(5,I)) *100;
% lTEMPLATE{ R{Tx_C9,1}+I-1 ,R{Tx_C9,5} }  = ((1-TxCoverage(5,I))*.8+TxCoverage(5,I)) *100;
% 
% % Baseline coverage for disability weight = same as Tx coverage
% lTEMPLATE{ R{Tx_C1,1}+I-1 ,R{Tx_C1,6} }  = TxCoverage(1,I)*100;
% lTEMPLATE{ R{Tx_C2,1}+I-1 ,R{Tx_C2,6} }  = TxCoverage(2,I)*100;   
% lTEMPLATE{ R{Tx_C3,1}+I-1 ,R{Tx_C3,6} }  = TxCoverage(3,I)*100;
% lTEMPLATE{ R{Tx_C4,1}+I-1 ,R{Tx_C4,6} }  = TxCoverage(4,I)*100;
% lTEMPLATE{ R{Tx_C5,1}+I-1 ,R{Tx_C5,6} }  = TxCoverage(5,I)*100;
% 
% lTEMPLATE{ R{Tx_C6,1}+I-1 ,R{Tx_C6,6} }  = TxCoverage(5,I)*100;
% lTEMPLATE{ R{Tx_C7,1}+I-1 ,R{Tx_C7,6} }  = TxCoverage(5,I)*100;
% 
% %Default intervention coverage for disability weight= current + 80% of those not
% lTEMPLATE{ R{Tx_C1,1}+I-1 ,R{Tx_C1,7} }  = ((1-TxCoverage(1,I))*.8+TxCoverage(1,I)) *100;
% lTEMPLATE{ R{Tx_C2,1}+I-1 ,R{Tx_C2,7} }  = ((1-TxCoverage(2,I))*.8+TxCoverage(2,I)) *100;   
% lTEMPLATE{ R{Tx_C3,1}+I-1 ,R{Tx_C3,7} }  = ((1-TxCoverage(3,I))*.8+TxCoverage(3,I)) *100;
% lTEMPLATE{ R{Tx_C4,1}+I-1 ,R{Tx_C4,7} }  = ((1-TxCoverage(4,I))*.8+TxCoverage(4,I)) *100;
% lTEMPLATE{ R{Tx_C5,1}+I-1 ,R{Tx_C5,7} }  = ((1-TxCoverage(5,I))*.8+TxCoverage(5,I)) *100;
% 
% lTEMPLATE{ R{Tx_C6,1}+I-1 ,R{Tx_C6,7} }  = ((1-TxCoverage(5,I))*.8+TxCoverage(5,I)) *100;
% lTEMPLATE{ R{Tx_C7,1}+I-1 ,R{Tx_C7,7} }  = ((1-TxCoverage(5,I))*.8+TxCoverage(5,I)) *100;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Treatment
% %impact on disability weights, see def of DW in start of function
% lTEMPLATE{ R{I1_CFR,1}+I-1 ,R{I1_CFR,4} }=-(DW{1,1}-DW{1,2})/DW{1,1}*100;
% lTEMPLATE{ R{I2_CFR,1}+I-1 ,R{I2_CFR,4} }=-(DW{2,1}-DW{2,2})/DW{2,1}*100;
% lTEMPLATE{ R{I3_CFR,1}+I-1 ,R{I3_CFR,4} }=-(DW{3,1}-DW{3,2})/DW{3,1}*100;
% lTEMPLATE{ R{I4_CFR,1}+I-1 ,R{I4_CFR,4} }=-(DW{4,1}-DW{4,2})/DW{4,1}*100;
% lTEMPLATE{ R{I5_CFR,1}+I-1 ,R{I5_CFR,4} }=-(DW{5,1}-DW{5,2})/DW{5,1}*100;
% 
% %Treatment+ Palliative Care
% lTEMPLATE{ R{I6_CFR_Function,1}+I-1 ,R{I6_CFR_Function,4} }=-(DW{5,1}-DW{6,2})/DW{5,1}*100;
% lTEMPLATE{ R{I7_CFR_Function,1}+I-1 ,R{I7_CFR_Function,4} }=-(DW{5,1}-DW{7,2})/DW{5,1}*100;
% 
% %Palliative Care
% lTEMPLATE{ R{I8_Function,1}+I-1 ,R{I8_Function,4} }=-(DW{5,1}-DW{6,2})/DW{5,1}*100;
% lTEMPLATE{ R{I9_Function,1}+I-1 ,R{I9_Function,4} }=-(DW{5,1}-DW{7,2})/DW{5,1}*100;
% 
% 
% %Bi annual screening of 50-69, age match not perfect
% lTEMPLATE{ R{I1_Screen,1}+I-1 ,R{I1_Screen,3} }=100*0/2;
% lTEMPLATE{ R{I1_Screen,1}+I-1 ,R{I1_Screen,4} }=100*0/2;
% lTEMPLATE{ R{I1_Screen,1}+I-1 ,R{I1_Screen,5} }=100*0/2;
% lTEMPLATE{ R{I1_Screen,1}+I-1 ,R{I1_Screen,6} }=100*0/2;
% %lTEMPLATE{ R{I1_Screen,1}+I-1 ,R{I1_Screen,7} }=100*0/2;
% if(ismember(I,[9,10]))
% lTEMPLATE{ R{I1_Screen,1}+I-1 ,R{I1_Screen,3} }=100*1/2;
% lTEMPLATE{ R{I1_Screen,1}+I-1 ,R{I1_Screen,4} }=100*1/2;
% lTEMPLATE{ R{I1_Screen,1}+I-1 ,R{I1_Screen,5} }=100*1/2;
% lTEMPLATE{ R{I1_Screen,1}+I-1 ,R{I1_Screen,6} }=100*1/2;
% %lTEMPLATE{ R{I1_Screen,1}+I-1 ,R{I1_Screen,7} }=100*1/2;
% end
% 
% %Bi annual screening of 40-69, age match not perfect
% lTEMPLATE{ R{I2_Screen,1}+I-1 ,R{I2_Screen,3} }=100*0/2;
% lTEMPLATE{ R{I2_Screen,1}+I-1 ,R{I2_Screen,4} }=100*0/2;
% lTEMPLATE{ R{I2_Screen,1}+I-1 ,R{I2_Screen,5} }=100*0/2;
% lTEMPLATE{ R{I2_Screen,1}+I-1 ,R{I2_Screen,6} }=100*0/2;
% %lTEMPLATE{ R{I2_Screen,1}+I-1 ,R{I2_Screen,7} }=100*0/2;
% if(ismember(I,[8,9,10]))
% lTEMPLATE{ R{I2_Screen,1}+I-1 ,R{I2_Screen,3} }=100*1/2;
% lTEMPLATE{ R{I2_Screen,1}+I-1 ,R{I2_Screen,4} }=100*1/2;
% lTEMPLATE{ R{I2_Screen,1}+I-1 ,R{I2_Screen,5} }=100*1/2;
% lTEMPLATE{ R{I2_Screen,1}+I-1 ,R{I2_Screen,6} }=100*1/2;
% %lTEMPLATE{ R{I2_Screen,1}+I-1 ,R{I2_Screen,7} }=100*1/2;
% end
% 
% %Screening-by clinical exam
% lTEMPLATE{ R{I3_Screen,1}+I-1 ,R{I3_Screen,3} }=100*0/2;
% lTEMPLATE{ R{I3_Screen,1}+I-1 ,R{I3_Screen,4} }=100*0/2;
% lTEMPLATE{ R{I3_Screen,1}+I-1 ,R{I3_Screen,5} }=100*0/2;
% lTEMPLATE{ R{I3_Screen,1}+I-1 ,R{I3_Screen,6} }=100*0/2;
% %lTEMPLATE{ R{I3_Screen,1}+I-1 ,R{I3_Screen,7} }=100*0/2;
% if(ismember(I,[8,9,10]))%assume it takes longer than bi-anual
% lTEMPLATE{ R{I3_Screen,1}+I-1 ,R{I3_Screen,3} }=100*1/4;
% lTEMPLATE{ R{I3_Screen,1}+I-1 ,R{I3_Screen,4} }=100*1/4;
% lTEMPLATE{ R{I3_Screen,1}+I-1 ,R{I3_Screen,5} }=100*1/4;
% lTEMPLATE{ R{I3_Screen,1}+I-1 ,R{I3_Screen,6} }=100*1/4;
% %lTEMPLATE{ R{I3_Screen,1}+I-1 ,R{I3_Screen,7} }=100*1/4;
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% %Disability weights
% %preclinical
% lTEMPLATE{ R{D_C1,1} , 3 }=DW{1,1};
% lTEMPLATE{ R{D_C2,1} , 3 }=DW{2,1};
% lTEMPLATE{ R{D_C3,1} , 3 }=DW{3,1};
% lTEMPLATE{ R{D_C4,1} , 3 }=DW{4,1};
% lTEMPLATE{ R{D_C5,1} , 3 }=DW{5,1};
% %clinical
% lTEMPLATE{ R{D_C1,2} , 3 }=DW{1,1};
% lTEMPLATE{ R{D_C2,2} , 3 }=DW{2,1};
% lTEMPLATE{ R{D_C3,2} , 3 }=DW{3,1};
% lTEMPLATE{ R{D_C4,2} , 3 }=DW{4,1};
% lTEMPLATE{ R{D_C5,2} , 3 }=DW{5,1};
% 
% %pre-cancer stages
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %onset to first stage
% %lTEMPLATE{ R{onset,1}+I-1 ,R{onset,3}}=OnsetRates(I,1);
% 
% %Baseline Prevalence
% lTEMPLATE{ R{PCH1_1618_Prev0,1}+I-1 ,R{PCH1_1618_Prev0,3} }=PrevG(I,1);
% lTEMPLATE{ R{PCH2_1618_Prev0,1}+I-1 ,R{PCH2_1618_Prev0,3} }=PrevG(I,2);
% lTEMPLATE{ R{PCH3_1618_Prev0,1}+I-1 ,R{PCH3_1618_Prev0,3} }=PrevG(I,3);
% lTEMPLATE{ R{PCL1_Prev0,1}+I-1 ,R{PCL1_Prev0,3} }=PrevG(I,4);
% lTEMPLATE{ R{PCL2_Prev0,1}+I-1 ,R{PCL2_Prev0,3} }=PrevG(I,5);
% lTEMPLATE{ R{PCL3_Prev0,1}+I-1 ,R{PCL3_Prev0,3} }=PrevG(I,6);
% lTEMPLATE{ R{PCR1_Prev0,1}+I-1 ,R{PCR1_Prev0,3} }=PrevG(I,7);
% lTEMPLATE{ R{PCM1_Prev0,1}+I-1 ,R{PCM1_Prev0,3} }=PrevG(I,8);
% lTEMPLATE{ R{PCR2_Prev0,1}+I-1 ,R{PCR2_Prev0,3} }=PrevG(I,13);
% lTEMPLATE{ R{PCH1_H_Prev0,1}+I-1 ,R{PCH1_H_Prev0,3} }=PrevG(I,10);
% lTEMPLATE{ R{PCH2_H_Prev0,1}+I-1 ,R{PCH2_H_Prev0,3} }=PrevG(I,11);
% lTEMPLATE{ R{PCH3_H_Prev0,1}+I-1 ,R{PCH3_H_Prev0,3} }=PrevG(I,12);
% lTEMPLATE{ R{PCR3_Prev0,1}+I-1 ,R{PCR3_Prev0,3} }=PrevG(I,9);
% %CFR rate
% lTEMPLATE{ R{PCH1_1618_CFR,1}+I-1 ,R{PCH1_1618_CFR,3} }=0;
% lTEMPLATE{ R{PCH2_1618_CFR,1}+I-1 ,R{PCH2_1618_CFR,3} }=0;
% lTEMPLATE{ R{PCH3_1618_CFR,1}+I-1 ,R{PCH3_1618_CFR,3} }=0;
% lTEMPLATE{ R{PCH1_H_CFR,1}+I-1 ,R{PCH1_H_CFR,3} }=0;
% lTEMPLATE{ R{PCH2_H_CFR,1}+I-1 ,R{PCH2_H_CFR,3} }=0;
% lTEMPLATE{ R{PCH3_H_CFR,1}+I-1 ,R{PCH3_H_CFR,3} }=0;
% lTEMPLATE{ R{PCL1_CFR,1}+I-1 ,R{PCL1_CFR,3} }=0;
% lTEMPLATE{ R{PCL2_CFR,1}+I-1 ,R{PCL2_CFR,3} }=0;
% lTEMPLATE{ R{PCL3_CFR,1}+I-1 ,R{PCL3_CFR,3} }=0;
% lTEMPLATE{ R{PCR1_CFR,1}+I-1 ,R{PCR1_CFR,3} }=0;
% lTEMPLATE{ R{PCR2_CFR,1}+I-1 ,R{PCR2_CFR,3} }=0;
% lTEMPLATE{ R{PCR3_CFR,1}+I-1 ,R{PCR3_CFR,3} }=0;
% lTEMPLATE{ R{PCM1_CFR,1}+I-1 ,R{PCM1_CFR,3} }=0;
% 
% %Progression
% lTEMPLATE{ R{PCH1_1618_Progr,1}+I-1 ,R{PCH1_1618_Progr,3} }=0.0931;
% lTEMPLATE{ R{PCH2_1618_Progr,1}+I-1 ,R{PCH2_1618_Progr,3} }=0.2107;
% lTEMPLATE{ R{PCH1_H_Progr,1}+I-1 ,R{PCH1_H_Progr,3} }=0.0931;
% lTEMPLATE{ R{PCH2_H_Progr,1}+I-1 ,R{PCH2_H_Progr,3} }=0.2107;
% lTEMPLATE{ R{PCL1_Progr,1}+I-1 ,R{PCL1_Progr,3} }=0.0568;
% lTEMPLATE{ R{PCL2_Progr,1}+I-1 ,R{PCL2_Progr,3} }=0.0921;
% lTEMPLATE{ R{PCM1_Progr,1}+I-1 ,R{PCM1_Progr,3} }=1/Immunity;
% 
% if I<7
%     lTEMPLATE{ R{PCH3_1618_Progr,1}+I-1 ,R{PCH3_1618_Progr,3} }=0.0292;
%     lTEMPLATE{ R{PCH3_H_Progr,1}+I-1 ,R{PCH3_H_Progr,3} }=0.0292;
%     lTEMPLATE{ R{PCL3_Progr,1}+I-1 ,R{PCL3_Progr,3} }=0.0070;
% end
% if (I<8 && I>=7)
%     lTEMPLATE{ R{PCH3_1618_Progr,1}+I-1 ,R{PCH3_1618_Progr,3} }=0.0506;
%     lTEMPLATE{ R{PCH3_H_Progr,1}+I-1 ,R{PCH3_H_Progr,3} }=0.0506;
%     lTEMPLATE{ R{PCL3_Progr,1}+I-1 ,R{PCL3_Progr,3} }=0.0140;
% end
% if (I<9 && I>=8)
%     lTEMPLATE{ R{PCH3_1618_Progr,1}+I-1 ,R{PCH3_1618_Progr,3} }=0.1344;
%     lTEMPLATE{ R{PCH3_H_Progr,1}+I-1 ,R{PCH3_H_Progr,3} }=0.1344;
%     lTEMPLATE{ R{PCL3_Progr,1}+I-1 ,R{PCL3_Progr,3} }=0.0221;
% end
% if I>=9
%     lTEMPLATE{ R{PCH3_1618_Progr,1}+I-1 ,R{PCH3_1618_Progr,3} }=0.1952; 
%     lTEMPLATE{ R{PCH3_H_Progr,1}+I-1 ,R{PCH3_H_Progr,3} }=0.1952; 
%     lTEMPLATE{ R{PCL3_Progr,1}+I-1 ,R{PCL3_Progr,3} }=0.0445;
% end
% 
% %Regression
% lTEMPLATE{ R{PCH1_1618_Rem,1}+I-1 ,R{PCH1_1618_Rem,3} }=0.2693;
% lTEMPLATE{ R{PCH2_1618_Rem1,1}+I-1 ,R{PCH2_1618_Rem1,3} }=0.1188;
% lTEMPLATE{ R{PCH2_1618_Rem2,1}+I-1 ,R{PCH2_1618_Rem2,3} }=0.1188;
% lTEMPLATE{ R{PCH3_1618_Rem1,1}+I-1 ,R{PCH3_1618_Rem1,3} }=0.01715;
% lTEMPLATE{ R{PCH3_1618_Rem2,1}+I-1 ,R{PCH3_1618_Rem2,3} }=0.01715;
% 
% lTEMPLATE{ R{PCH1_H_Rem,1}+I-1 ,R{PCH1_H_Rem,3} }=0.2693;
% lTEMPLATE{ R{PCH2_H_Rem1,1}+I-1 ,R{PCH2_H_Rem1,3} }=0.1188;
% lTEMPLATE{ R{PCH2_H_Rem2,1}+I-1 ,R{PCH2_H_Rem2,3} }=0.1188;
% lTEMPLATE{ R{PCH3_H_Rem1,1}+I-1 ,R{PCH3_H_Rem1,3} }=0.01715;
% lTEMPLATE{ R{PCH3_H_Rem2,1}+I-1 ,R{PCH3_H_Rem2,3} }=0.01715;
% 
% 
% 
% lTEMPLATE{ R{PCL1_Rem,1}+I-1 ,R{PCL1_Rem,3} }=0.2693;
% lTEMPLATE{ R{PCL2_Rem1,1}+I-1 ,R{PCL2_Rem1,3} }=0.1059;
% lTEMPLATE{ R{PCL2_Rem2,1}+I-1 ,R{PCL2_Rem2,3} }=0.1059;
% lTEMPLATE{ R{PCL3_Rem1,1}+I-1 ,R{PCL3_Rem1,3} }=0.0704;
% lTEMPLATE{ R{PCL3_Rem2,1}+I-1 ,R{PCL3_Rem2,3} }=0.0704;
% lTEMPLATE{ R{PCR1_Rem,1}+I-1 ,R{PCR1_Rem,3} }=0.0363;
% lTEMPLATE{ R{PCR2_Rem,1}+I-1 ,R{PCR2_Rem,3} }=0.0363;
% lTEMPLATE{ R{PCR3_Rem,1}+I-1 ,R{PCR3_Rem,3} }=0.0363;
% end
% 
% end
% 
% end
% 
% function [lTEMPLATE2]=WriteTemplateM(RES2,TEMPLATEM)
% lTEMPLATE2=TEMPLATEM;
% PrevM=RES2{1};
% T=RES2{2};
% OnsetM=RES2{3};
% 
% OnsetR_1618=1;OnsetR_HR=2;OnsetR_LR=3;
% Prev_1618=4;Prev_HR=5;Prev_LR=6;Prev_Im=7;
% Reg_1618=8;Reg_HR=9;Reg_LR=10;Reg_Im=11;
% 
% %Onset rate of HPV
% R{OnsetR_1618,1}=9;R{OnsetR_1618,2}=20;R{OnsetR_1618,3}=7;
% R{OnsetR_HR,1}=9;R{OnsetR_HR,2}=20;R{OnsetR_HR,3}=10;
% R{OnsetR_LR,1}=9;R{OnsetR_LR,2}=20;R{OnsetR_LR,3}=13;
% 
% %Regression rate
% R{Reg_1618,1}=75;R{Reg_1618,2}=86;R{Reg_1618,3}=9;
% R{Reg_HR,1}=97;R{Reg_HR,2}=108;R{Reg_HR,3}=9;
% R{Reg_LR,1}=119;R{Reg_LR,2}=130;R{Reg_LR,3}=9;
% R{Reg_Im,1}=141;R{Reg_Im,2}=152;R{Reg_Im,3}=9;
% 
% %Prevalence
% P{Prev_1618,1}=75;P{Prev_1618,2}=86;P{Prev_1618,3}=4;
% P{Prev_HR,1}=97;P{Prev_HR,2}=108;P{Prev_HR,3}=4;
% P{Prev_LR,1}=119;P{Prev_LR,2}=130;P{Prev_LR,3}=4;
% P{Prev_Im,1}=141;P{Prev_Im,2}=152;P{Prev_Im,3}=4;
% 
% for I=1:12
% lTEMPLATE2{ R{OnsetR_1618,1}+I-1 ,R{OnsetR_1618,3}}=OnsetM(I,1);  
% lTEMPLATE2{ R{OnsetR_HR,1}+I-1 ,R{OnsetR_HR,3}}=OnsetM(I,2);  
% lTEMPLATE2{ R{OnsetR_LR,1}+I-1 ,R{OnsetR_LR,3}}=OnsetM(I,3);  
% 
% lTEMPLATE2{ R{Reg_1618,1}+I-1 ,R{Reg_1618,3}}=0.0363; 
% lTEMPLATE2{ R{Reg_HR,1}+I-1 ,R{Reg_HR,3}}=0.0363; 
% lTEMPLATE2{ R{Reg_LR,1}+I-1 ,R{Reg_LR,3}}=0.0363; 
% lTEMPLATE2{ R{Reg_Im,1}+I-1 ,R{Reg_Im,3}}=0.1; 
% 
% lTEMPLATE2{ P{Prev_1618,1}+I-1 ,P{Prev_1618,3}}=PrevM(I,1); 
% lTEMPLATE2{ P{Prev_HR,1}+I-1 ,P{Prev_HR,3}}=PrevM(I,2); 
% lTEMPLATE2{ P{Prev_LR,1}+I-1 ,P{Prev_LR,3}}=PrevM(I,3); 
% lTEMPLATE2{ P{Prev_Im,1}+I-1 ,P{Prev_Im,3}}=PrevM(I,4); 
% 
% end
% end
% 
