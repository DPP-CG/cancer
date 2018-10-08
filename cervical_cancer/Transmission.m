
function [OnsetW,RES2] =Transmission(HPVonsetG,PrevW,PrevHPV,PAR)
AgeArrayLower=[1    5	10	15  20	25	30	40	50	60	70	80];
AgeArrayUpper=[4	9	14  19 	24	29	39	49	59	69	79	100];
AgeArrayMid=((AgeArrayUpper+AgeArrayLower)/2);

%exponential function
% beta1 =  ones(1,2);
% beta1 = beta1*0;
% LB =[0 0];
% UB = [400 10];



%step function
beta1 =  ones(1,12);
beta1 = beta1*0;
LB = [0 0 0 0 0 0 0 0 0 0 0 0];
UB = [0.00001 0.00001 0.00001 1 1 1 1 1 1 1 1 1];

%power function
%  beta1 =  ones(1,3);
%  beta1 = beta1*0;
%  LB = [0 0 1];
%  UB = [4000 100 1000];



X =AgeArrayMid;

%OPTIONS = optimset('largescale','off','LevenbergMarquardt','on','display','off','Diagnostics','on','MaxFunEvals',700,'TolFun',1e-8,'TolCon',1e-8);
OPTIONS = optimset('MaxFunEvals',2000,'Diagnostics','on','TolFun',1e-6,'TolCon',1e-6);
Y=[(HPVonsetG(1:12,1)+HPVonsetG(1:12,2));HPVonsetG(1:12,3)]*1000;
FitArray=[(HPVonsetG(1:12,1)+HPVonsetG(1:12,2));HPVonsetG(1:12,3)]*1000;
[beta,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@(x,X)...
    CalOnsetW(x,X,FitArray,0,PrevW,PAR),...
    beta1,X,Y,LB,UB,OPTIONS);


[OnsetW,PrevM,T,OMT,M1,M2] = CalOnsetW(beta,X,FitArray,1,PrevW,PAR);
T2=T;
T2(5:12,1) =smooth(T(5:12,1) );
T2(5:12,2) =smooth(T(5:12,2) );

figure
hold on
plot(1:12,OnsetW(1:12),'r-');
plot(1:12,FitArray(1:12),'k*');
title('Onset Rate for HR');

figure
hold on
plot(1:12,OnsetW(13:24),'r-');
plot(1:12,FitArray(13:24),'k*');
title('Onset Rate for LR');

% figure
% hold on
% plot(1:12,T(1,1:12),'r-');
% plot(1:12,T(2,1:12),'k*');
% title('Onset Rate for LR');

RES2{1}=PrevM;
RES2{2}=T2;
RES2{3}=OMT;
% female mixing mat
RES2{4}=M1;
% male mixing mat
RES2{5}=M2;
end



%Compartmental model for men of HPV infection
function [OnsetW,PrevM,T,OMT,M1,M2] =CalOnsetW(beta1,X,HPVonsetG,plotInd2,PrevW,PAR)
AgeArrayLower=[1    5	10	15  20	25	30	40	50	60	70	80];
AgeArrayUpper=[4	9	14  19 	24	29	39	49	59	69	79	100];
AgeArrayMid=((AgeArrayUpper+AgeArrayLower)/2);
p=1;
rm=1/10; %immunity

% Population proportion by individual age
probAgeA=PAR{1};probAgeA(81:end)=probAgeA(81)/length(probAgeA(81:end));
probAgeB=PAR{15};probAgeB(81:end)=probAgeB(81)/length(probAgeB(81:end));

%Natural mortality by individual age
mu=PAR{16}';mu(81:end)=mu(80)+1/5;

%Population proportion by age group
for i = 1 : 12
    indC=AgeArrayLower(i):AgeArrayUpper(i);
    POP_Dist(1,i)=sum(probAgeA(1,indC));
    POP_Dist(2,i)=sum(probAgeB(1,indC));
end

m_a=length(probAgeB);
pop0=100000;
Stages = 5;
POPM=zeros(m_a,Stages);
POPM(1:m_a,1)=probAgeB*pop0;
t_max=500;

rc=0.0363;
%PERU mixing matrix: women (age group in row ) mixing with men (age group in column)
%%%% Mixing matrix was updated on Feb 6th, 2018. 01:28pm.
if  PAR{13} == 0  
    M1=[
        0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0
0	0	0.000000	0.380738	0.432729	0.180930	0.005601	0.000002	0.000000	0.000000	0.000000	0.000000
0	0	0.000000	0.055943	0.276534	0.502872	0.163553	0.001098	0.000000	0.000000	0.000000	0.000000
0	0	0.000000	0.004632	0.064654	0.331974	0.568303	0.030406	0.000030	0.000000	0.000000	0.000000
0	0	0.000000	0.000007	0.000558	0.015307	0.382802	0.584946	0.016371	0.000008	0.000000	0.000000
0	0	0.000000	0.000000	0.000000	0.000007	0.004320	0.376646	0.601428	0.017590	0.000009	0.000000
0	0	0.000000	0.000000	0.000000	0.000000	0.000001	0.004320	0.376652	0.601437	0.017590	0.000000
0	0	0.000000	0.000000	0.000000	0.000000	0.000000	0.000001	0.004395	0.383134	0.611788	0.000683
0	0	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000002	0.007847	0.684134	0.308017
0	0	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000274	0.999726
];
M2=[
    0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0.57359186	0.34790105	0.07762722	0.00087978	0.00000009	0.00000000	0.00000000	0.00000000	0.00000000
0	0	0	0.46072364	0.40536926	0.13120991	0.00269662	0.00000057	0.00000000	0.00000000	0.00000000	0.00000000
0	0	0	0.29794762	0.44537548	0.24491659	0.01175311	0.00000721	0.00000000	0.00000000	0.00000000	0.00000000
0	0	0	0.02639341	0.18812648	0.49329853	0.28815961	0.00402094	0.00000103	0.00000000	0.00000000	0.00000000
0	0	0	0.00005150	0.00251402	0.04514524	0.57288493	0.37490958	0.00449374	0.00000099	0.00000000	0.00000000
0	0	0	0.00000000	0.00000045	0.00005852	0.01758881	0.60140148	0.37662980	0.00432003	0.00000091	0.00000000
0	0	0	0.00000000	0.00000000	0.00000000	0.00000942	0.01758969	0.60143183	0.37664880	0.00432025	0.00000000
0	0	0	0.00000000	0.00000000	0.00000000	0.00000000	0.00000946	0.01766418	0.60397881	0.37824386	0.00010368
0	0	0	0.00000000	0.00000000	0.00000000	0.00000000	0.00000000	0.00000009	0.00122729	0.31007246	0.68870017
    ];
% M1=[
%     0	0	0	0	0	0	0	0	0	0	0	0
% 0	0	0	0	0	0	0	0	0	0	0	0
% 0	0	0	0	0	0	0	0	0	0	0	0
% 0	0	0.000000	0.304379	0.479302	0.257838	0.011945	0.000007	0.000000	0.000000	0.000000	0.000000
% 0	0	0.000000	0.082461	0.327776	0.479302	0.109984	0.000477	0.000000	0.000000	0.000000	0.000000
% 0	0	0.000000	0.009512	0.102774	0.408517	0.464301	0.014888	0.000009	0.000000	0.000000	0.000000
% 0	0	0.000000	0.000035	0.001848	0.036382	0.534902	0.420769	0.006062	0.000002	0.000000	0.000000
% 0	0	0.000000	0.000000	0.000000	0.000035	0.012784	0.549054	0.431902	0.006223	0.000002	0.000000
% 0	0	0.000000	0.000000	0.000000	0.000000	0.000005	0.012784	0.549072	0.431915	0.006223	0.000000
% 0	0	0.000000	0.000000	0.000000	0.000000	0.000000	0.000005	0.012862	0.552417	0.434547	0.000168
% 0	0	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000008	0.019749	0.848201	0.132042
% 0	0	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000792	0.999208
% 
%    ];
% %mixing matrix: men (age group in row ) mixing with women (age group in column)
% M2=[
%     0	0	0	0	0	0	0	0	0	0	0	0
% 0	0	0	0	0	0	0	0	0	0	0	0
% 0	0	0	0	0	0	0	0	0	0	0	0
% 0	0	0	0.57359186	0.34790105	0.07762722	0.00087978	0.00000009	0.00000000	0.00000000	0.00000000	0.00000000
% 0	0	0	0.53813518	0.36800988	0.09258339	0.00127139	0.00000016	0.00000000	0.00000000	0.00000000	0.00000000
% 0	0	0	0.23715084	0.44084634	0.30147780	0.02050557	0.00001946	0.00000000	0.00000000	0.00000000	0.00000000
% 0	0	0	0.01439399	0.13253027	0.44890491	0.39497051	0.00919639	0.00000392	0.00000000	0.00000000	0.00000000
% 0	0	0	0.00001158	0.00078800	0.01972204	0.42570077	0.54117174	0.01260050	0.00000537	0.00000000	0.00000000
% 0	0	0	0.00000000	0.00000006	0.00001175	0.00622279	0.43191019	0.54906546	0.01278430	0.00000545	0.00000000
% 0	0	0	0.00000000	0.00000000	0.00000000	0.00000164	0.00622288	0.43191693	0.54907403	0.01278450	0.00000003
% 0	0	0	0.00000000	0.00000000	0.00000000	0.00000000	0.00000166	0.00630069	0.43731747	0.55593946	0.00044072
% 0	0	0	0.00000000	0.00000000	0.00000000	0.00000000	0.00000000	0.00000001	0.00026258	0.13466768	0.86506973
% 
%     ];
elseif  PAR{13} == 2 || PAR{13} ==12 || PAR{13} ==17 || PAR{13} ==16

% %%%% mixing matrix: AFRE: women (age group in row ) mixing with men (age group in column)
M1=[
    0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0.29292	0.53267	0.17324	0.00116	0	0	0	0
0	0	0	0	0.06496	0.33352	0.57095	0.03055	0	0	0	0
0	0	0	0	0.00056	0.01531	0.38281	0.58495	0.01637	0	0	0
0	0	0	0	0	0	0.00432	0.37665	0.60143	0.01759	0	0
0	0	0	0	0	0	0	0.00432	0.37665	0.60144	0.01759	0
0	0	0	0	0	0	0	0	0.00439	0.38313	0.61179	0.00068
0	0	0	0	0	0	0	0	0	0.00785	0.68413	0.30802
0	0	0	0	0	0	0	0	0	0	0.00027	0.99973

   ];
%mixing matrix: men (age group in row ) mixing with women (age group in column)
M2=[
    0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0.46072	0.40537	0.13121	0.0027	0	0	0	0	0
0	0	0	0.29795	0.44538	0.24492	0.01175	0	0	0	0	0
0	0	0	0.02639	0.18813	0.4933	0.28816	0.00402	0	0	0	0
0	0	0	0	0.00251	0.04515	0.57288	0.37491	0.00449	0	0	0
0	0	0	0	0	0	0.01759	0.6014	0.37663	0.00432	0	0
0	0	0	0	0	0	0	0.01759	0.60143	0.37665	0.00432	0
0	0	0	0	0	0	0	0	0.01766	0.60398	0.37824	0.0001
0	0	0	0	0	0	0	0	0	0.00123	0.31007	0.6887
  ];
else
    'No partnership mixing data for this country'
end
 T=zeros(2,12);
 Tbyage=zeros(1,100);
 %estimating turnover rate in men (exponential function) 
% for J = 20:100
%     Tbyage(1,J)=beta1(1)*exp(-beta1(2)*J);
% end
% Tbyage(1,1:19) = 0;
% 
% for i = 1 : 12
%     indC=AgeArrayLower(i):AgeArrayUpper(i);
%     T(2,i)=sum(Tbyage(1,indC).*(probAgeB(1,indC)/sum(probAgeB(1,indC))));
% end

%estimating turnover rate in men (step function)
T(2,:)= beta1(1,1:12);

%estimating  turnover rate in men (power function)
% for J = 20:100
%     Tbyage(1,J)=beta1(1)*beta1(3)^(-beta1(2)*J);
% end
% Tbyage(1,1:19) = 0;
% for i = 1 : 12
%     indC=AgeArrayLower(i):AgeArrayUpper(i);
%     T(2,i)=sum(Tbyage(1,indC).*(probAgeB(1,indC)/sum(probAgeB(1,indC))));
% end

%Calculating  turnover rate in women as function of T in men 
for j=1:12
T(1,j)=(POP_Dist(2,:).*T(2,:))*M2(:,j);
end
T(1,:)= T(1,:)./POP_Dist(1,:);
T=T';
%  T=[
% 0	0
% 0	0
% 0	0
% 0.392298347	0
% 0.611470608	0.62550531
% 0.656352533	0.582196892
% 0.32928124	0.525204351
% 0.326174194	0.455746013
% 0.259652763	0.395021842
% 0.190562783	0.342386354
% 0.111967601	0.298280157
% 0.121110378	0.242191862];

OMT=zeros(12,3);
OM1618=zeros(101,1);
OMHR=zeros(101,1);
OMLR=zeros(101,1);
%Estimating onset rate of HPV in men
C=0;%Vaccination coverage
OMT(1:12,1)=p*T(:,2).*(M2*PrevW(:,1));
OM1618(1:101,1)= interp1(AgeArrayMid,OMT(1:12,1),[1:101],'PCHIP');
for intI=1:101
    if OM1618(intI,1)<0
        OM1618(intI,1)=0;
    end
end

OMT(1:12,2)=p*T(:,2).*(M2*PrevW(:,2));
OMHR(1:101,1)= interp1(AgeArrayMid,OMT(1:12,2),[1:101],'PCHIP');
for intI=1:101
    if OMHR(intI,1)<0
        OMHR(intI,1)=0;
    end
end

OMT(1:12,3)=p*T(:,2).*(M2*PrevW(:,3));
OMLR(1:101,1)= interp1(AgeArrayMid,OMT(1:12,3),[1:101],'PCHIP');
for intI=1:101
    if OMLR(intI,1)<0
        OMLR(intI,1)=0;
    end
end



dt = 1/12;
alpha=1;

% Progression model in men 
for tm=1:t_max
    % OnsetW=zeros(1,12);
    for t1=1:1/dt
        %immunity
        POPM(:,3)=POPM(:,3) + (-mu.*POPM(:,3) + rc.*POPM(:,2) + rc.*POPM(:,4) + rc.*POPM(:,5) - rm.*POPM(:,3))*dt;
        %HPV 16/18
        POPM(:,2)=POPM(:,2) + (-mu.*POPM(:,2) + OM1618.*POPM(:,1) - rc.*POPM(:,2))*dt;
        %HPV High Risk
        POPM(:,4)=POPM(:,4) + (-mu.*POPM(:,4) + OMHR.*POPM(:,1) - rc.*POPM(:,4))*dt;
        %HPV Low Risk
        POPM(:,5)=POPM(:,5) + (-mu.*POPM(:,5) + OMLR.*POPM(:,1) - rc.*POPM(:,5))*dt;
        %HPV Disease Free
        POPM(:,1)=POPM(:,1) + (-mu.*POPM(:,1) - OM1618.*POPM(:,1) - OMHR.*POPM(:,1) - OMLR.*POPM(:,1) + rm.*POPM(:,3))*dt;
    end
    for I=1:5
        POPM(2:end,I)=POPM(1:end-1,I);
        POPM(1,I)= 0;
    end
    Births=probAgeB(1,1)* sum(sum(POPM));%42;%InsituOnsetRate.*probAgeA' *pop0;%BirthsF*sum(POP(:));
    POPM(1,1)=Births;
end

AgePOPM=sum(POPM,2);
%
for i = 1 : 12
    indC=AgeArrayLower(i):AgeArrayUpper(i);
    PrevM(i,1)=sum(POPM(indC,2))./sum(AgePOPM(indC,1));
    PrevM(i,2)=sum(POPM(indC,4))./sum(AgePOPM(indC,1));
    PrevM(i,3)=sum(POPM(indC,5))./sum(AgePOPM(indC,1));
    PrevM(i,4)=sum(POPM(indC,3))./sum(AgePOPM(indC,1));
end
%Calculate onset rate in women
OnsetTW(1:12,1)=alpha*p*T(:,1).*(M1*PrevM(1:12,1))*1000;
OnsetTW(1:12,2)=alpha*p*T(:,1).*(M1*PrevM(1:12,2))*1000;
OnsetTW(1:12,3)=alpha*p*T(:,1).*(M1*PrevM(1:12,3))*1000;

OnsetW(1:12,1)=OnsetTW(1:12,1)+OnsetTW(1:12,2);
OnsetW(13:24,1)=OnsetTW(1:12,3);
end