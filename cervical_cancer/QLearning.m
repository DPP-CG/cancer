function QLearning(PAR,steadyStatePopW, steadyStatePopM, RES2)
T=RES2{2};

% female mixing mat
M1=RES2{4};
% male mixing mat
M2=RES2{5};
numStates=20;%HPV prevalence from 0 to 10% in intervals of 0.5%
numActions = 20*(20);% Vaccine covergae at age 9 in 5% intervals * screening coverage at age 35.
numStages= 100;
Q=zeros(numStates, numActions, numStages);
V=zeros(numStates, numActions, numStages);
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

t_max=500;

rc=0.0363;
%PERU mixing matrix: women (age group in row ) mixing with men (age group in column)
%%%% Mixing matrix was updated on Feb 6th, 2018. 01:28pm.

OMT=zeros(12,3);
OM1618=zeros(101,1);
OMHR=zeros(101,1);
OMLR=zeros(101,1);

%DATA for women
AgeGroupL=[1 5 10 15 (ageD+1) 25 30 40 50 60 70 80];
AgeGroupU=[4 9 11 ageD 24 29 39 49 59 69 79 100];

r2=0.0931; %HPV 16/18 to LSIL
r3=0.2107; %LSIL to HSIL

r5      = 0.01715; %0;%0.01715; %HSIL to Regression
r4      = 0.1188;  %0.064;%0.1188; %LSIL to HPV Regression
r11     = 0.0704;  %0.0233;%0.0704; %HSIL to Regression
r12     = 0.1059;  %0.0888;%0.1059; %LSIL to Regression
r26     = r4;
r27     = r5;

r6= 0.01715;  %r5; %HSIL to Immunity
r7= 0.1188; %r4; %LSIL to Immunity

r9=0.0568; %Lowrisk HPV to LSIL
r10=0.0921; %LSIL to HSIL

r13= 0.0704; %r11; %HSIL to Immunity
r14= 0.1059; %r12; %LSIL to Immunity

r17=0.2693; %HPV 16/18 to Immunity
r18=0.2693; %LowRisk HPV to Immunity

%% Regression to Immunity

r20=0.2693; %0.0363;%0.2693; %Goldie et al 2003 TABLE I HPV clearance probabilities converted to rates
r21=r20;
r22=r20;

%%
%Highrisk HPV, parameter same as HPV 16/18
r24=r2;
r25=r3;
r28=r6;
r29=r7;
r30=r17;


rc=r20;
%InSituOnset(60:end) =InSituOnset(60);

probAgeA=PAR{1};probAgeA(81:end)=probAgeA(81)/length(probAgeA(81:end));
mu=PAR{2}';mu(81:end)=mu(80)+1/5;
BirthsF=PAR{3};%female pop only

AgeHPVprev = PAR{11};
% HPVprev = PAR{12};

agemax=101;
m_a=length(probAgeA);
%age,state, duration of disease
pop0=95000000;%work with pop of standard size

Stages = 15;
%POP=zeros(m_a,Stages);
%POP(1:m_a,1)=probAgeA*pop0;
POP=steadyStatePopW;
POP(:,16)=0;
POPM =steadyStatePopM;
%POPM=POP(:,1:3);


r15 = zeros(size(POP,1),size(POP,1));%HSIL of HPV 16/18 to CIS, varied by ages
r16 = zeros(size(POP,1),size(POP,1));%HSIL of Lowrisk HPV to CIS, varied by ages
A = [0  0.0292 0.0506 0.1344 0.1952];
B = [0  0.0070 0.0140 0.0221 0.0445];

Midage = [(1+14)/2 (15+30)/2 (31+39)/2 (40+49)/2 (50+101)/2];
r15interp = interp1(Midage,A,[0:100],'PCHIP');
r16interp = interp1(Midage,B,[0:100],'PCHIP');
for i=1:size(POP,1)
    r15(i,i)=r15interp(i);
end

for i=1:size(POP,1)
    r16(i,i)=r16interp(i);
end

r31=r15; %HSIL of HPV 16/18 to CIS, varied by ages, same as HPV 16/18
t_max=numStages;

HPVonsetRate = zeros(1, m_a);
AgeArrayM=(AgeGroupU+AgeGroupL)/2;

beta(5:12)=smooth(beta(5:12));
% beta(17:24)=smooth(beta(17:24));
HPVonsetRate(1,:)= interp1(AgeArrayM,beta(1:12),[0:100],'PCHIP');
% HPVonsetRate(2,:)= interp1(AgeArrayM,beta(13:24),[0:100],'PCHIP');
HPVonsetRate(1,1:15) = 0;
% HPVonsetRate(2,1:19) = 0;
for intI=1:101
    if HPVonsetRate(1,intI)<0
        HPVonsetRate(1,intI)=0;
    end
    %     if HPVonsetRate(2,intI)<0
    %         HPVonsetRate(2,intI)=0;
    %     end
end


Immunity = 10;%beta(end); %Source of Immunity = 10 years is Insinga et al 2009,
r19=1./Immunity; %Immunity to Normal

dt=1/12;
v=zeros(101,t_max);
b=zeros(101,t_max);



% Progression model
 %Onset Preclinical0 from HPV 16/18 CIN2_3
        value1618= (r15*POP(:,4) )*dt;
        %Onset Preclinical0 from HPV High Risk CIN2_3
        valueHigh= (r31*POP(:,14) )*dt;
        %Onset Preclinical0 from HPV Low Risk CIN2_3
        valueLow= ( r16*POP(:,7))*dt;
        
    value1618 =  value1618 ./ sum(POP,2);
    valueLow = valueLow./ sum(POP,2);
    valueHigh = valueHigh./ sum(POP,2);
    InsituOnsetRateHRLR(1,1:101) =value1618'+valueHigh';
    InsituOnsetRateHRLR(1,102:202)= valueLow';%1000*value ./ sum(POP,2);
    InsituOnsetRate(1,1:101) = InsituOnsetRateHRLR(1,1:101) + InsituOnsetRateHRLR(1,102:202);
    Result = sum([InsituOnsetRate(1,1:101).*probAgeA(1,1:101)]);
for tm=1:t_max
    %Determine Q-learning state
    fromState = ceil(Result*200)+1;
    %Action to take in stage
    vaccinceCov=zeros(1,101);
    screeningCov=zeros(1,101);
    randAction=randi(numActions);
    vaccinceCov(1,9)=ceil(randAction/ 20)*0.05;
    screeningCov(1,35)=mod(randAction, 20)*0.05;
    
    %Estimating onset rate of HPV in men
    C=0;%Vaccination coverage
    
    PrevG=zeros(12,15);
    AgePOP=sum(POP,2);
    
    for i = 1 : 12
        indB=AgeArrayLower(i):AgeArrayUpper(i);
        for j = 1:15
            PrevG(i,j)=sum(POP(indB,j))./sum(AgePOP(indB,1));
        end
    end
    PrevW(:,1)=PrevG(:,2)+PrevG(:,3)+PrevG(:,4)+PrevG(:,10);%HPV 16/18 prevalence
    PrevW(:,2)=PrevG(:,12)+PrevG(:,13)+PrevG(:,14)+PrevG(:,15);%HPV other HR
    PrevW(:,3)=PrevG(:,5)+PrevG(:,6)+PrevG(:,7)+PrevG(:,11);%HPV LR
    
    OMT(1:12,1)=T(:,2).*(M2*PrevW(:,1));
    OM1618(1:101,1)= interp1(AgeArrayMid,OMT(1:12,1),[1:101],'PCHIP');
    for intI=1:101
        if OM1618(intI,1)<0
            OM1618(intI,1)=0;
        end
    end
    
    OMT(1:12,2)=T(:,2).*(M2*PrevW(:,2));
    OMHR(1:101,1)= interp1(AgeArrayMid,OMT(1:12,2),[1:101],'PCHIP');
    for intI=1:101
        if OMHR(intI,1)<0
            OMHR(intI,1)=0;
        end
    end
    
    OMT(1:12,3)=T(:,2).*(M2*PrevW(:,3));
    OMLR(1:101,1)= interp1(AgeArrayMid,OMT(1:12,3),[1:101],'PCHIP');
    for intI=1:101
        if OMLR(intI,1)<0
            OMLR(intI,1)=0;
        end
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
    OnsetTW(1:12,1)=T(:,1).*(M1*PrevM(1:12,1));
    OnsetTW(1:12,2)=T(:,1).*(M1*PrevM(1:12,2));
    OnsetTW(1:12,3)=T(:,1).*(M1*PrevM(1:12,3));
    
    HPVonsetG2(:,1)=interp1(AgeArrayM,OnsetTW(1:12,1),[0:100],'PCHIP');
    HPVonsetG2(:,2)=interp1(AgeArrayM,OnsetTW(1:12,2),[0:100],'PCHIP');
    HPVonsetG2(:,3)=interp1(AgeArrayM,OnsetTW(1:12,3),[0:100],'PCHIP');
    r1 =  HPVonsetG2(:,1)';
    r23 = HPVonsetG2(:,2)';
    r8 =  HPVonsetG2(:,3)';
    
    
    
    % OnsetW=zeros(1,12);
    %HPV progression in men
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
    
    %%%HPV progression in women
    %fast time loop
    value1618 = zeros(101,1);
    valueLow = zeros(101,1);
    valueHigh = zeros(101,1);
    v15 = zeros(101,1);
    v16 = zeros(101,1);
    v31 = zeros(101,1);
    
    for t1=1:1/dt
        %Vaccinated
        POPt(:,16)=POP(:,16)+sum(POP(:,1:15)).* vaccinceCov(1,:);
        %immunity%
        POPt(:,9) = POP(:,9)-POP(:,9).* vaccinceCov(1,:) + (-mu.*POP(:,9) -r19.*POP(:,9) + r7.*POP(:,3) + r6.*POP(:,4) + r14.*POP(:,6) + r13.*POP(:,7) + r20.*POP(:,10) + r17.*POP(:,2) + r18.*POP(:,5) + r22.*POP(:,15) + r21.*POP(:,11)+ r28.*POP(:,14) + r29.*POP(:,13) + r30.*POP(:,12))*dt;
        
        %HPV 16/18 Regression
        POPt(:,10)=POP(:,10) -POP(:,10).* vaccinceCov(1,:) +(POP(:,3)+POP(:,4)).* screeningCov(1,:)+ (-mu.*POP(:,10) + r4.*POP(:,3) + r5.*POP(:,4) - r20.*POP(:,10))*dt;
        %HPV Low Risk Regression
        POPt(:,11)=POP(:,11) -POP(:,11).* vaccinceCov(1,:) + (POP(:,6)+POP(:,7)).* screeningCov(1,:)+(-mu.*POP(:,11) + r12.*POP(:,6) + r11.*POP(:,7) - r21.*POP(:,11))*dt;
        %HPV High Risk Regression
        POPt(:,15)=POP(:,15)-POP(:,15).* vaccinceCov(1,:)  + (POP(:,13)+POP(:,14)).* screeningCov(1,:)+(-mu.*POP(:,15) + r26.*POP(:,13) + r27.*POP(:,14) - r22.*POP(:,15))*dt;
        % POP(:,8)=POP(:,8) + (-mu.*POP(:,8) +r15*POP(:,4) + r16*POP(:,7))*dt;
        
        %Onset Preclinical0 from HPV 16/18 CIN2_3
        value1618=value1618 + (r15*POP(:,4) )*dt;
        %Onset Preclinical0 from HPV High Risk CIN2_3
        valueHigh=valueHigh + (r31*POP(:,14) )*dt;
        %Onset Preclinical0 from HPV Low Risk CIN2_3
        valueLow=valueLow + ( r16*POP(:,7))*dt;
        
        v15=v15+r15*POP(:,4)*dt;
        v31=v31+r31*POP(:,14)*dt;
        v16=v16+r16*POP(:,7)*dt;
        %CIN2/3-16/18
        POPt(:,4)=POP(:,4)-POP(:,4).* vaccinceCov(1,:) -POP(:,4).* screeningCov(1,:) + (-mu.*POP(:,4) -r5.*POP(:,4) - r6.*POP(:,4) - r15*POP(:,4) + r3.*POP(:,3))*dt;
        %CIN1-16/18
        POPt(:,3)=POP(:,3) -POP(:,3).* vaccinceCov(1,:)-POP(:,3).* screeningCov(1,:) + (-mu.*POP(:,3) -r7.*POP(:,3) - r4.*POP(:,3) - r3.*POP(:,3) + r2.*POP(:,2))*dt;
        %HPV16/18
        POPt(:,2)=POP(:,2) -POP(:,2).* vaccinceCov(1,:) + (-mu.*POP(:,2) -r2.*POP(:,2) + r1'.*POP(:,1) - r17.*POP(:,2))*dt;
        
        %CIN2/3-LR
        POPt(:,7)=POP(:,7) -POP(:,7).* vaccinceCov(1,:) -POP(:,7).* screeningCov(1,:)+ (-mu.*POP(:,7) -r13.*POP(:,7) - r11.*POP(:,7) - r16*POP(:,7) + r10.*POP(:,6))*dt;
        %CIN1-LR
        POPt(:,6)=POP(:,6) -POP(:,6).* vaccinceCov(1,:) -POP(:,6).* screeningCov(1,:)+ (-mu.*POP(:,6) -r14.*POP(:,6) - r10.*POP(:,6) - r12.*POP(:,6) + r9.*POP(:,5))*dt;
        %HPV LR
        POPt(:,5)=POP(:,5) -POP(:,5).* vaccinceCov(1,:) + (-mu.*POP(:,5) -r9.*POP(:,5) + r8'.*POP(:,1) - r18.*POP(:,5))*dt;
        
        %CIN2/3-HR
        POPt(:,14)=POP(:,14) -POP(:,14).* vaccinceCov(1,:)-POP(:,14).* screeningCov(1,:) + (-mu.*POP(:,14) -r28.*POP(:,14) - r27.*POP(:,14) - r31*POP(:,14) + r25.*POP(:,13))*dt;
        %CIN1-HR
        POPt(:,13)=POP(:,13) -POP(:,13).* vaccinceCov(1,:)-POP(:,13).* screeningCov(1,:)+ (-mu.*POP(:,13) -r29.*POP(:,13) - r25.*POP(:,13) - r26.*POP(:,13) + r24.*POP(:,12))*dt;
        %HPV HR
        POPt(:,12)=POP(:,12) -POP(:,12).* vaccinceCov(1,:) + (-mu.*POP(:,12) -r24.*POP(:,12) + r23'.*POP(:,1) - r30.*POP(:,12))*dt;
        
        %Susceptible
        POPt(:,1)=POP(:,1) -POP(:,1).* vaccinceCov(1,:) + (-mu.*POP(:,1) -r1'.*POP(:,1) - r8'.*POP(:,1) - r23'.*POP(:,1) + r19.*POP(:,9))*dt;
        
    end
    POP=POPt;
    %      for m=1:10
    %      SUMPOP(t,m)=sum(POP(20:30,m));
    %      end
    value1618 =  value1618 ./ sum(POP,2);
    valueLow = valueLow./ sum(POP,2);
    valueHigh = valueHigh./ sum(POP,2);
    InsituOnsetRateHRLR(1,1:101) =value1618'+valueHigh';
    InsituOnsetRateHRLR(1,102:202)= valueLow';%1000*value ./ sum(POP,2);
    InsituOnsetRate(1,1:101) = InsituOnsetRateHRLR(1,1:101) + InsituOnsetRateHRLR(1,102:202);
    ISO1(1,:)=v15./sum(POP,2);
    ISO2(1,:)=v16./sum(POP,2);
    % for z=1:101;
    v(:,t)=InsituOnsetRate(1,:);
    b(:,t)=POP(:,1);
    %  end
    %valueN = sum(sum(POP));
    %aging
    %lAge=POP(end,1)+POP(end,2)+POP(end,3)+POP(end,4);
    for I=1:Stages
        POP(2:end,I)=POP(1:end-1,I);
        POP(1,I)= 0;
    end
    %POP(end,4) = POP(end,4) +lAge;
    %must decide if births are to be based on population size
    Births=probAgeA(1,1)* sum(sum(POP));%42;%InsituOnsetRate.*probAgeA' *pop0;%BirthsF*sum(POP(:));
    POP(1,1)=Births;
    
    PrevalenceHPVByAge(1,:) =  (POP(:,2)+POP(:,5)+POP(:,12)+POP(:,10)+POP(:,11)+POP(:,15))./sum(POP,2);
    
    HPVageDist= zeros(1,size(AgeHPVprev,2));
    for J = 1:size(AgeHPVprev,2)
        Ind = AgeHPVprev(1,J):AgeHPVprev(2,J);
        HPVageDist(J)=(sum(PrevalenceHPVByAge(1,Ind).*probAgeA(1,Ind)))/sum(probAgeA(1,Ind));
        
    end
    %Result = [InsituOnsetRate(1,:), HPVageDist *10]';
    
    Result = sum([InsituOnsetRate(1,1:101).*probAgeA(1,1:101)]);
    dt = 1/12;
    V(fromState, randAction, tm)=V(fromState, randAction, tm)+1;
    alpha=0.1/V(fromState, randAction, tm);
   
    toState= ceil(Result*200)+1;
    Q(fromState, randAction, tm)=(1-alpha)*Q(from-state, randAction, tm)+alpha*(Result + max(Q(toState, :, tm)));
   
    
end

end

