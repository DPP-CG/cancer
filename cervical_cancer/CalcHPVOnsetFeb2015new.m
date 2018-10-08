function [InsituOnsetRate,HPVonsetRate,PrevG,HPVonsetG,Immunity,Prate,beta0] = CalcHPVOnsetFeb2015new(InSituOnset, AgeArrayLower,AgeArrayUpper,PAR,Mark)

X =[1:101];
distHPV = PAR{10};
distHigh = PAR{12};

% Sexual debut age for country/region;

ageD=15;


if Mark==2
beta = ones(1,12);

HPVonsetG=InSituOnset;
FitArray=HPVonsetG;
[InsituOnsetRate, HPVonsetRate,PrevG,HPVonsetG,Immunity,Prate]=PopSim(beta,X,FitArray,AgeArrayLower,AgeArrayUpper,PAR,1,distHPV,ageD,distHigh,Mark);

else


beta0 = ones(1,12);
beta0 = beta0.*rand(1,12);
beta0(1,1:3)=0;
LB = [0 0 0 0 0 0 0 0 0 0 0 0];%zeros(1,size(AgeArrayLower,2)+1);
UB = [0.00001 0.00001 0.00001 1 1 1 1 1 1 1 1 1];%ones(1,size(AgeArrayLower,2)+1);



%OPTIONS = optimset('largescale','off','LevenbergMarquardt','on','display','off','Diagnostics','on','MaxFunEvals',700,'TolFun',1e-8,'TolCon',1e-8);
OPTIONS = optimset('MaxFunEvals',2000,'Diagnostics','on','TolFun',1e-6,'TolCon',1e-6);
Y =[InSituOnset]';
FitArray = [InSituOnset]' ;
[beta,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@(x,X)...
    PopSim(x,X,FitArray,AgeArrayLower,AgeArrayUpper,PAR,0, distHPV,ageD,distHigh,Mark),...
    beta0,X,Y,LB,UB,OPTIONS);

%plot the results
[InsituOnsetRate, HPVonsetRate,PrevG,HPVonsetG,Immunity,Prate]=PopSim(beta,X,FitArray,AgeArrayLower,AgeArrayUpper,PAR,1,distHPV,ageD,distHigh,Mark);

end
end

function [Result,HPVonsetRate,PrevG,HPVonsetG,Immunity,Prate] = PopSim(beta,X,InSituOnset,AgeArrayLower,AgeArrayUpper,PAR,plotInd, distHPV,ageD,distHigh,Mark)  

AgeGroupL=[1 5 10 15 (ageD+1) 25 30 40 50 60 70 80];
AgeGroupU=[4 9 11 ageD 24 29 39 49 59 69 79 100];

r2=0.0931; %HPV 16/18 to LSIL
r3=0.2107; %LSIL to HSIL
%% Modifying regression rates (r4, r5, r26, r27, r11, r12)

r5      = 0.01715; %0;%0.01715; %HSIL to Regression
r4      = 0.1188;  %0.064;%0.1188; %LSIL to HPV Regression
r11     = 0.0704;  %0.0233;%0.0704; %HSIL to Regression
r12     = 0.1059;  %0.0888;%0.1059; %LSIL to Regression
r26     = r4;
r27     = r5;

% refer function for more details
%[r4, r5, r11, r12, r26, r27] = get_modified_rate(r4, r5, r11, r12, r26, r27);

%%

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
pop0=100000;%work with pop of standard size

Stages = 15;
POP=zeros(m_a,Stages);
POP(1:m_a,1)=probAgeA*pop0;
POPM=POP(:,1:3);

% r15=0.1024;
% r16=0.0219;
r15 = zeros(size(POP,1),size(POP,1));%HSIL of HPV 16/18 to CIS, varied by ages
r16 = zeros(size(POP,1),size(POP,1));%HSIL of Lowrisk HPV to CIS, varied by ages
A = [0  0.0292 0.0506 0.1344 0.1952];
B = [0  0.0070 0.0140 0.0221 0.0445];

% A = [0  0.0292 0.0506 0.0506 0.0506];
% B = [0  0.0070 0.0140 0.0140 0.0140];

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
%
% for i=1:30
%     r15(i,i)=0.0292;
% end
% for i=31:39
%     r15(i,i)=0.0506;
% end
% for i=40:49
%     r15(i,i)=0.1344;
% end
% for i=50:size(POP,1)
%     r15(i,i)=0.1952;
% end
%
% for i=1:30
%     r16(i,i)=0.0070;
% end
% for i=31:39
%     r16(i,i)=0.0140;
% end
% for i=40:49
%     r16(i,i)=0.0221;
% end
% for i=50:size(POP,1)
%     r16(i,i)=0.0445;
% end

t_max=500;

% beta2=zeros(1,12);
% beta2(1,4:8)=beta;
% CFR = zeros(m_a,t_max);
% AgeArrayMid=((AgeArrayUpper+AgeArrayLower)/2);
%  HPVonsetRate = interp1(AgeArrayMid,beta2,[0:100],'PCHIP');
HPVonsetRate = zeros(1, m_a);
% for J = 15:100
%     HPVonsetRate(1,J)=(beta(1)*beta(2)*(beta(1)*(J-15))^(1/beta(2)))/(1+(beta(1)*(J-15))^(beta(2)));%log logsitic
%     %  HPVonsetRate(1,J)=beta(1)/J; %inverse
%     % %    %HPVonsetRate(1,J)=beta(1)/(beta(2)^J); %logistic
%     % %    %HPVonsetRate(1,J)=beta(1)*exp(-beta(2)*J);%exponetial
% end
% HPVonsetRate(1,1:14) = 0;
%  HPVonsetRate(1,61:end) = HPVonsetRate(1,60);
%
%  for J = 1:12
%      Ind = (AgeArrayLower(J)):(AgeArrayUpper(J));
%      HPVonsetRate(Ind)=beta(J);
%  end
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

% if plotInd == 1
%
%     modelfun = @(beta,J)(beta(1)*beta(2)*((beta(1).*J(:,1)).^(1/beta(2))))./(1+((beta(1).*J(:,1)).^(beta(2))));
%     beta0=[0.3 5];
%     mdl = NonLinearModel.fit(X(15:100,1),HPVonsetRate(1,15:100),modelfun,beta0)
%     HPVonsetRate(1,15:100)=feval(mdl,(X(15:100,1)));
% end
if Mark==1
    % HPVonsetRate(60:end)=HPVonsetRate(59);
    %r1 = eye(size(POP,1)).*HPVonsetRate * distHPV(1); %Normal to Highrisk HPV, to be calculated
    %r8 = eye(size(POP,1)).*HPVonsetRate * distHPV(2);%Normal to Lowrisk HPV, to be calculated
    %r1 = interp1(AgeArrayM,[4.68514E-15	4.68514E-15	5.92841E-15	2.11E-06	0.032319988	0.012890128	0.012089773	0.008998394	0.008716173	0.01079305	0.010052524	0.010387282],[0:100],'PCHIP');
    r1 = HPVonsetRate(1,:)*distHigh(1); %Normal to Highrisk HPV, to be calculated
    r8 = 0;%HPVonsetRate(2,:);%Normal to Lowrisk HPV, to be calculated
    %r23 = interp1(AgeArrayM,[1.75193E-14	1.75193E-14	2.21683E-14	7.89E-06	0.120855309	0.048200526	0.04520773	0.033648023	0.032592705	0.040358846	0.037589769	0.038841543],[0:100],'PCHIP');
    r23 = HPVonsetRate(1,:)*distHigh(2);
    %r8 = interp1(AgeArrayM,[2.22045E-14	2.22045E-14	2.87119E-14	1E-05	0.110608186	0.057244706	0.041869869	0.021836054	0.016636076	0.01802241	0.011777121	0.011742254],[0:100],'PCHIP');

end
if Mark==0
    r8 = HPVonsetRate(1,:)*distHPV(2);
    r1 = 0;
    r23 = 0;
end

if Mark==2
    HPVonsetG=InSituOnset;
    HPVonsetG2(:,1)=interp1(AgeArrayM,HPVonsetG(1:12,1),[0:100],'PCHIP');
    HPVonsetG2(:,2)=interp1(AgeArrayM,HPVonsetG(1:12,2),[0:100],'PCHIP');
    HPVonsetG2(:,3)=interp1(AgeArrayM,HPVonsetG(1:12,3),[0:100],'PCHIP');
    r1 =  HPVonsetG2(:,1)';
    r23 = HPVonsetG2(:,2)';
    r8 =  HPVonsetG2(:,3)';    
end

Immunity = 10;%beta(end); %Source of Immunity = 10 years is Insinga et al 2009,
r19=1./Immunity; %Immunity to Normal

dt=1/12;
v=zeros(101,t_max);
b=zeros(101,t_max);

for t=1:t_max
    
    %fast time loop
    value1618 = zeros(101,1);
    valueLow = zeros(101,1);
    valueHigh = zeros(101,1);
    v15 = zeros(101,1);
    v16 = zeros(101,1);
    v31 = zeros(101,1);
    for t1=1:1/dt
        
        %immunity%
        POP(:,9) = POP(:,9) + (-mu.*POP(:,9) -r19.*POP(:,9) + r7.*POP(:,3) + r6.*POP(:,4) + r14.*POP(:,6) + r13.*POP(:,7) + r20.*POP(:,10) + r17.*POP(:,2) + r18.*POP(:,5) + r22.*POP(:,15) + r21.*POP(:,11)+ r28.*POP(:,14) + r29.*POP(:,13) + r30.*POP(:,12))*dt;
        
        %HPV 16/18 Regression
        POP(:,10)=POP(:,10) + (-mu.*POP(:,10) + r4.*POP(:,3) + r5.*POP(:,4) - r20.*POP(:,10))*dt;
        %HPV Low Risk Regression
        POP(:,11)=POP(:,11) + (-mu.*POP(:,11) + r12.*POP(:,6) + r11.*POP(:,7) - r21.*POP(:,11))*dt;
        %HPV High Risk Regression
        POP(:,15)=POP(:,15) + (-mu.*POP(:,15) + r26.*POP(:,13) + r27.*POP(:,14) - r22.*POP(:,15))*dt;
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
        POP(:,4)=POP(:,4) + (-mu.*POP(:,4) -r5.*POP(:,4) - r6.*POP(:,4) - r15*POP(:,4) + r3.*POP(:,3))*dt;
        POP(:,3)=POP(:,3) + (-mu.*POP(:,3) -r7.*POP(:,3) - r4.*POP(:,3) - r3.*POP(:,3) + r2.*POP(:,2))*dt;
        POP(:,2)=POP(:,2) + (-mu.*POP(:,2) -r2.*POP(:,2) + r1'.*POP(:,1) - r17.*POP(:,2))*dt;
        
        POP(:,7)=POP(:,7) + (-mu.*POP(:,7) -r13.*POP(:,7) - r11.*POP(:,7) - r16*POP(:,7) + r10.*POP(:,6))*dt;
        POP(:,6)=POP(:,6) + (-mu.*POP(:,6) -r14.*POP(:,6) - r10.*POP(:,6) - r12.*POP(:,6) + r9.*POP(:,5))*dt;
        POP(:,5)=POP(:,5) + (-mu.*POP(:,5) -r9.*POP(:,5) + r8'.*POP(:,1) - r18.*POP(:,5))*dt;
        
        POP(:,14)=POP(:,14) + (-mu.*POP(:,14) -r28.*POP(:,14) - r27.*POP(:,14) - r31*POP(:,14) + r25.*POP(:,13))*dt;
        POP(:,13)=POP(:,13) + (-mu.*POP(:,13) -r29.*POP(:,13) - r25.*POP(:,13) - r26.*POP(:,13) + r24.*POP(:,12))*dt;
        POP(:,12)=POP(:,12) + (-mu.*POP(:,12) -r24.*POP(:,12) + r23'.*POP(:,1) - r30.*POP(:,12))*dt;
        
        POP(:,1)=POP(:,1) + (-mu.*POP(:,1) -r1'.*POP(:,1) - r8'.*POP(:,1) - r23'.*POP(:,1) + r19.*POP(:,9))*dt;
        
    end
    
    %      for m=1:10
    %      SUMPOP(t,m)=sum(POP(20:30,m));
    %      end
    value1618 = 1000* value1618 ./ sum(POP,2);
    valueLow = 1000*valueLow./ sum(POP,2);
    valueHigh = 1000*valueHigh./ sum(POP,2);
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
    
    Result = [InsituOnsetRate(1,1:101)];
    
end
dlmwrite('resultPrashant.csv',Result,'-append');

if plotInd==1
    
    PrevalenceHPVByAge = ((POP(:,2)+POP(:,5)+POP(:,12)+POP(:,10)+POP(:,11)+POP(:,15))./sum(POP,2));
    PrevalenceHPV = sum(PrevalenceHPVByAge (20:69,1).*probAgeA(1,20:69)')/sum(probAgeA(1,20:69));
    PrevG=zeros(12,13);
    AgePOP=sum(POP,2);
    %
    HPVonsetG=zeros(12,2);
    
    for i = 1 : 12
        
        indB=AgeArrayLower(i):AgeArrayUpper(i);
        
        for j = 1:13
            if j < 9
                if j == 7
                    PrevG(i,j)=sum(POP(indB,j+3))./sum(AgePOP(indB,1));
                else
                    PrevG(i,j)=sum(POP(indB,j+1))./sum(AgePOP(indB,1));
                end
            else
                PrevG(i,j)=sum(POP(indB,j+2))./sum(AgePOP(indB,1));
            end
        end
        
        TempCIS(i,1)=sum([InsituOnsetRate(1,indB)]');
        TCIS15(i,1)=sum([ISO1(1,indB)]');
        TCIS16(i,1)=sum([ISO2(1,indB)]');
        TempHPV(i,1)=PrevG(i,1)+PrevG(i,4)+PrevG(i,7);
        
        
        Prate(i,1)=sum(r15interp(indB)'.*AgePOP(indB,1))/sum(AgePOP(indB,1));
        Prate(i,2)=sum(r16interp(indB)'.*AgePOP(indB,1))/sum(AgePOP(indB,1));
        %          HPVonsetG(i,1)=sum(r1(1,indB)'.*sum(POP(indB,:),2))/sum(AgePOP(indB,1));
        %          HPVonsetG(i,2)=sum(r8(1,indB)'.*sum(POP(indB,:),2))/sum(AgePOP(indB,1));
    end
    
    if Mark==1
        HPVonsetG(1:12,1) = beta(1:12)*distHigh(1);
        HPVonsetG(1:12,2) = beta(1:12)*distHigh(2);
        HPVonsetG(1:12,3) = 0;
    end
    if Mark==0
        HPVonsetG(1:12,1) = 0;
        HPVonsetG(1:12,2) = 0;
        HPVonsetG(1:12,3) = beta(1:12);
    end
%  
%     Factor = ones(12,1)-sum(PrevG(1:12,:),2);
%     
%     HPVonset1618(1:12,1)=HPVonsetG(1:12,1).*Factor;
%     
%     PrevW=PrevG(:,1)+PrevG(:,2)+PrevG(:,3)+PrevG(:,7);
%     
%     
%     
%     beta1 = [0.5];
%     LB =[0];
%     UB = [1];
%     %LB(end) = 0;
%     %UB(end) = 0.00001;
%     % beta0 = [0 0 0];
%     % LB =[0 0 0];
%     % UB = [1 1 1];
%     X =[0:100]';
%     
%     %OPTIONS = optimset('largescale','off','LevenbergMarquardt','on','display','off','Diagnostics','on','MaxFunEvals',700,'TolFun',1e-8,'TolCon',1e-8);
%     OPTIONS = optimset('MaxFunEvals',2000,'Diagnostics','on','TolFun',1e-6,'TolCon',1e-6);
%     Y=[HPVonset1618(1:12,1)]' ;
%     FitArray=[HPVonset1618(1:12,1)];
%     [beta,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@(x,X)...
%         CalOnsetW(x,X,FitArray, t_max, AgeGroupL, AgeGroupU, mu, POPM, rc, r19, 0,Births,PrevW),...
%         beta1,X,Y,LB,UB,OPTIONS);
%     
%     [OnsetW] = CalOnsetW(beta, X, FitArray, t_max,  AgeGroupL, AgeGroupU, mu, POPM, rc, r19, 1,Births,PrevW);
%     %
%     %
%     %  %   pause
    % end
    
%     
%     figure
%     hold on
%     plot(1:12, FitArray,'b');
%     hold on
%     plot(1:12, OnsetW,'r');
%     title('HPV onset rates fit');
%     
%     
    xlswrite('POP.xlsx',POP);
    
    %     BetaValue=xlsread('BetaValue.xlsx','BetaValue');
    %     BetaValue(PAR{13},1)=beta(1);
    %     BetaValue(PAR{13},2)=beta(2);
    %     xlswrite('BetaValue.xlsx',BetaValue,'BetaValue');
    
    
    xlswrite('value.xlsx',v);
    
    
    fid = fopen('data.xls','a');
    fprintf(fid,'%6.8f \t',PrevalenceHPVByAge);
    fprintf(fid,'\n');
    fprintf(fid,'%6.8f \n',PrevalenceHPV);
    fprintf(fid,'%6.8f \t',HPVonsetRate);
    fprintf(fid,'\n');
    fprintf(fid,'%6.8f \t',InSituOnset);
    fprintf(fid,'\n');
    fprintf(fid,'%6.8f \t',Result);
    fprintf(fid,'\n');
    fclose(fid);
    
    
    
    close all
    figure
    hold on
    plot(0:100, PrevalenceHPVByAge(:,1),'r');
    title('HPV prevalence rates');
    
    figure
    hold on
    plot(0:100, HPVonsetRate(1,:),'r');
    title('HPV onset rates');
    
    if Mark~=2

    figure
    hold on
    plot(0:100, [InsituOnsetRate(1,:)]','r');
    plot(0:100, InSituOnset,'r*');
    legend('Estimated','Actual');
    title('Insitu rates');
    end
end
end

function [y1, y2, y3, y4, y5, y6] = get_modified_rate(x1, x2, x3, x4, x5, x6)

%% data from WHO/PAHO
% reference:screening guidelines since 2012
% 1. Luque, J. S., Maupin, J. N., Ferris, D. G., & Condorhuaman, W. S. G. (2016). Reaching women in the Peruvian Andes through cervical cancer screening campaigns: Assessing attitudes of stakeholders and patients. Patient preference and adherence, 10, 2107.
% 2. PDF file "Regression rate data": BOX\Vijeta\Activity log\Data

% s_c = 0.503;                    % screening coverage
% rate = [0.2; 1];                % first is HPV = 1/5 and second is pap = 1/1
% sensitivity = [0.88; 0.55];     % these values are directly from spectrum

%reference:screening guidelines from 2000 - 2012
%Who is getting Pap smears in urban Peru?
% Valerie A Paz Soldan,1* Frank H Lee,2 Cesar Carcamo,3 King K Holmes,4 Geoff P Garnett5
% and Patricia Garcia3
% International Journal of Epidemiology 2008;37:862–869
% doi:10.1093/ije/dyn118

 s_c = 0.309;                    % screening coverage for population age 18-29 years in 20 cities of Peru
 rate = [0; 1/3];                % first is HPV = 0 and second is pap = 1/3 (however there no data on frequency:propotion screened does significantly vary by number of children)
 sensitivity = [0.88; 0.55];     % these values are directly from spectrum
x = [x1, x2, x3, x4, x5, x6];


%% Calculations
% for fumula, refer "Effect of Screening on regression rate": BOX - Vijeta - Activity log - Manuscripts

%y = ((1 - s_c)*x) + ( ones(1,size(x,2)) * (s_c * (rate(2)*sensitivity(2) + rate(2)*(1 - sensitivity(2)) * rate(1)*sensitivity(1))));
y = ((1 - s_c)*x) + ( ones(1,size(x,2)) * (s_c * (rate'*sensitivity - prod(rate.*sensitivity))));
y1 = y(1); y2 = y(2); y3 = y(3); y4 = y(4); y5 = y(5); y6 = y(6);
end



