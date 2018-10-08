%updated March 2018 to incoporate following hanges:  mortality rates for Tx and NoTx are now available(see Gopalappa, MDM, 2018); Tx coverage can be assumed to be 100%  
function [Mu_Tx, Impact, BaseMortality, TxCoverage, StgDist, RMortalityTx, RMortalityNoTx, ModDwellTime,A,B,C,D,E,F,G] = Mortality(DwellTime, MortDFbyAge, Prevalence, StageDist, CancerMort,PAR)
%Okonkwo et. al, JNCI 2008		
RNoTx = [0.079, 0.136, 0.231, 0.500];
%Zelle et. al (original Groot et. al., 2006)	
RTx = [0.006, .020666667, 0.084, .27];

%Below Commented by Prashant on 10 Jul 2016. 
% %Groot, The Breast Journal, Volume 12 Suppl. 1, 2006; For African region
% value = mean( [0.042, 0.093]);
% RTx = [0.006, 0.006, value, 0.275];
% value = mean([0.063,.15]);
% RNoTx = [0.02, 0.02, value, 0.3];
MaxAge = 101;
RTxByNoTx = RTx./RNoTx;
for i = 1:MaxAge
StgDist(i,:) =  StageDist;
end

ModDwellTime = DwellTime;

RMortalityTx = zeros(4,MaxAge);
RMortalityNoTx = zeros(4,MaxAge);
 t = 10;
for i = 1:MaxAge - 10
    
   DwellRate(1,i:i+9)= 1 ./ DwellTime(1,i:i+9);
   DwellRate(2,i:i+9)= 1 ./ DwellTime(2,i:i+9);
   DwellRate(3,i:i+9)= 1 ./ DwellTime(3,i:i+9);
   DwellRate(4,i:i+9)= 1 ./ DwellTime(4,i:i+9);
  
% Ten-year survival without treatment
   Healthy(1,i) = TenYearSurvival(zeros(1,10),MortDFbyAge(1,i:i+9));
   Insitu(1,i) = TenYearSurvival(DwellRate(1:4,i:i+9),MortDFbyAge(1,i:i+9));
   Local(1,i) = TenYearSurvival(DwellRate(2:4,i:i+9),MortDFbyAge(1,i:i+9));
   Regional(1,i) = TenYearSurvival(DwellRate(3:4,i:i+9),MortDFbyAge(1,i:i+9));
   Distant(1,i) = TenYearSurvival(DwellRate(4:4,i:i+9),MortDFbyAge(1,i:i+9));
   
   
% Ten-year RELATIVE survival without treatment   
   RRInsitu(1,i) = Insitu(1,i) / Healthy(1,i);
   RRLocal(1,i) = Local(1,i) / Healthy(1,i);
   RRRegional(1,i) = Regional(1,i) / Healthy(1,i);
   RRDistant(1,i) = Distant(1,i) / Healthy(1,i);
  
  %Yearly mortality rate without treatment = mortality from disease + natural mortality(Pr(mortality > RRInsitu) = exp^-rate 
   Mu_NoTx(1,i) = (MortDFbyAge(1,i)*t - log(RRInsitu(1,i)))/t;
   Mu_NoTx(2,i) = (MortDFbyAge(1,i)*t - log(RRLocal(1,i)))/t;
   Mu_NoTx(3,i) = (MortDFbyAge(1,i)*t - log(RRRegional(1,i)))/t;
   Mu_NoTx(4,i) = (MortDFbyAge(1,i)*t - log(RRDistant(1,i)))/t;
     
  
end

Mu_NoTx(1,91:101) = interp1(1:90,Mu_NoTx(1,1:90),(91:101),'cubic','extrap');
Mu_NoTx(2,91:101) = interp1(1:90,Mu_NoTx(2,1:90),(91:101),'cubic','extrap');
Mu_NoTx(3,91:101) = interp1(1:90,Mu_NoTx(3,1:90),(91:101),'cubic','extrap');
Mu_NoTx(4,91:101) = interp1(1:90,Mu_NoTx(4,1:90),(91:101),'cubic','extrap');

for i = 1:MaxAge
  Mu_Tx(:,i) = max(0,(Mu_NoTx(:,i) - MortDFbyAge(1,i))) .* RTxByNoTx(:) + MortDFbyAge(1,i);
end

 Mu_Tx(1,:) = max(0,Mu_Tx(1,:) - MortDFbyAge(1,:));
    Mu_Tx(2,:) = max(0,Mu_Tx(2,:) - MortDFbyAge(1,:));
    Mu_Tx(3,:) = max(0,Mu_Tx(3,:) - MortDFbyAge(1,:));
    Mu_Tx(4,:) = max(0,Mu_Tx(4,:) - MortDFbyAge(1,:));
    
    Mu_NoTx(1,:) = max(0,Mu_NoTx(1,:) - MortDFbyAge(1,:));
    Mu_NoTx(2,:) = max(0,Mu_NoTx(2,:) - MortDFbyAge(1,:));
    Mu_NoTx(3,:) = max(0,Mu_NoTx(3,:) - MortDFbyAge(1,:));
    Mu_NoTx(4,:) = max(0,Mu_NoTx(4,:) - MortDFbyAge(1,:));
    
   
    RMortalityTx(1,:)= Mu_Tx(1,:)./ MortDFbyAge(1,:);
    RMortalityTx(2,:)= Mu_Tx(2,:)./ MortDFbyAge(1,:);
    RMortalityTx(3,:)= Mu_Tx(3,:)./ MortDFbyAge(1,:);
    RMortalityTx(4,:)= Mu_Tx(4,:)./ MortDFbyAge(1,:);
    
    RMortalityNoTx(1,:)= Mu_NoTx(1,:)./ MortDFbyAge(1,:);
    RMortalityNoTx(2,:)= Mu_NoTx(2,:)./ MortDFbyAge(1,:);
    RMortalityNoTx(3,:)= Mu_NoTx(3,:)./ MortDFbyAge(1,:);
    RMortalityNoTx(4,:)= Mu_NoTx(4,:)./ MortDFbyAge(1,:); 
Coverage = ones(MaxAge, 4)*1;
BaseMortality = zeros(size(Mu_Tx,1), size(Mu_Tx,2));
BaseMortality(1,:) = (Coverage(:,1)' .*  Mu_Tx(1,:) + (1- Coverage(:,1)').* Mu_NoTx(1,:))  ;
BaseMortality(2,:) = (Coverage(:,2)' .*  Mu_Tx(2,:) + (1- Coverage(:,2)').* Mu_NoTx(2,:))  ;
BaseMortality(3,:) = (Coverage(:,3)' .*  Mu_Tx(3,:) + (1- Coverage(:,3)').* Mu_NoTx(3,:))  ;
BaseMortality(4,:) = (Coverage(:,4)' .*  Mu_Tx(4,:) + (1- Coverage(:,4)').* Mu_NoTx(4,:)) ;

% estimating increase in weighted mortality for unit change in coverage
Impact (1,:) = (-Mu_Tx(1,:) + Mu_NoTx(1,:)) ./ BaseMortality(1,:);
Impact (2,:) = (-Mu_Tx(2,:) + Mu_NoTx(2,:)) ./ BaseMortality(2,:);
Impact (3,:) = (-Mu_Tx(3,:) + Mu_NoTx(3,:)) ./ BaseMortality(3,:);
Impact (4,:) = (-Mu_Tx(4,:) + Mu_NoTx(4,:)) ./ BaseMortality(4,:);

TxCoverage(:,:) = max(0,Coverage(:,1:4));

calcCancerMort =  StageDist(1,1) * BaseMortality(1,:);
for i = 2: 4
calcCancerMort = calcCancerMort + StageDist(1,i) *BaseMortality(i,:);
end

close all
figure
hold on
plot(1:MaxAge-10, Healthy(1,1:MaxAge-10), 'r');
plot(1:MaxAge-10, Insitu(1,1:MaxAge-10), 'y');
plot(1:MaxAge-10, Local(1,1:MaxAge-10), 'b');
plot(1:MaxAge-10, Regional(1,1:MaxAge-10), 'g');
plot(1:MaxAge-10, Distant(1,1:MaxAge-10), 'c');
legend('healthy', 'in-situ', 'local', 'regional', 'distant');
% plot(30:80, T_Local(1,30:80), 'b-');
% plot(30:80, T_Regional(1,30:80), 'g-');
% plot(30:80, T_Distant(1,30:80), 'c-');
title('10-year Survival')

% % % for row=1:length(Healthy)
% % %     for col=1:MaxAge-10
% % %         if row == 1
% % %             OUTPUT_DATA{row,col} = {row,col};
% % %         else
% % %             OUTPUT_DATA{row,col} = [row,Healthy(col)];
% % %         end
% % %         
% % %     end
% % % end


% % % %title('10-year Survival')
% % % for i = 1 : MaxAge - 10
% % %     OUTPUT_DATA{i,1}=i;
% % %     OUTPUT_DATA{i,2}=Healthy(1,i);
% % %     OUTPUT_DATA{i,3}=Insitu(1,i);
% % %     OUTPUT_DATA{i,4}=Local(1,i);
% % %     OUTPUT_DATA{i,5}=Regional(1,i);
% % %     OUTPUT_DATA{i,6}=Distant(1,i);    
% % % end
% % % 
% % % for i = MaxAge - 10 : MaxAge
% % %     OUTPUT_DATA{i,1}=i;
% % % end

% OUTPUT_DATA

figure
hold on
plot(1:MaxAge, RMortalityNoTx(1,1:MaxAge), 'y');
plot(1:MaxAge, RMortalityNoTx(2,1:MaxAge), 'b');
plot(1:MaxAge, RMortalityNoTx(3,1:MaxAge), 'g');
plot(1:MaxAge,RMortalityNoTx(4,1:MaxAge), 'c');

legend('in-situ', 'local', 'regional', 'distant');
plot(1:MaxAge, RMortalityTx(1,1:MaxAge), 'y+');
plot(1:MaxAge, RMortalityTx(2,1:MaxAge), 'b+');
plot(1:MaxAge, RMortalityTx(3,1:MaxAge), 'g+');
plot(1:MaxAge, RMortalityTx(4,1:MaxAge), 'c+');

title('yearly relative mortality in cancer comaed to non-cancer persons ("-"without Tx "+" with Tx)')

% % % for i = 1 : MaxAge
% % %     
% % %     %legend('in-situ', 'local', 'regional', 'distant');
% % %     OUTPUT_DATA{i,8}=RMortalityNoTx(1,i);
% % %     OUTPUT_DATA{i,9}=RMortalityNoTx(2,i);
% % %     OUTPUT_DATA{i,10}=RMortalityNoTx(3,i);
% % %     OUTPUT_DATA{i,11}=RMortalityNoTx(4,i); 
% % %     
% % %     OUTPUT_DATA{i,13}=RMortalityTx(1,i);
% % %     OUTPUT_DATA{i,14}=RMortalityTx(2,i);
% % %     OUTPUT_DATA{i,15}=RMortalityTx(3,i);
% % %     OUTPUT_DATA{i,16}=RMortalityTx(4,i); 
% % %     
% % % end

%OUTPUT_DATA

% figure
% hold on
% plot(1:MaxAge-10, RRInsitu(1,1:MaxAge-10), 'y');
% plot(1:MaxAge-10, RRLocal(1,1:MaxAge-10), 'b');
% plot(1:MaxAge-10, RRRegional(1,1:MaxAge-10), 'g');
% plot(1:MaxAge-10,RRDistant(1,1:MaxAge-10), 'c');
% legend('in-situ', 'local', 'regional', 'distant');
% plot(1:MaxAge-10, RI(1,1:MaxAge-10), 'y+');
% plot(1:MaxAge-10, RL(1,1:MaxAge-10), 'b+');
% plot(1:MaxAge-10, RR(1,1:MaxAge-10), 'g+');
% plot(1:MaxAge-10, RD(1,1:MaxAge-10), 'c+');
% title('10-year relative survival in cancer compared to non-cancer persons ("-"without Tx "+" with Tx)')




% % % for i = 1 : MaxAge
% % %     
% % % %   title('10-year relative survival in cancer compared to non-cancer persons ("-"without Tx "+" with Tx)')
% % %     OUTPUT_DATA{i,18}=RRInsitu(1,i);
% % %     OUTPUT_DATA{i,19}=RRLocal(1,i);
% % %     OUTPUT_DATA{i,20}=RRRegional(1,i);
% % %     OUTPUT_DATA{i,21}=RRDistant(1,i); 
% % %     
% % %     OUTPUT_DATA{i,23}=RI(1,i);
% % %     OUTPUT_DATA{i,24}=RL(1,i);
% % %     OUTPUT_DATA{i,25}=RR(1,i);
% % %     OUTPUT_DATA{i,26}=RD(1,i); 
% % %     
% % % end

%OUTPUT_DATA


% figure
% hold on
% plot(1:MaxAge, Coverage(1:MaxAge,1), 'b');
% plot(1:MaxAge, Coverage(1:MaxAge,2), 'g');
% plot(1:MaxAge,Coverage(1:MaxAge,3), 'c');
% plot(1:MaxAge,Coverage(1:MaxAge,4), 'r');
% legend('in-situ', 'local', 'regional', 'distant');
% title(' Treatment Coverage')
% 
% figure
% hold on
% plot(1:MaxAge, Coverage(1:MaxAge,5), 'b');
% plot(1:MaxAge, Coverage(1:MaxAge,6), 'g');
% plot(1:MaxAge,Coverage(1:MaxAge,7), 'c');
% plot(1:MaxAge,Coverage(1:MaxAge,8), 'r');
% legend('in-situ', 'local', 'regional', 'distant');
% title(' Stage Dist')

figure
hold on
plot(1:MaxAge, BaseMortality(1,1:MaxAge), 'b');
plot(1:MaxAge, BaseMortality(2,1:MaxAge), 'g');
plot(1:MaxAge,BaseMortality(3,1:MaxAge), 'c');
plot(1:MaxAge,BaseMortality(4,1:MaxAge), 'r');
% plot(1:MaxAge, CancerMort(1:MaxAge,1), 'c+');
legend('in-situ', 'local', 'regional', 'distant');
title(' Cancert mortality by stage')
figure
hold on
plot(1:MaxAge, Mu_Tx(1,1:MaxAge), 'b+');
plot(1:MaxAge, Mu_Tx(2,1:MaxAge), 'g+');
plot(1:MaxAge, Mu_Tx(3,1:MaxAge), 'y+');
plot(1:MaxAge, Mu_Tx(4,1:MaxAge), 'r+');

legend('in-situ', 'local', 'regional', 'distant');
plot(1:MaxAge, Mu_NoTx(1,1:MaxAge), 'b');
plot(1:MaxAge, Mu_NoTx(2,1:MaxAge), 'g');
plot(1:MaxAge, Mu_NoTx(3,1:MaxAge), 'y');
plot(1:MaxAge, Mu_NoTx(4,1:MaxAge), 'r');

title(' Mu-Treatment(+) and No Treatment(-)')

figure
hold on
plot(1:80, CancerMort(1:80,1) - MortDFbyAge(1,1:80)' , 'c');
plot(1:80, CancerMort(1:80,1), 'c+');
plot(1:80, calcCancerMort(1,1:80), 'r+');
legend('', 'cancerMortGLOBOCAN', 'calcCancerMort');

% title(' Mu-Cancer')
% 
% 
% 
% 
% figure
% hold on
%plot(1:MaxAge, MortalityArray(1,1:MaxAge), 'r');
title(' Mortalities')

CountryIndex = PAR{9};
CountryName = PAR{10};
% 
%  [A,B,C,D,E,F,G,H] = MikeWriteTemplateMort(TEMPMORT1,TEMPMORT2,TEMPMORT3,TEMPMORT4,TEMPMORT5,TEMPMORT6,TEMPMORT7,TEMPMORT8,Healthy,Insitu,Local,Regional,Distant,RMortalityNoTx,RMortalityTx,RRInsitu,RRLocal,RRRegional,RRDistant,RI,RL,RR,RD,Coverage,BaseMortality,Mu_Tx,Mu_NoTx,CancerMort,MortalityArray,CountryIndex,CountryName)
% 
% 
%   xlswrite('Mike Mort Results.xls',A,'10-year Survial');
%   xlswrite('Mike Mort Results.xls',B,'yrly reltv mort comp non cancer');
%   xlswrite('Mike Mort Results.xls',C,'10-yr reltv mort comp non cance');
%   xlswrite('Mike Mort Results.xls',D,'Treatment Coverage');
%    xlswrite('Mike Mort Results.xls',E,'Stage Dist');
%    xlswrite('Mike Mort Results.xls',F,'Cancer Mort by stage');
%    xlswrite('Mike Mort Results.xls',G,'Mu Treat and No Treat');
%    xlswrite('Mike Mort Results.xls',H,'Mortalities');
%   
% 
%   'Mike Mort Results Written'
  
%pause

end

function Proportion = TenYearSurvival(RatesArray, Mu)
RatesArray = RatesArray /12;
Mu = Mu / 12;
N = zeros (1,size(RatesArray,1));
N(1,1) = 10000;
for t = 1:10 * 12
    N(1,1) = N(1,1) - N(1,1) * (RatesArray(1,max(1,round(t/12)))+ Mu(1,max(1,round(t/12))));
    for i =  size(RatesArray):-1:2
    N(1,i)= N(1,i) + N(1,i-1) * RatesArray(i-1,max(1,round(t/12))) - N(1,i) * (RatesArray(i, max(1,round(t/12))) + Mu(1, max(1,round(t/12))));
    end
end


Proportion = sum(sum(N))/10000;
end





%% OLD CODE: 2015 version
%function [Mu_Tx, Impact, BaseMortality, TxCoverage, StgDist, RMortalityTx, RMortalityNoTx, ModDwellTime,A,B,C,D,E,F,G] = Mortality(DwellTime, MortalityArray, Prevalence, StageDist, CancerMort,PAR)
% global CancerMortality 
% global CancerMortalityCalc 
% %global StageDistribution
% MaxAge = 101;
% 
% % [~,~,TEMPMARKOV]=xlsread('C:\Users\Mike\Documents\Umass\Spring 2015\Cancer Research\excel\Template for OutputData Markov.xls');
% [~,~,TEMPMORT1]=xlsread('Template for OutputData Mort.xls','10-year Survial');
% [~,~,TEMPMORT2]=xlsread('Template for OutputData Mort.xls','yrly reltv mort comp non cancer');
% [~,~,TEMPMORT3]=xlsread('Template for OutputData Mort.xls','10-yr reltv mort comp non cance');
% [~,~,TEMPMORT4]=xlsread('Template for OutputData Mort.xls','Treatment Coverage');
%  [~,~,TEMPMORT5]=xlsread('Template for OutputData Mort.xls','Stage Dist');
%  [~,~,TEMPMORT6]=xlsread('Template for OutputData Mort.xls','Cancer Mort by stage');
%  [~,~,TEMPMORT7]=xlsread('Template for OutputData Mort.xls','Mu Treat and No Treat');
%  [~,~,TEMPMORT8]=xlsread('Template for OutputData Mort.xls','Mortalities');
% 
% %StageDistribution = StageDist;
% Healthy = zeros(1,MaxAge);
% Insitu = zeros(1,MaxAge);
% Local = zeros(1,MaxAge);
% Regional = zeros(1,MaxAge);
% Distant = zeros(1,MaxAge);
% RRInsitu = zeros(1,MaxAge);
% RRLocal = zeros(1,MaxAge);
% RRRegional = zeros(1,MaxAge);
% RRDistant = zeros(1,MaxAge);
% 
% global Mu_TxVal 
% Mu_TxVal = zeros(4,1);
% Mu_Tx = zeros(4,MaxAge);
% global Mu_NoTxVal 
% Mu_NoTxVal = zeros(4,1);
% Mu_NoTx = zeros(4,MaxAge);
% Coverage = zeros(MaxAge,8);
% 
% RInsitu_Tx = .99;
% RLocal_Tx = 0.94;
% RRegional_Tx = 0.60;
% RDistant_Tx = 0.12;
% 
% %Okonkwo et. al, JNCI 2008		
% RNoTx = [0.079, 0.136, 0.231, 0.500];
% %Zelle et. al (original Groot et. al., 2006)	
% RTx = [0.006, .020666667, 0.084, .27];
% 
% %Below Commented by Prashant on 10 Jul 2016. 
% % %Groot, The Breast Journal, Volume 12 Suppl. 1, 2006; For African region
% % value = mean( [0.042, 0.093]);
% % RTx = [0.006, 0.006, value, 0.275];
% % value = mean([0.063,.15]);
% % RNoTx = [0.02, 0.02, value, 0.3];
% 
% RTxByNoTx = RTx./RNoTx;
% 
% 
% ModDwellTime = DwellTime;
% 
% RMortalityTx = zeros(4,MaxAge);
% RMortalityNoTx = zeros(4,MaxAge);
%  t = 10;
% for i = 1:MaxAge - 10
%     
%    DwellRate(1,i:i+9)= 1 ./ DwellTime(1,i:i+9);
%    DwellRate(2,i:i+9)= 1 ./ DwellTime(2,i:i+9);
%    DwellRate(3,i:i+9)= 1 ./ DwellTime(3,i:i+9);
%    DwellRate(4,i:i+9)= 1 ./ DwellTime(4,i:i+9);
%   
% % Ten-year survival without treatment
%    Healthy(1,i) = TenYearSurvival(zeros(1,10),MortalityArray(1,i:i+9));
%    Insitu(1,i) = TenYearSurvival(DwellRate(1:4,i:i+9),MortalityArray(1,i:i+9));
%    Local(1,i) = TenYearSurvival(DwellRate(2:4,i:i+9),MortalityArray(1,i:i+9));
%    Regional(1,i) = TenYearSurvival(DwellRate(3:4,i:i+9),MortalityArray(1,i:i+9));
%    Distant(1,i) = TenYearSurvival(DwellRate(4:4,i:i+9),MortalityArray(1,i:i+9));
%    
%    
% % Ten-year RELATIVE survival without treatment   
%    RRInsitu(1,i) = Insitu(1,i) / Healthy(1,i);
%    RRLocal(1,i) = Local(1,i) / Healthy(1,i);
%    RRRegional(1,i) = Regional(1,i) / Healthy(1,i);
%    RRDistant(1,i) = Distant(1,i) / Healthy(1,i);
%   
%   %Yearly mortality rate without treatment = mortality from disease + natural mortality(Pr(mortality > RRInsitu) = exp^-rate 
%    Mu_NoTx(1,i) = (MortalityArray(1,i)*t - log(RRInsitu(1,i)))/t;
%    Mu_NoTx(2,i) = (MortalityArray(1,i)*t - log(RRLocal(1,i)))/t;
%    Mu_NoTx(3,i) = (MortalityArray(1,i)*t - log(RRRegional(1,i)))/t;
%    Mu_NoTx(4,i) = (MortalityArray(1,i)*t - log(RRDistant(1,i)))/t;
%      
%   
% end
% 
% Mu_NoTx(1,91:101) = interp1(1:90,Mu_NoTx(1,1:90),(91:101),'cubic','extrap');
% Mu_NoTx(2,91:101) = interp1(1:90,Mu_NoTx(2,1:90),(91:101),'cubic','extrap');
% Mu_NoTx(3,91:101) = interp1(1:90,Mu_NoTx(3,1:90),(91:101),'cubic','extrap');
% Mu_NoTx(4,91:101) = interp1(1:90,Mu_NoTx(4,1:90),(91:101),'cubic','extrap');
% 
% 
%   
% %%%%
% % for i = MaxAge - 10 : MaxAge
% %    RRInsitu(1,i) =  RRInsitu(1,i - 1); 
% %    RRLocal(1,i) = RRLocal(1,i - 1);
% %    RRRegional(1,i) = RRRegional(1,i - 1);
% %    RRDistant(1,i) = RRDistant(1,i - 1);
% %    
% %    Mu_NoTx(1,i) = Mu_NoTx(1,i-1) ;
% %    Mu_NoTx(2,i) = Mu_NoTx(2,i-1) ;
% %    Mu_NoTx(3,i) = Mu_NoTx(3,i-1) ;
% %    Mu_NoTx(4,i) = Mu_NoTx(4,i-1) ;
%   
% %end
% 
% 'Ten year'
% % InsituCF = sum(dot (RRInsitu (1:MaxAge), Prevalence(1:MaxAge)) )/ sum(Prevalence(1:MaxAge))
% % LocalCF = sum(dot (RRLocal (1:MaxAge), Prevalence(1:MaxAge)) )/ sum(Prevalence(1:MaxAge))
% % RegionalCF =  sum(dot (RRRegional (1:MaxAge), Prevalence(1:MaxAge)) )/ sum(Prevalence(1:MaxAge))
% % DistantCF =  sum(dot (RRDistant(1:MaxAge), Prevalence(1:MaxAge)) )/ sum(Prevalence(1:MaxAge))
% 
% 
% InsituCF = mean(RRInsitu);%(1,30:80));
% LocalCF = mean(RRLocal);%(1,30:80));
% RegionalCF = mean(RRRegional);%(1,30:80));
% DistantCF = mean(RRDistant);%(1,30:80));
% 
% TotalDFMu = sum(dot(MortalityArray(1,1:MaxAge),Prevalence(1:MaxAge)))/ sum(Prevalence(1:MaxAge))
% 
% RI_F = RInsitu_Tx / InsituCF;
% RL_F = RLocal_Tx / LocalCF;
% RR_F = RRegional_Tx / RegionalCF;
% RD_F = RDistant_Tx / DistantCF; 
% 
% 
% for i = 1:MaxAge
% %     % 10-year Relative surival with_Tx
% %     RI(1,i) = min(1,RI_F * RRInsitu(1,i));
% %     RL(1,i) = min(1,RL_F * RRLocal(1,i));
% %     RR(1,i) = min(1,RR_F * RRRegional(1,i));
% %     RD(1,i) = min(1,RD_F * RRDistant(1,i));
% %     
% %   %  Yearly Mortality rate with treatment
% %     Mu_Tx(1,i) = (MortalityArray(1,i)*t - log(RI(1,i)))/t;
% %     Mu_Tx(2,i) = (MortalityArray(1,i)*t - log(RL(1,i)))/t;
% %     Mu_Tx(3,i) = (MortalityArray(1,i)*t - log(RR(1,i)))/t;
% %     Mu_Tx(4,i) = (MortalityArray(1,i)*t - log(RD(1,i)))/t;
%    
% 
%         Mu_Tx(:,i) = max(0,(Mu_NoTx(:,i) - MortalityArray(1,i))) .* RTxByNoTx(:) + MortalityArray(1,i);
%         
%   
%     
%      %Adjusting mortality with treatment to match mortality rate from
% %     %country data
% %     R1 =(1- RLocal_Tx) /  (1 - RRegional_Tx);% Mu_NoTx(2,i)/Mu_NoTx(3,i);
% %     R2 = (1- RLocal_Tx) /  (1 - RDistant_Tx);% Mu_NoTx(2,i)/Mu_NoTx(4,i);
% %     R3 =  (1- RLocal_Tx) /  (1 - RInsitu_Tx);%Mu_NoTx(2,i)/Mu_NoTx(1,i);
%     
%         R1 = RTx(2) / RTx(3) ;
%         R2 = RTx(2) / RTx(4) ;
%         R3 = RTx(2) / RTx(1) ;
%    %If total mortality greater estimate mortalit adjust (using stage 2 as refernce, all otehr stages represented relative to stage 2)     
% %     if ((CancerMort(i,1) + MortalityArray(1,i)) / (StageDist(1,2)+ (StageDist(1,3)/R1)+ (StageDist(1,4)/R2 ) + (StageDist(1,1)/R3 ) )) > 1.2 * Mu_Tx(2,i)
% %         Mu_Tx(2,i) = ((CancerMort(i,1) + MortalityArray(1,i)) / (StageDist(1,2)+ (StageDist(1,3)/R1)+ (StageDist(1,4)/R2 ) + (StageDist(1,1)/R3 ) ));
% %         Mu_Tx(3,i) = Mu_Tx(2,i)/ R1;
% %         Mu_Tx(4,i) = Mu_Tx(2,i)/ R2;
% %         Mu_Tx(1,i) = Mu_Tx(2,i)/ R3;
% %          
% %         ValueF(1:4,1) = (Mu_Tx(:,i) ./ RTxByNoTx(:))./ Mu_NoTx(:,i);
% %         
% %         Mu_NoTx(:,i) = Mu_Tx(:,i) ./ RTxByNoTx(:);
% %         
% %         ModDwellTime(1:4,i) =  ModDwellTime(1:4,i)./ValueF(1:4,1);
% %        
% %    
% %     end
%         
% 
%      % 10-year Relative surival with_Tx
%    
%         RI(1,i) =  min(1,exp(MortalityArray(1,i)*t- Mu_Tx(1,i)*t));
%         RL(1,i) =  min(1,exp(MortalityArray(1,i)*t- Mu_Tx(2,i)*t));
%         RR(1,i) =  min(1,exp(MortalityArray(1,i)*t- Mu_Tx(3,i)*t));
%         RD(1,i) =  min(1,exp(MortalityArray(1,i)*t- Mu_Tx(4,i)*t));
%       
%     
%     CancerMortality = CancerMort(i,1) + MortalityArray(1,i);
%     Mu_TxVal(:,1) = Mu_Tx(:,i);
%     Mu_NoTxVal(:,1) = Mu_NoTx(:,i);
%     
%     if StageDist(1,1) == 0
%         uVal = 0;
%     else
%         uVal = 1;
%     end
%     
%     
%     % Assuming Tx coverage is 100% 
%     LowerBound = [uVal 1 1 1  0 0 0 0  ];
%     UpperBound = [uVal 1 1 1  uVal 1 1 1 ];
%     options = optimoptions(@fmincon,'MaxIter', 1000, 'TolX', 0.001);
%     %Covergae: 1:4= coverage of Tx in  stage insitu, local, regional, and distant; 4:6 =
%     %Stage distribution
%     [Coverage(i, 1:8), error(i)] = fmincon(@FuncCancerMort, [0 1 1 1  0 0.3 0.4 .3 ], [], [], [0 0 0 0  1 1 1 1 ], [1],LowerBound,UpperBound,[] );
%     %Coverage(6,i) = 1- Coverage(5,i) - Coverage (4,i);
%     Cancer(1,i) = CancerMortalityCalc;
% 
%     for j = 1: 4
%     if Coverage(i, j) < 0.01
%        Coverage(i, j) = 0; 
%     end
%     end
% end
%     %Subtracting disease free mortality from cancer mortality as it is added
%     %back in Spectrum
% %   
%     
%     Mu_Tx(1,:) = max(0,Mu_Tx(1,:) - MortalityArray(1,:));
%     Mu_Tx(2,:) = max(0,Mu_Tx(2,:) - MortalityArray(1,:));
%     Mu_Tx(3,:) = max(0,Mu_Tx(3,:) - MortalityArray(1,:));
%     Mu_Tx(4,:) = max(0,Mu_Tx(4,:) - MortalityArray(1,:));
%     
%     Mu_NoTx(1,:) = max(0,Mu_NoTx(1,:) - MortalityArray(1,:));
%     Mu_NoTx(2,:) = max(0,Mu_NoTx(2,:) - MortalityArray(1,:));
%     Mu_NoTx(3,:) = max(0,Mu_NoTx(3,:) - MortalityArray(1,:));
%     Mu_NoTx(4,:) = max(0,Mu_NoTx(4,:) - MortalityArray(1,:));
%     
% %     Mu_Tx(1,80:end) = Mu_Tx(1,80);
% %     Mu_Tx(2,80:end) = Mu_Tx(2,80);
% %     Mu_Tx(3,80:end) = Mu_Tx(3,80);
% %     Mu_Tx(4,80:end) = Mu_Tx(4,80);
% %     
% %     Mu_NoTx(1,80:end) = Mu_NoTx(1,80);
% %     Mu_NoTx(2,80:end) = Mu_NoTx(2,80);
% %     Mu_NoTx(3,80:end) = Mu_NoTx(3,80);
% %     Mu_NoTx(4,80:end) = Mu_NoTx(4,80);
%    
%     RMortalityTx(1,:)= Mu_Tx(1,:)./ MortalityArray(1,:);
%     RMortalityTx(2,:)= Mu_Tx(2,:)./ MortalityArray(1,:);
%     RMortalityTx(3,:)= Mu_Tx(3,:)./ MortalityArray(1,:);
%     RMortalityTx(4,:)= Mu_Tx(4,:)./ MortalityArray(1,:);
%     
%     RMortalityNoTx(1,:)= Mu_NoTx(1,:)./ MortalityArray(1,:);
%     RMortalityNoTx(2,:)= Mu_NoTx(2,:)./ MortalityArray(1,:);
%     RMortalityNoTx(3,:)= Mu_NoTx(3,:)./ MortalityArray(1,:);
%     RMortalityNoTx(4,:)= Mu_NoTx(4,:)./ MortalityArray(1,:); 
%    
% 
% % for i = 1: MaxAge
% %     for j = 1: 4
% %         if ((MortalityArray(1,i)*t - log(RRInsitu(1,i)))/t) < Mu_NoTx(1,i) 
% %             
% %                LowerBound = [0 ];
% %                UpperBound = [40];
% %                 options = optimoptions(@fmincon,'MaxIter', 1000, 'TolX', 0.001);
% %                [ModDwellTime(j,i)] = fmincon(@FuncCancerMort, [], [], [], [0 0 0 0 1 1 1 1], [1],LowerBound,UpperBound,[] );
% % 
% % end
% 
% % Coverage(80:100,1) = Coverage(80,1);
% % Coverage(80:100,2) = Coverage(80,2);
% % Coverage(80:100,3) = Coverage(80,3);
% % Coverage(80:100,4) = Coverage(80,4);
% % Calculating basecase mortality, i.e., the mortality weighted by
% % treatment coverage
% % 
% 
% BaseMortality = zeros(size(Mu_Tx,1), size(Mu_Tx,2));
% BaseMortality(1,:) = (Coverage(:,1)' .*  Mu_Tx(1,:) + (1- Coverage(:,1)').* Mu_NoTx(1,:))  ;
% BaseMortality(2,:) = (Coverage(:,2)' .*  Mu_Tx(2,:) + (1- Coverage(:,2)').* Mu_NoTx(2,:))  ;
% BaseMortality(3,:) = (Coverage(:,3)' .*  Mu_Tx(3,:) + (1- Coverage(:,3)').* Mu_NoTx(3,:))  ;
% BaseMortality(4,:) = (Coverage(:,4)' .*  Mu_Tx(4,:) + (1- Coverage(:,4)').* Mu_NoTx(4,:)) ;
% 
% % estimating increase in weighted mortality for unit change in coverage
% Impact (1,:) = (Mu_Tx(1,:) - Mu_NoTx(1,:)) ./ BaseMortality(1,:);
% Impact (2,:) = (Mu_Tx(2,:) - Mu_NoTx(2,:)) ./ BaseMortality(2,:);
% Impact (3,:) = (Mu_Tx(3,:) - Mu_NoTx(3,:)) ./ BaseMortality(3,:);
% Impact (4,:) = (Mu_Tx(4,:) - Mu_NoTx(4,:)) ./ BaseMortality(4,:);
% 
% TxCoverage(:,:) = max(0,Coverage(:,1:4));
% StgDist(:,:) =  max(0,Coverage(:,5:8));
% close all
% figure
% hold on
% plot(1:MaxAge, Healthy(1,1:MaxAge), 'r');
% plot(1:MaxAge, Insitu(1,1:MaxAge), 'y');
% plot(1:MaxAge, Local(1,1:MaxAge), 'b');
% plot(1:MaxAge, Regional(1,1:MaxAge), 'g');
% plot(1:MaxAge, Distant(1,1:MaxAge), 'c');
% legend('healthy', 'in-situ', 'local', 'regional', 'distant');
% % plot(30:80, T_Local(1,30:80), 'b-');
% % plot(30:80, T_Regional(1,30:80), 'g-');
% % plot(30:80, T_Distant(1,30:80), 'c-');
% title('10-year Survival')
% 
% % % % for row=1:length(Healthy)
% % % %     for col=1:MaxAge-10
% % % %         if row == 1
% % % %             OUTPUT_DATA{row,col} = {row,col};
% % % %         else
% % % %             OUTPUT_DATA{row,col} = [row,Healthy(col)];
% % % %         end
% % % %         
% % % %     end
% % % % end
% 
% 
% % % % %title('10-year Survival')
% % % % for i = 1 : MaxAge - 10
% % % %     OUTPUT_DATA{i,1}=i;
% % % %     OUTPUT_DATA{i,2}=Healthy(1,i);
% % % %     OUTPUT_DATA{i,3}=Insitu(1,i);
% % % %     OUTPUT_DATA{i,4}=Local(1,i);
% % % %     OUTPUT_DATA{i,5}=Regional(1,i);
% % % %     OUTPUT_DATA{i,6}=Distant(1,i);    
% % % % end
% % % % 
% % % % for i = MaxAge - 10 : MaxAge
% % % %     OUTPUT_DATA{i,1}=i;
% % % % end
% 
% % OUTPUT_DATA
% 
% figure
% hold on
% plot(1:MaxAge, RMortalityNoTx(1,1:MaxAge), 'y');
% plot(1:MaxAge, RMortalityNoTx(2,1:MaxAge), 'b');
% plot(1:MaxAge, RMortalityNoTx(3,1:MaxAge), 'g');
% plot(1:MaxAge,RMortalityNoTx(4,1:MaxAge), 'c');
% 
% legend('in-situ', 'local', 'regional', 'distant');
% plot(1:MaxAge, RMortalityTx(1,1:MaxAge), 'y+');
% plot(1:MaxAge, RMortalityTx(2,1:MaxAge), 'b+');
% plot(1:MaxAge, RMortalityTx(3,1:MaxAge), 'g+');
% plot(1:MaxAge, RMortalityTx(4,1:MaxAge), 'c+');
% 
% title('yearly relative mortality in cancer comaed to non-cancer persons ("-"without Tx "+" with Tx)')
% 
% % % % for i = 1 : MaxAge
% % % %     
% % % %     %legend('in-situ', 'local', 'regional', 'distant');
% % % %     OUTPUT_DATA{i,8}=RMortalityNoTx(1,i);
% % % %     OUTPUT_DATA{i,9}=RMortalityNoTx(2,i);
% % % %     OUTPUT_DATA{i,10}=RMortalityNoTx(3,i);
% % % %     OUTPUT_DATA{i,11}=RMortalityNoTx(4,i); 
% % % %     
% % % %     OUTPUT_DATA{i,13}=RMortalityTx(1,i);
% % % %     OUTPUT_DATA{i,14}=RMortalityTx(2,i);
% % % %     OUTPUT_DATA{i,15}=RMortalityTx(3,i);
% % % %     OUTPUT_DATA{i,16}=RMortalityTx(4,i); 
% % % %     
% % % % end
% 
% %OUTPUT_DATA
% 
% figure
% hold on
% plot(1:MaxAge, RRInsitu(1,1:MaxAge), 'y');
% plot(1:MaxAge, RRLocal(1,1:MaxAge), 'b');
% plot(1:MaxAge, RRRegional(1,1:MaxAge), 'g');
% plot(1:MaxAge,RRDistant(1,1:MaxAge), 'c');
% legend('in-situ', 'local', 'regional', 'distant');
% plot(1:MaxAge, RI(1,1:MaxAge), 'y+');
% plot(1:MaxAge, RL(1,1:MaxAge), 'b+');
% plot(1:MaxAge, RR(1,1:MaxAge), 'g+');
% plot(1:MaxAge, RD(1,1:MaxAge), 'c+');
% title('10-year relative survival in cancer compared to non-cancer persons ("-"without Tx "+" with Tx)')
% 
% 
% 
% 
% % % % for i = 1 : MaxAge
% % % %     
% % % % %   title('10-year relative survival in cancer compared to non-cancer persons ("-"without Tx "+" with Tx)')
% % % %     OUTPUT_DATA{i,18}=RRInsitu(1,i);
% % % %     OUTPUT_DATA{i,19}=RRLocal(1,i);
% % % %     OUTPUT_DATA{i,20}=RRRegional(1,i);
% % % %     OUTPUT_DATA{i,21}=RRDistant(1,i); 
% % % %     
% % % %     OUTPUT_DATA{i,23}=RI(1,i);
% % % %     OUTPUT_DATA{i,24}=RL(1,i);
% % % %     OUTPUT_DATA{i,25}=RR(1,i);
% % % %     OUTPUT_DATA{i,26}=RD(1,i); 
% % % %     
% % % % end
% 
% %OUTPUT_DATA
% 
% 
% figure
% hold on
% plot(1:MaxAge, Coverage(1:MaxAge,1), 'b');
% plot(1:MaxAge, Coverage(1:MaxAge,2), 'g');
% plot(1:MaxAge,Coverage(1:MaxAge,3), 'c');
% plot(1:MaxAge,Coverage(1:MaxAge,4), 'r');
% legend('in-situ', 'local', 'regional', 'distant');
% title(' Treatment Coverage')
% 
% figure
% hold on
% plot(1:MaxAge, Coverage(1:MaxAge,5), 'b');
% plot(1:MaxAge, Coverage(1:MaxAge,6), 'g');
% plot(1:MaxAge,Coverage(1:MaxAge,7), 'c');
% plot(1:MaxAge,Coverage(1:MaxAge,8), 'r');
% legend('in-situ', 'local', 'regional', 'distant');
% title(' Stage Dist')
% 
% figure
% hold on
% plot(1:MaxAge, BaseMortality(1,1:MaxAge), 'b');
% plot(1:MaxAge, BaseMortality(2,1:MaxAge), 'g');
% plot(1:MaxAge,BaseMortality(3,1:MaxAge), 'c');
% plot(1:MaxAge,BaseMortality(4,1:MaxAge), 'r');
% % plot(1:MaxAge, CancerMort(1:MaxAge,1), 'c+');
% legend('in-situ', 'local', 'regional', 'distant');
% title(' Cancert mortality by stage')
% figure
% hold on
% plot(1:MaxAge, Mu_Tx(1,1:MaxAge), 'b+');
% plot(1:MaxAge, Mu_Tx(2,1:MaxAge), 'g+');
% plot(1:MaxAge, Mu_Tx(3,1:MaxAge), 'y+');
% plot(1:MaxAge, Mu_Tx(4,1:MaxAge), 'r+');
% 
% legend('in-situ', 'local', 'regional', 'distant');
% plot(1:MaxAge, Mu_NoTx(1,1:MaxAge), 'b');
% plot(1:MaxAge, Mu_NoTx(2,1:MaxAge), 'g');
% plot(1:MaxAge, Mu_NoTx(3,1:MaxAge), 'y');
% plot(1:MaxAge, Mu_NoTx(4,1:MaxAge), 'r');
% 
% title(' Mu-Treatment(+) and No Treatment(-)')
% 
% figure
% hold on
% plot(1:80, CancerMort(1:80,1) - MortalityArray(1,1:80)' , 'c');
% plot(1:80, CancerMort(1:80,1), 'c+');
% % title(' Mu-Cancer')
% % 
% % 
% % 
% % 
% % figure
% % hold on
% %plot(1:MaxAge, MortalityArray(1,1:MaxAge), 'r');
% title(' Mortalities')
% 
% CountryIndex = PAR{9};
% CountryName = PAR{10};
% % 
% %  [A,B,C,D,E,F,G,H] = MikeWriteTemplateMort(TEMPMORT1,TEMPMORT2,TEMPMORT3,TEMPMORT4,TEMPMORT5,TEMPMORT6,TEMPMORT7,TEMPMORT8,Healthy,Insitu,Local,Regional,Distant,RMortalityNoTx,RMortalityTx,RRInsitu,RRLocal,RRRegional,RRDistant,RI,RL,RR,RD,Coverage,BaseMortality,Mu_Tx,Mu_NoTx,CancerMort,MortalityArray,CountryIndex,CountryName)
% % 
% % 
% %   xlswrite('Mike Mort Results.xls',A,'10-year Survial');
% %   xlswrite('Mike Mort Results.xls',B,'yrly reltv mort comp non cancer');
% %   xlswrite('Mike Mort Results.xls',C,'10-yr reltv mort comp non cance');
% %   xlswrite('Mike Mort Results.xls',D,'Treatment Coverage');
% %    xlswrite('Mike Mort Results.xls',E,'Stage Dist');
% %    xlswrite('Mike Mort Results.xls',F,'Cancer Mort by stage');
% %    xlswrite('Mike Mort Results.xls',G,'Mu Treat and No Treat');
% %    xlswrite('Mike Mort Results.xls',H,'Mortalities');
% %   
% % 
% %   'Mike Mort Results Written'
%   
% %pause
% end

%function Error = FuncCancerMort(coverage)
% global Mu_TxVal 
% global Mu_NoTxVal
% %global StageDistribution 
% global CancerMortality
% global CancerMortalityCalc
% CancerMortalityCalc =  (coverage(1,1) *  Mu_TxVal(1,1) + (1- coverage(1,1))* Mu_NoTxVal(1,1)) * coverage(1,5) + ...
%                        (coverage(1,2) *  Mu_TxVal(2,1) + (1- coverage(1,2))* Mu_NoTxVal(2,1)) * coverage(1,6) + ...
%                        (coverage(1,3) *  Mu_TxVal(3,1) + (1- coverage(1,3))* Mu_NoTxVal(3,1)) * coverage(1,7) + ...
%                        (coverage(1,4) *  Mu_TxVal(4,1) + (1- coverage(1,4))* Mu_NoTxVal(4,1)) * coverage(1,8);
%        
% Error = CancerMortality - CancerMortalityCalc;
% Error = Error * Error;
% end


