function [A,B,C,D,E,F,G,H] = MikeWriteTemplateMort(TEMPMORT1,TEMPMORT2,TEMPMORT3,TEMPMORT4,TEMPMORT5,TEMPMORT6,TEMPMORT7,TEMPMORT8,Healthy,Insitu,Local,Regional,Distant,RMortalityNoTx,RMortalityTx,RRInsitu,RRLocal,RRRegional,RRDistant,RI,RL,RR,RD,Coverage,BaseMortality,Mu_Tx,Mu_NoTx,CancerMort,MortalityArray,CountryIndex,CountryName)

%  [~,~,TEMPMORT]=xlsread('C:\Users\Mike\Documents\Umass\Spring 2015\Cancer Research\excel\Template for OutputData Markov.xls');
%  [~,~,TEMPMARKOV]=xlsread('C:\Users\Mike\Documents\Umass\Spring 2015\Cancer Research\excel\Template for OutputData Mort.xls');

            A = TEMPMORT1;
            B = TEMPMORT2;
            C = TEMPMORT3;
            D = TEMPMORT4;
            E = TEMPMORT5;
            F = TEMPMORT6;
            G = TEMPMORT7;
            H = TEMPMORT8;
            
            
%             CountryIndex = 1;
%             CountryIndex=cell2mat(CountryIndex);
            CountryIndex = CountryIndex - 1

            
            
%             
% PrevS=RES{1};
% PrevC=RES{2};
% Incidence=RES{3}
% OnsetRates=RES{4}
% OnsetRates(end)=OnsetRates(end-1);
% DiagnosticRate=1./RES{5}
% DwellRate=1./RES{6}
% MortRate= RES{7};%Relative Mortality without treatment
% MortRateDF= RES{8};% Disease free mortality
% Impact = RES{9};%Relative Mortality with treatment
% TxCoverage = RES{10};% Current coverage of treatment
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Start filling in mortality %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if  CountryIndex == 0
%     P = 0;
% else
%     P = 1;
% end
    
%     sA = (CountryIndex*5)+ P
sA = CountryIndex .* 6;
A{4,5+sA} = CountryName;

for icol=5:9
    if icol ==5
        for irow=7:107
            A{irow,(icol+sA)}=Healthy(1,irow-6);
        end
    elseif icol==6 %+(CountryIndex*5)+1
        for irow=7:107
            A{irow,(icol+sA)}=Insitu(1,irow-6);
        end
    elseif icol==7 %+(CountryIndex*5)+1
        for irow=7:107
            A{irow,(icol+sA)}=Local(1,irow-6);
        end
    elseif icol==8 %+(CountryIndex*5)+1
        for irow=7:107
            A{irow,(icol+sA)}=Regional(1,irow-6);
        end
    elseif icol==9
        for irow=7:107
            A{irow,(icol+sA)}=Distant(1,irow-6);
        end
    end
end

size(Healthy);
Healthy;

A;


%%%%%%%%%%%%%%%%%%%%%
% sBG = (CountryIndex*9)+P;
sBG = (CountryIndex .* 10);
B{4,5+sBG} = CountryName;
G{4,5+sBG} = CountryName;

for icol=5:8
        for irow=8:108
            B{irow,icol+sBG}=RMortalityNoTx(icol-4,irow-7);
            G{irow,icol+sBG}=Mu_Tx(icol-4,irow-7);
        end
end

for icol=10:13
        for irow=8:108
            B{irow,icol+sBG}=RMortalityTx(icol-9,irow-7);
            G{irow,icol+sBG}=Mu_NoTx(icol-9,irow-7);
        end
end

B;
G;
%%%%%%%%%%%%%%%%%%%%%
% sC = (CountryIndex*9)+P;
sC = (CountryIndex .* 10);
C{4,5+sC} = CountryName;

% for icol=5:8
%     if icol ==5
%         for irow=8:108
%             C{irow,icol+sC}=RRInsitu(1,irow-7);
%         end
%     elseif icol==6
%         for irow=8:108
%             C{irow,icol+sC}=RRLocal(1,irow-7);
%         end
%     elseif icol==7
%         for irow=8:108
%             C{irow,icol+sC}=RRRegional(1,irow-7);
%         end
%     else
%         for irow=8:108
%             C{irow,icol+sC}=RRDistant(1,irow-7);
%         end
%     end
% end


%trial stuff
for icol=5:8
    if icol ==5
        for irow=8:108
            C{irow,icol+sC}=RRInsitu(1,irow-7);
        end
    elseif icol==6
        for irow=8:108
            C{irow,icol+sC}=RRLocal(1,irow-7);
        end
    elseif icol==7
        for irow=8:108
            C{irow,icol+sC}=RRRegional(1,irow-7);
        end
    elseif icol==8
        for irow=8:108
            C{irow,icol+sC}=RRDistant(1,irow-7);
        end
    end
end

%%%%%%

for icol=10:13
    if icol ==10
        for irow=8:108
            C{irow,icol+sC}=RI(1,irow-7);
        end
    elseif icol==11
        for irow=8:108
            C{irow,icol+sC}=RL(1,irow-7);
        end
    elseif icol==12
        for irow=8:108
            C{irow,icol+sC}=RR(1,irow-7);
        end
    elseif icol==13
        for irow=8:108
            C{irow,icol+sC}=RD(1,irow-7);
        end
    end
end

C;
%%%%%%%%%%%%%%%%%%%%%
% sDEF = (CountryIndex*4)+P;
sDEF = (CountryIndex .* 5);
D{4,5+sDEF} = CountryName;
E{4,5+sDEF} = CountryName;
F{4,5+sDEF} = CountryName;

for icol=5:8
        for irow=8:108
            D{irow,icol+sDEF}=Coverage(irow-7,icol-4);
            E{irow,icol+sDEF}=Coverage(irow-7,icol);
            F{irow,icol+sDEF}=BaseMortality(icol-4,irow-7);
        end
end

D;
E;
F;
%%%%%%%%%%%%%%%%%%%%
% sH = (CountryIndex*2)+P;
sH = (CountryIndex .* 3);
H{4,5+sH} = CountryName;

for irow=7:86
    H{irow,5+sH}=CancerMort(irow-6,1) - MortalityArray(1,irow-6);
    H{irow,6+sH}=CancerMort(irow-6,1);
end
H;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% End filling in mortality %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Start filling in Markov %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% End filling in Markov %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end
    
    
% % %     function [OUTPUT_DATA]=WriteOutputData(RES)
% % % 
% % % OUTPUT_DATA = [];
% % % 
% % % %Data from Markov and Cohort model
% % % 
% % % PrevS=RES{1};
% % % PrevC=RES{2};
% % % Incidence=RES{3}
% % % OnsetRates=RES{4}
% % % OnsetRates(end)=OnsetRates(end-1);
% % % DiagnosticRate=1./RES{5}
% % % DwellRate=1./RES{6}
% % % MortRate= RES{7};%Relative Mortality without treatment
% % % MortRateDF= RES{8};% Disease free mortality
% % % Impact = RES{9};%Relative Mortality with treatment
% % % TxCoverage = RES{10};% Current coverage of treatment
% % % 
% % % PrevS(1:4,:)= 0;
% % % PrevC(1:4,:) =0;
% % % OnsetRates(1:4)= 0;
% % % %Impact(:,1:6)=0;
% % % 
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % 
% % % % % % %save results to be mapped to associaitons file
% % % % % % RES{1}=PrevS;
% % % % % % RES{2}=PrevC;
% % % % % % RES{3}=Incidence;
% % % % % % RES{4}=OnsetRates';%onsetRateSmooth';
% % % % % % RES{5}=DiagnosticRate;
% % % % % % RES{6}=DwellTimes';
% % % % % % RES{7}=BaseMortalityC; %Relative Mortality without treatmet
% % % % % % RES{8}=MortA';% Disease free mortality
% % % % % % RES{9} = Impact;% %Relative Mortality with treatment
% % % % % % RES{10} = TxCoverage_A;% Current coverage of treatment
% % % 
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % 
% % % %Figures from mortality.m (8 Figures
% % % % % % figure
% % % % % % hold on
% % % % % % plot(1:MaxAge, Healthy(1,1:MaxAge), 'r');
% % % % % % plot(1:MaxAge, Insitu(1,1:MaxAge), 'y');
% % % % % % plot(1:MaxAge, Local(1,1:MaxAge), 'b');
% % % % % % plot(1:MaxAge, Regional(1,1:MaxAge), 'g');
% % % % % % plot(1:MaxAge, Distant(1,1:MaxAge), 'c');
% % % % % % legend('healthy', 'in-situ', 'local', 'regional', 'distant');
% % % % % % % plot(30:80, T_Local(1,30:80), 'b-');
% % % % % % % plot(30:80, T_Regional(1,30:80), 'g-');
% % % % % % % plot(30:80, T_Distant(1,30:80), 'c-');
% % % % % % title('10-year Survival')
% % % 
% % % MaxAge = 91
% % % 
% % % for row=1:length(Healthy)
% % %     for col=1:MaxAge
% % %         if row == 1
% % %             OUTPUT_DATA{row,col} = {row,col};
% % %         else
% % %             OUTPUT_DATA{row,col} = [row,Healthy(col)];
% % %         end
% % %         
% % %     end
% % % end
% % % 
% % % OUTPUT_DATA
% % % 
% % % pause
% % % 
% % % figure
% % % hold on
% % % plot(1:MaxAge, RMortalityNoTx(1,1:MaxAge), 'y');
% % % plot(1:MaxAge, RMortalityNoTx(2,1:MaxAge), 'b');
% % % plot(1:MaxAge, RMortalityNoTx(3,1:MaxAge), 'g');
% % % plot(1:MaxAge,RMortalityNoTx(4,1:MaxAge), 'c');
% % % 
% % % legend('in-situ', 'local', 'regional', 'distant');
% % % plot(1:MaxAge, RMortalityTx(1,1:MaxAge), 'y+');
% % % plot(1:MaxAge, RMortalityTx(2,1:MaxAge), 'b+');
% % % plot(1:MaxAge, RMortalityTx(3,1:MaxAge), 'g+');
% % % plot(1:MaxAge, RMortalityTx(4,1:MaxAge), 'c+');
% % % 
% % % title('yearly relative mortality in cancer compared to non-cancer persons ("-"without Tx "+" with Tx)')
% % % 
% % % figure
% % % hold on
% % % plot(1:MaxAge, RRInsitu(1,1:MaxAge), 'y');
% % % plot(1:MaxAge, RRLocal(1,1:MaxAge), 'b');
% % % plot(1:MaxAge, RRRegional(1,1:MaxAge), 'g');
% % % plot(1:MaxAge,RRDistant(1,1:MaxAge), 'c');
% % % legend('in-situ', 'local', 'regional', 'distant');
% % % plot(1:MaxAge, RI(1,1:MaxAge), 'y+');
% % % plot(1:MaxAge, RL(1,1:MaxAge), 'b+');
% % % plot(1:MaxAge, RR(1,1:MaxAge), 'g+');
% % % plot(1:MaxAge, RD(1,1:MaxAge), 'c+');
% % % title('10-year relative survival in cancer compared to non-cancer persons ("-"without Tx "+" with Tx)')
% % % 
% % % figure
% % % hold on
% % % plot(1:MaxAge, Coverage(1:MaxAge,1), 'b');
% % % plot(1:MaxAge, Coverage(1:MaxAge,2), 'g');
% % % plot(1:MaxAge,Coverage(1:MaxAge,3), 'c');
% % % plot(1:MaxAge,Coverage(1:MaxAge,4), 'r');
% % % legend('in-situ', 'local', 'regional', 'distant');
% % % title(' Treatment Coverage')
% % % 
% % % figure
% % % hold on
% % % plot(1:MaxAge, Coverage(1:MaxAge,5), 'b');
% % % plot(1:MaxAge, Coverage(1:MaxAge,6), 'g');
% % % plot(1:MaxAge,Coverage(1:MaxAge,7), 'c');
% % % plot(1:MaxAge,Coverage(1:MaxAge,8), 'r');
% % % legend('in-situ', 'local', 'regional', 'distant');
% % % title(' Stage Dist')
% % % 
% % % figure
% % % hold on
% % % plot(1:MaxAge, BaseMortality(1,1:MaxAge), 'b');
% % % plot(1:MaxAge, BaseMortality(2,1:MaxAge), 'g');
% % % plot(1:MaxAge,BaseMortality(3,1:MaxAge), 'c');
% % % plot(1:MaxAge,BaseMortality(4,1:MaxAge), 'r');
% % % % plot(1:MaxAge, CancerMort(1:MaxAge,1), 'c+');
% % % legend('in-situ', 'local', 'regional', 'distant');
% % % title(' Cancert mortality by stage')
% % % figure
% % % hold on
% % % plot(1:MaxAge, Mu_Tx(1,1:MaxAge), 'b+');
% % % plot(1:MaxAge, Mu_Tx(2,1:MaxAge), 'g+');
% % % plot(1:MaxAge, Mu_Tx(3,1:MaxAge), 'y+');
% % % plot(1:MaxAge, Mu_Tx(4,1:MaxAge), 'r+');
% % % 
% % % legend('in-situ', 'local', 'regional', 'distant');
% % % plot(1:MaxAge, Mu_NoTx(1,1:MaxAge), 'b');
% % % plot(1:MaxAge, Mu_NoTx(2,1:MaxAge), 'g');
% % % plot(1:MaxAge, Mu_NoTx(3,1:MaxAge), 'y');
% % % plot(1:MaxAge, Mu_NoTx(4,1:MaxAge), 'r');
% % % 
% % % title(' Mu-Treatment(+) and No Treatment(-)')
% % % 
% % % figure
% % % hold on
% % % plot(1:80, CancerMort(1:80,1) - MortalityArray(1,1:80)' , 'c');
% % % plot(1:80, CancerMort(1:80,1), 'c+');
% % % % title(' Mu-Cancer')
% % % % 
% % % % 
% % % % 
% % % % 
% % % % figure
% % % % hold on
% % % %plot(1:MaxAge, MortalityArray(1,1:MaxAge), 'r');
% % % title(' Mortalities')
% % % 
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % 
% % % %Figures from MarkovChainCountry2.m (8 figures)
% % % figure
% % % subplot(4,1,1)
% % % plot(POP(:,1))
% % % title('total pop')
% % % 
% % % subplot(4,1,2)
% % % plot(sum(POP(:,[2:5]),2) )
% % % title('cancer not diagnosed')
% % % 
% % % subplot(4,1,3)
% % % plot(sum(POP(:,[6:9]),2) )
% % % title('total pop')
% % % title('cancer diagnosed')
% % % 
% % % subplot(4,1,4)
% % % plot(OnsetRatesF )
% % % title('total pop')
% % % title('Onset Rate')
% % % 
% % % 
% % % PrevA=Z(1:length(prev))
% % % figure
% % % hold on
% % % plot(AgeArrayMid,prev*100,'k*')
% % % plot(AgeArrayMid,PrevA,'r-')
% % % title('Prevalence of Diagnosed cases by age')
% % % 
% % % figure
% % % hold on
% % % plot(0:100,OnsetRatesF )
% % % 
% % % title('Onset Rate')
% % % 
% % % figure
% % % hold on
% % % plot(0:100,CFR(:,t))
% % % plot(0:100,CFR(:,t-100))
% % % plot(0:100,CFR(:,t-200))
% % % plot(0:100,CFR(:,t-300))
% % % title('CFR')
% % % 
% % % 
% % % 
% % % figure
% % % hold on
% % % X=0:100;
% % % plot(X,diagrT(:,1),'b')
% % % plot(X,diagrT(:,2),'r')
% % % plot(X,diagrT(:,3),'g')
% % % plot(X,diagrT(:,4),'c')
% % % legend('in-situ', 'local', 'regional', 'distant');
% % % title('Time to diagnosis')
% % % 
% % % figure
% % % hold on
% % % X=0:100;
% % % plot(X,DwellTimesF(1,:),'b')
% % % plot(X,DwellTimesF(2,:),'r')
% % % plot(X,DwellTimesF(3,:),'g')
% % % plot(X,DwellTimesF(4,:),'c')
% % % legend('in-situ', 'local', 'regional', 'distant');
% % % 
% % % title('Dwell Time')
% % % 
% % % figure
% % % hold on
% % % plot(1:4,DiagsA(1:4,t),'r-')
% % % plot(1:4,StageDist(1:4),'k*')
% % % title('Stage at diagnosis')
% % % 
% % % 
% % % figure
% % % hold on
% % % subplot(4,1,1)
% % % plot(1:500,DiagsA(2,:),'r-')
% % % title('Local')
% % % 
% % % subplot(4,1,2)
% % % plot(1:500,DiagsA(3,:),'r-')
% % % title('Regional')
% % % 
% % % subplot(4,1,3)
% % % plot(1:500,DiagsA(4,:),'r-')
% % % title('Distant')
% % % 
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % 
% % % %  OUTPUT_DATA = [
% % % 
% % % 
% % % end
