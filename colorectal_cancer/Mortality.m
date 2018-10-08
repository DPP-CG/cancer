function [Mu_Tx, Impact, BaseMortality, TxCoverage, StgDist, RMortalityTx, RMortalityNoTx, ModDwellTime] = Mortality(DwellTime, MortalityArray, Prevalence, StageDist, CancerMort,PAR)
    global CancerMortality 
    global CancerMortalityCalc 
    %global StageDistribution
    MaxAge = 101;

    %StageDistribution = StageDist;
    Healthy = zeros(1,MaxAge);
    Insitu = zeros(1,MaxAge);
    Local = zeros(1,MaxAge);
    Regional = zeros(1,MaxAge);
    Distant = zeros(1,MaxAge);
    StageIV = zeros(1,MaxAge);
    RRInsitu = zeros(1,MaxAge);
    RRLocal = zeros(1,MaxAge);
    RRRegional = zeros(1,MaxAge);
    RRDistant = zeros(1,MaxAge);
    RRStageIV = zeros(1,MaxAge);

    global Mu_TxVal 
    Mu_TxVal = zeros(5,1);
    Mu_Tx = zeros(5,MaxAge);
    
    global Mu_NoTxVal 
    Mu_NoTxVal = zeros(5,1);
    Mu_NoTx = zeros(5,MaxAge);
    Coverage = zeros(MaxAge,10);

    DwellTime(5,:)=0;
    ModDwellTime = DwellTime;

    RMortalityTx = zeros(5,MaxAge);
    RMortalityNoTx = zeros(5,MaxAge);
    t = 10;
    
    for i = 1:MaxAge - 10

       DwellRate(1,i:i+9)= 1 ./ DwellTime(1,i:i+9);
       DwellRate(2,i:i+9)= 1 ./ DwellTime(2,i:i+9);
       DwellRate(3,i:i+9)= 1 ./ DwellTime(3,i:i+9);
       DwellRate(4,i:i+9)= 1 ./ DwellTime(4,i:i+9);

       if StageDist(1,5) == 0
           DwellRate(5,i:i+9) = 0;
       else
           DwellRate(5,i:i+9)= 1 ./ DwellTime(5,i:i+9);
       end

        % Ten-year survival without treatment
       if StageDist(1,5) == 0
            stg = 4;
       else
            stg = 5;
       end

       Healthy(1,i) = TenYearSurvival(MortalityArray(1,i:i+9),zeros(1,10));
       Insitu(1,i) = TenYearSurvival(DwellRate(1:stg,i:i+9),MortalityArray(1,i:i+9));
       Local(1,i) = TenYearSurvival(DwellRate(2:stg,i:i+9),MortalityArray(1,i:i+9));
       Regional(1,i) = TenYearSurvival(DwellRate(3:stg,i:i+9),MortalityArray(1,i:i+9));
       Distant(1,i) = TenYearSurvival(DwellRate(4:stg,i:i+9),MortalityArray(1,i:i+9));
       StageIV(1,i) = TenYearSurvival(DwellRate(5:5,i:i+9),MortalityArray(1,i:i+9));

    % Ten-year RELATIVE survival without treatment   
       Insitu(1,i) = max(Insitu(1,i),Local(1,i));% very high rate for in-situ in older age

       RRInsitu(1,i) = Insitu(1,i) / Healthy(1,i);
       RRLocal(1,i) = Local(1,i) / Healthy(1,i);
       RRRegional(1,i) = Regional(1,i) / Healthy(1,i);
       RRDistant(1,i) = Distant(1,i) / Healthy(1,i);
       if StageDist(1,5) == 0
           RRStageIV(1,i) = 0;
       else
           RRStageIV(1,i) = StageIV(1,i) / Healthy(1,i);
       end

      %Yearly mortality rate without treatment 
       Mu_NoTx(1,i) = (MortalityArray(1,i)*t - log(RRInsitu(1,i)))/t;
       Mu_NoTx(2,i) = (MortalityArray(1,i)*t - log(RRLocal(1,i)))/t;
       Mu_NoTx(3,i) = (MortalityArray(1,i)*t - log(RRRegional(1,i)))/t;
       Mu_NoTx(4,i) = (MortalityArray(1,i)*t - log(RRDistant(1,i)))/t;
       if StageDist(1,5) == 0
          Mu_NoTx(5,i) = 0;
       else
        Mu_NoTx(5,i) = (MortalityArray(1,i)*t - log(RRStageIV(1,i)))/t;
       end


    end

    for i = MaxAge - 10 : MaxAge
       RRInsitu(1,i) =  RRInsitu(1,i - 1); 
       RRLocal(1,i) = RRLocal(1,i - 1);
       RRRegional(1,i) = RRRegional(1,i - 1);
       RRDistant(1,i) = RRDistant(1,i - 1);
       RRStageIV(1,i) = RRStageIV(1,i - 1);

       Mu_NoTx(1,i) = Mu_NoTx(1,i-1) ;
       Mu_NoTx(2,i) = Mu_NoTx(2,i-1) ;
       Mu_NoTx(3,i) = Mu_NoTx(3,i-1) ;
       Mu_NoTx(4,i) = Mu_NoTx(4,i-1) ;
       Mu_NoTx(5,i) = Mu_NoTx(5,i-1) ;
    end

    'Ten year'
    %%IF relative survival is available
    % InsituCF = mean(RRInsitu);%(1,30:80));
    % LocalCF = mean(RRLocal);%(1,30:80));
    % RegionalCF = mean(RRRegional);%(1,30:80));
    % DistantCF = mean(RRDistant);%(1,30:80));
    % StageIVCF = mean(RRStageIV);%(1,30:80));

    
    %MuTxAll = [0 0.002 0.03 0.55 0]; %Leshno et al., 2003, converted to rate from annual CRC mortality in Epidemilogy file.
    %Above Commented by Prashant 12 July
    %Using below Mortality with treatment as evaluated from the cancer
    %assumptions file ~\Dropbox\Cancer manuscripts\Breast and CRC_ cancer impact model assumptions - 7-11-2016.xlsx
    MuTxAll = [0.0107 0.0107 0.052 0.5660 0];


    for i = 2:5
        MuRatio(i) = (sum(Mu_NoTx(i,:).*Prevalence)/sum(Prevalence)) / MuTxAll(i);
        Mu_Tx(i,:) = Mu_NoTx(i,:)./MuRatio(i);
    end


    Mu_Tx(1,:) = Mu_NoTx(1,:).*Mu_Tx(2,:)./Mu_NoTx(2,:);
    for i = 1:MaxAge
            if StageDist(1,5) == 0
               Mu_Tx(5,i) = 0;
            end

         % 10-year Relative surival with_Tx

            RI(1,i) =  min(1,exp(MortalityArray(1,i)*t- Mu_Tx(1,i)*t));
            RL(1,i) =  min(1,exp(MortalityArray(1,i)*t- Mu_Tx(2,i)*t));
            RR(1,i) =  min(1,exp(MortalityArray(1,i)*t- Mu_Tx(3,i)*t));
            RD(1,i) =  min(1,exp(MortalityArray(1,i)*t- Mu_Tx(4,i)*t));
            RIV(1,i) =  min(1,exp(MortalityArray(1,i)*t- Mu_Tx(5,i)*t));


        CancerMortality = CancerMort(i,1) + MortalityArray(1,i);
        Mu_TxVal(:,1) = Mu_Tx(:,i);
        Mu_NoTxVal(:,1) = Mu_NoTx(:,i);

        if StageDist(1,1) == 0
            uVal = 0;
        else
            uVal = 1;
        end

        if StageDist(1,5) == 0
            uVal5 = 0;
        else
            uVal5 = 1;
        end
        % Assuming Tx coverage is 100% 
        LowerBound = [uVal 1 1 1 uVal5 0 0 0 0 0 ];
        UpperBound = [uVal 1 1 1 uVal5 uVal 1 1 1 uVal5 ];
        options = optimoptions(@fmincon,'MaxIter', 1000, 'TolX', 0.001);
        %Covergae: 1:4= coverage of Tx in  stage insitu, local, regional, and distant; 4:6 =
        %Stage distribution
        [Coverage(i, 1:10), error(i)] = fmincon(@FuncCancerMort, [0 1 1 1 0 0 0.3 0.4 .3 0], [], [], [0 0 0 0 0 1 1 1 1 1], [1],LowerBound,UpperBound,[] );
        %Coverage(6,i) = 1- Coverage(5,i) - Coverage (4,i);
        Cancer(1,i) = CancerMortalityCalc;

        for j = 1: 10
        if Coverage(i, j) < 0.01
           Coverage(i, j) = 0; 
        end
        end
    end
        %Subtracting disease free mortality from cancer mortality as it is added
        %back in Spectrum
    %   

        Mu_Tx(1,:) = max(0,Mu_Tx(1,:) - MortalityArray(1,:));
        Mu_Tx(2,:) = max(0,Mu_Tx(2,:) - MortalityArray(1,:));
        Mu_Tx(3,:) = max(0,Mu_Tx(3,:) - MortalityArray(1,:));
        Mu_Tx(4,:) = max(0,Mu_Tx(4,:) - MortalityArray(1,:));
        Mu_Tx(5,:) = max(0,Mu_Tx(5,:) - MortalityArray(1,:));

        Mu_NoTx(1,:) = max(0,Mu_NoTx(1,:) - MortalityArray(1,:));
        Mu_NoTx(2,:) = max(0,Mu_NoTx(2,:) - MortalityArray(1,:));
        Mu_NoTx(3,:) = max(0,Mu_NoTx(3,:) - MortalityArray(1,:));
        Mu_NoTx(4,:) = max(0,Mu_NoTx(4,:) - MortalityArray(1,:));
        Mu_NoTx(5,:) = max(0,Mu_NoTx(5,:) - MortalityArray(1,:));

        Mu_Tx(1,80:end) = Mu_Tx(1,80);
        Mu_Tx(2,80:end) = Mu_Tx(2,80);
        Mu_Tx(3,80:end) = Mu_Tx(3,80);
        Mu_Tx(4,80:end) = Mu_Tx(4,80);
        Mu_Tx(5,80:end) = Mu_Tx(5,80);

        Mu_NoTx(1,80:end) = Mu_NoTx(1,80);
        Mu_NoTx(2,80:end) = Mu_NoTx(2,80);
        Mu_NoTx(3,80:end) = Mu_NoTx(3,80);
        Mu_NoTx(4,80:end) = Mu_NoTx(4,80);
        Mu_NoTx(5,80:end) = Mu_NoTx(5,80);

        RMortalityTx(1,:)= Mu_Tx(1,:)./ MortalityArray(1,:);
        RMortalityTx(2,:)= Mu_Tx(2,:)./ MortalityArray(1,:);
        RMortalityTx(3,:)= Mu_Tx(3,:)./ MortalityArray(1,:);
        RMortalityTx(4,:)= Mu_Tx(4,:)./ MortalityArray(1,:);
        RMortalityTx(5,:)= Mu_Tx(5,:)./ MortalityArray(1,:);

        RMortalityNoTx(1,:)= Mu_NoTx(1,:)./ MortalityArray(1,:);
        RMortalityNoTx(2,:)= Mu_NoTx(2,:)./ MortalityArray(1,:);
        RMortalityNoTx(3,:)= Mu_NoTx(3,:)./ MortalityArray(1,:);
        RMortalityNoTx(4,:)= Mu_NoTx(4,:)./ MortalityArray(1,:); 
        RMortalityNoTx(5,:)= Mu_NoTx(5,:)./ MortalityArray(1,:); 


    BaseMortality = zeros(size(Mu_Tx,1), size(Mu_Tx,2));
    BaseMortality(1,:) = (Coverage(:,1)' .*  Mu_Tx(1,:) + (1- Coverage(:,1)').* Mu_NoTx(1,:))  ;
    BaseMortality(2,:) = (Coverage(:,2)' .*  Mu_Tx(2,:) + (1- Coverage(:,2)').* Mu_NoTx(2,:))  ;
    BaseMortality(3,:) = (Coverage(:,3)' .*  Mu_Tx(3,:) + (1- Coverage(:,3)').* Mu_NoTx(3,:))  ;
    BaseMortality(4,:) = (Coverage(:,4)' .*  Mu_Tx(4,:) + (1- Coverage(:,4)').* Mu_NoTx(4,:)) ;
    BaseMortality(5,:) = (Coverage(:,5)' .*  Mu_Tx(5,:) + (1- Coverage(:,5)').* Mu_NoTx(5,:)) ;

    % estimating increase in weighted mortality for unit change in coverage
    Impact (1,:) = (Mu_Tx(1,:) - Mu_NoTx(1,:)) ./ BaseMortality(1,:);
    Impact (2,:) = (Mu_Tx(2,:) - Mu_NoTx(2,:)) ./ BaseMortality(2,:);
    Impact (3,:) = (Mu_Tx(3,:) - Mu_NoTx(3,:)) ./ BaseMortality(3,:);
    Impact (4,:) = (Mu_Tx(4,:) - Mu_NoTx(4,:)) ./ BaseMortality(4,:);
    Impact (5,:) = (Mu_Tx(5,:) - Mu_NoTx(5,:)) ./ BaseMortality(5,:);

    TxCoverage(:,:) = max(0,Coverage(:,1:5));
    StgDist(:,:) =  max(0,Coverage(:,6:10));
    close all
    figure
    hold on
    plot(1:MaxAge, Healthy(1,1:MaxAge), 'r');
    plot(1:MaxAge, Insitu(1,1:MaxAge), 'y');
    plot(1:MaxAge, Local(1,1:MaxAge), 'b');
    plot(1:MaxAge, Regional(1,1:MaxAge), 'g');
    plot(1:MaxAge, Distant(1,1:MaxAge), 'c');
    %plot(1:MaxAge, StageIV(1,1:MaxAge), 'm');
    legend('in-situ', 'local', 'regional', 'distant');

    % plot(30:80, T_Local(1,30:80), 'b-');
    % plot(30:80, T_Regional(1,30:80), 'g-');
    % plot(30:80, T_Distant(1,30:80), 'c-');
    title('5-year Survival')


    figure
    hold on
    plot(1:MaxAge, RMortalityNoTx(1,1:MaxAge), 'y');
    plot(1:MaxAge, RMortalityNoTx(2,1:MaxAge), 'b');
    plot(1:MaxAge, RMortalityNoTx(3,1:MaxAge), 'g');
    plot(1:MaxAge,RMortalityNoTx(4,1:MaxAge), 'c');
    %plot(1:MaxAge,RMortalityNoTx(5,1:MaxAge), 'm');

    legend('in-situ', 'local', 'regional', 'distant');
    plot(1:MaxAge, RMortalityTx(1,1:MaxAge), 'y+');
    plot(1:MaxAge, RMortalityTx(2,1:MaxAge), 'b+');
    plot(1:MaxAge, RMortalityTx(3,1:MaxAge), 'g+');
    plot(1:MaxAge, RMortalityTx(4,1:MaxAge), 'c+');
    %plot(1:MaxAge, RMortalityTx(5,1:MaxAge), 'm+');

    title('yearly relative mortality in cancer compared to non-cancer persons ("-"without Tx "+" with Tx)')

    figure
    hold on
    plot(1:MaxAge, RRInsitu(1,1:MaxAge), 'y');
    plot(1:MaxAge, RRLocal(1,1:MaxAge), 'b');
    plot(1:MaxAge, RRRegional(1,1:MaxAge), 'g');
    plot(1:MaxAge,RRDistant(1,1:MaxAge), 'c');
    %plot(1:MaxAge,RRStageIV(1,1:MaxAge), 'm');

    legend('in-situ', 'local', 'regional', 'distant');
    plot(1:MaxAge, RI(1,1:MaxAge), 'y+');
    plot(1:MaxAge, RL(1,1:MaxAge), 'b+');
    plot(1:MaxAge, RR(1,1:MaxAge), 'g+');
    plot(1:MaxAge, RD(1,1:MaxAge), 'c+');
    %plot(1:MaxAge, RD(1,1:MaxAge), 'm+');
    title('10-year relative survival in cancer compared to non-cancer persons ("-"without Tx "+" with Tx)')

    figure
    hold on
    plot(1:MaxAge, Coverage(1:MaxAge,1), 'b');
    plot(1:MaxAge, Coverage(1:MaxAge,2), 'g');
    plot(1:MaxAge,Coverage(1:MaxAge,3), 'c');
    plot(1:MaxAge,Coverage(1:MaxAge,4), 'r');
    %plot(1:MaxAge,Coverage(1:MaxAge,5), 'm');
    legend('in-situ', 'local', 'regional', 'distant');

    title(' Treatment Coverage')

    figure
    hold on
    plot(1:MaxAge, Coverage(1:MaxAge,6), 'b');
    plot(1:MaxAge, Coverage(1:MaxAge,7), 'g');
    plot(1:MaxAge,Coverage(1:MaxAge,8), 'c');
    plot(1:MaxAge,Coverage(1:MaxAge,9), 'r');
    %plot(1:MaxAge,Coverage(1:MaxAge,10), 'm');
    legend('in-situ', 'local', 'regional', 'distant');

    title('Not Relevant: Stage Dist')

    figure
    hold on
    plot(1:MaxAge, BaseMortality(1,1:MaxAge), 'b');
    plot(1:MaxAge, BaseMortality(2,1:MaxAge), 'g');
    plot(1:MaxAge,BaseMortality(3,1:MaxAge), 'c');
    plot(1:MaxAge,BaseMortality(4,1:MaxAge), 'y');
    %plot(1:MaxAge,BaseMortality(5,1:MaxAge), 'm');
    legend('in-situ', 'local', 'regional', 'distant');

    % plot(1:MaxAge,BaseMortality(4,1:MaxAge), 'r');
    % plot(1:MaxAge, CancerMort(1:MaxAge,1), 'c+');

    title('Not Relevant: Cancer mortality by stage')
    figure
    hold on
    plot(1:MaxAge, Mu_Tx(1,1:MaxAge), 'b+');
    plot(1:MaxAge, Mu_Tx(2,1:MaxAge), 'g+');
    plot(1:MaxAge, Mu_Tx(3,1:MaxAge), 'y+');
    plot(1:MaxAge, Mu_Tx(4,1:MaxAge), 'r+');
    %plot(1:MaxAge, Mu_Tx(5,1:MaxAge), 'm+');
    legend('in-situ', 'local', 'regional', 'distant');

    plot(1:MaxAge, Mu_NoTx(1,1:MaxAge), 'b');
    plot(1:MaxAge, Mu_NoTx(2,1:MaxAge), 'g');
    plot(1:MaxAge, Mu_NoTx(3,1:MaxAge), 'y');
    plot(1:MaxAge, Mu_NoTx(4,1:MaxAge), 'r');
    %plot(1:MaxAge, Mu_NoTx(5,1:MaxAge), 'm');

    title('Very Relevant: Relative Mortality compared to disease free: with Treatment (+) & without Treatment (-)')

    figure
    hold on
    plot(1:80, CancerMort(1:80,1) - MortalityArray(1,1:80)' , 'c');
    plot(1:80, CancerMort(1:80,1), 'c+');

    title(' Mortalities')


      'Mort Results Written'

    %pause
end

function Error = FuncCancerMort(coverage)
global Mu_TxVal 
global Mu_NoTxVal
%global StageDistribution 
global CancerMortality
global CancerMortalityCalc
CancerMortalityCalc =  (coverage(1,1) *  Mu_TxVal(1,1) + (1- coverage(1,1))* Mu_NoTxVal(1,1)) * coverage(1,6) + ...
                       (coverage(1,2) *  Mu_TxVal(2,1) + (1- coverage(1,2))* Mu_NoTxVal(2,1)) * coverage(1,7) + ...
                       (coverage(1,3) *  Mu_TxVal(3,1) + (1- coverage(1,3))* Mu_NoTxVal(3,1)) * coverage(1,8) + ...
                       (coverage(1,4) *  Mu_TxVal(4,1) + (1- coverage(1,4))* Mu_NoTxVal(4,1)) * coverage(1,9) + ...
                       (coverage(1,5) *  Mu_TxVal(5,1) + (1- coverage(1,5))* Mu_NoTxVal(5,1)) * coverage(1,10);
       
Error = CancerMortality - CancerMortalityCalc;
Error = Error * Error;
end


function Proportion = TenYearSurvival(RatesArray, Mu)
    dt =100;
    RatesArray = RatesArray /dt;
    Mu = Mu / dt;
    N = zeros (1,size(RatesArray,1));
    N(1,1) = 100000;
    for t = 1:10 * dt
        N(1,1) = N(1,1) - N(1,1) * (RatesArray(1,max(1,round(t/dt)))+ Mu(1,max(1,round(t/dt))));
        for i = size(RatesArray,1):-1:2
        N(1,i)= N(1,i) + N(1,i-1) * RatesArray(i-1,max(1,round(t/dt))) - N(1,i) * (RatesArray(i, max(1,round(t/dt))) + Mu(1, max(1,round(t/dt))));
        end
    end


    Proportion = sum(sum(N))/100000;
end



