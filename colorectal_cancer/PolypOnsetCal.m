function [InsituOnsetRateAllAge,polypOnsetRateAgeGroup, PrevPolys] = PolypOnsetCal(OnsetRatesF,AgeArrayLower,AgeArrayUpper,PAR)

%X =[1:101];
X =[20:80];%CGJan23_16- fitting to only 20:80 ages 
OnsetRate=OnsetRatesF;
OnsetRate(1:20)=0;
% OnsetRate(80:end)=OnsetRate(79);    %assuming that we have constant or stable onset beyond age 80 instead of dip indicated by interpolation

Mark = 2;

for i = 1 : 1%CGJan23_16- using exponential fit. Tried different starting points, all generate same solution indicating global optimal
    if Mark == 1  %Step function
        
        betaMatrix = zeros(12, 5);
        
        beta0 = ones(1,12);
        r = 0.1*rand(1);
        beta0 = beta0*r;
        
        beta0(1,1:4)=0;
        LB = [0 0 0 0 0 0 0 0 0 0 0 0];%zeros(1,size(AgeArrayLower,2)+1);
        UB = [0.00000000000001 0.00000000000001 0.00000000000001 0.00000000000001 1 1 1 1 1 1 1 1];%ones(1,size(AgeArrayLower,2)+1);
        
        %         %OPTIONS = optimset('largescale','off','LevenbergMarquardt','on','display','off','Diagnostics','on','MaxFunEvals',700,'TolFun',1e-8,'TolCon',1e-8);
        %         OPTIONS = optimset('MaxFunEvals',2000,'Diagnostics','on','TolFun',1e-6,'TolCon',1e-6);
        %         Y =OnsetRate';
        %         FitArray = OnsetRate';
        %
        %         %finding the best beta fit using lsqcurvefit
        %         [beta,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@(x,X)...
        %         PopSim(x,X,FitArray,AgeArrayLower,AgeArrayUpper,PAR,0, Mark),...
        %         beta0,X,Y,LB,UB,OPTIONS);
        %
        
        
    elseif Mark == 2    %Exponential function
        
        beta0 = ones(1, 2);
        r = 0.01*rand(1);
        beta0 = beta0*r;
        LB = [0 0 ];%zeros(1,size(AgeArrayLower,2)+1);
        UB = [0.1 0.1];%ones(1,size(AgeArrayLower,2)+1);
        
    elseif Mark ==3    %Logistic function
        beta0 = [1 3 0];
        
        LB = [-100 -100 -1 ];%zeros(1,size(AgeArrayLower,2)+1);
        UB = [100 100 1];%ones(1,size(AgeArrayLower,2)+1);
        
    elseif Mark == 4    %Logistic function2
        beta0 = [50 0.1 0];
        % beta0 = beta0*0.7;
        LB = [40 0 0 ];%zeros(1,size(AgeArrayLower,2)+1);
        UB = [100 1 1];%ones(1,size(AgeArrayLower,2)+1);
        
    elseif Mark == 5    %Polynomial Function
        beta0 = [0.5 0.4 0.3 0.2 0.1];
        % beta0 = beta0*0.7;
        LB = [-20 -20 -20 -20 -20];%zeros(1,size(AgeArrayLower,2)+1);
        UB = [20 20 20 20 20];%ones(1,size(AgeArrayLower,2)+1);
        Y =[OnsetRate];
        beta0 = polyfit(X, Y, 4);
        beta0;
        
    elseif Mark == 6    %Logistic function2
        beta0 = [1 1 1 1 1 1];
        % beta0 = beta0*0.7;
        LB = [-1 -1 -1 -1 -1 -1 ];%zeros(1,size(AgeArrayLower,2)+1);
        UB = [100 100 100 100 100 100];%ones(1,size(AgeArrayLower,2)+1);
    end
    
    
    %OPTIONS = optimset('largescale','off','LevenbergMarquardt','on','display','off','Diagnostics','on','MaxFunEvals',700,'TolFun',1e-8,'TolCon',1e-8);
    OPTIONS = optimset('MaxFunEvals',2000,'Diagnostics','on','TolFun',1e-6,'TolCon',1e-6);
    Y =[OnsetRate(20:80)]';%CGJan23_16- fitting to only 20:80 ages 
    FitArray = [OnsetRate(20:80)]';%CGJan23_16- fitting to only 20:80 ages 
    
    %finding the best beta fit using lsqcurvefit
    [beta,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@(x,X)...
        PopSim(x,X,FitArray,AgeArrayLower,AgeArrayUpper,PAR,0, Mark),...
        beta0,X,Y,LB,UB,OPTIONS);
    
    %    [beta,resnorm,residual,exitflag,output,lambda,jacobian] = fmincon(@(x,X)PopSim(x,X,FitArray,AgeArrayLower,AgeArrayUpper,PAR,0, Mark, Y),beta0,[],[],[],[],@(Z)polyval(-x, X));
    
    %when it at its the best beta value, use the beta in PopSim. We know this
    %by the 1 or 0 at the end
    [InsituOnsetRate,polypOnsetRateAgeGroup, PrevPolys]=PopSim(beta,X,OnsetRate,AgeArrayLower,AgeArrayUpper,PAR,1, Mark);
   InsituOnsetRateAllAge(1:101)=0;
   InsituOnsetRateAllAge(20:80) = InsituOnsetRate;%CGJan23_16- compensating for fitting to only 20:80 ages 
   InsituOnsetRateAllAge(1:20)=0;InsituOnsetRateAllAge(80:101)=InsituOnsetRateAllAge(80);%CGJan23_16- compensating fitting to only 20:80 ages 
    
   betaMatrix(:, i) = polypOnsetRateAgeGroup(:,1);
    
end

figure
hold on
plot (1:12, betaMatrix(:, 1), 'r');

end


function [Result,polypOnsetRateAgeGroup, PrevPolys] = PopSim(beta,~,OnsetRate,AgeArrayLower,AgeArrayUpper,PAR,plotInd, Mark)

    AgeGroupL=[1 5 10 14 20 25 30 40 50 60 70 80];
    AgeGroupU=[4 9 15 19 24 29 39 49 59 69 79 100];

    probAgeA=PAR{1};probAgeA(81:end)=probAgeA(81)/length(probAgeA(81:end));
    mu=PAR{2}';mu(81:end)=mu(80)+1/5;
    BirthsF=PAR{3};

    agemax=101;
    X =1:101;
    m_a=length(probAgeA);

    pop0=100000;%work with pop of standard size

    Stages = 5;
    POP=zeros(m_a,Stages);
    POP(1:m_a,1)=probAgeA*pop0;

    t_max=500;

    PolyponsetRate = zeros(5, m_a);

    AgeArrayM=(AgeGroupU+AgeGroupL)/2;


    if Mark == 1  %Step function 
        %beta(5:12)=smooth(beta(5:12));
        PolyponsetRate(1,:)= interp1(AgeArrayM,beta(1:12),0:100,'PCHIP');
        PolyponsetRate(1,1:20) = 0;
        PolyponsetRate=PolyponsetRate';

    elseif Mark == 2   %Exponential function 

        PolyponsetRate(1,:)= beta(1)*(exp(beta(2)*X));
        PolyponsetRate(1,1:20) = 0;   
        PolyponsetRate=PolyponsetRate';

    elseif Mark ==3    %Logistic function 
        PolyponsetRate(1,:)= beta(3)./(1+(beta(1)*(exp(-beta(2)*X))));
        PolyponsetRate(1,1:20) = 0;
        PolyponsetRate=PolyponsetRate';
        
    elseif Mark ==4   %Logistic function2
        PolyponsetRate(1,20:end)= beta(3)./(1+(exp(-beta(2)*(X(20:end)-beta(1)))));
        PolyponsetRate(1,1:20) = 0;
        PolyponsetRate=PolyponsetRate';
        
    elseif Mark ==5   %Polynomial function
        %PolyponsetRate(1,20:end)= (beta(4).*(X(20:end).^3)) +(beta(3).*(X(20:end).^2)) + (beta(2).*(X(20:end).^1)) + beta(1);
        %PolyponsetRate(1,20:end)= (beta(1).*(X(20:end).^4)) +(beta(2).*(X(20:end).^3)) + (beta(3).*(X(20:end).^2))+ (beta(4).*(X(20:end).^1)) + beta(5);
        PolyponsetRate(1,:)= polyval(beta,X);
        %PolyponsetRate(1,1:20) = 0;
        PolyponsetRate=PolyponsetRate';
        
    elseif Mark == 6   %Generalized Logistic function
        %PolyponsetRate(1,20:end)= (beta(4).*(X(20:end).^3)) +(beta(3).*(X(20:end).^2)) + (beta(2).*(X(20:end).^1)) + beta(1);
        %PolyponsetRate(1,20:end)= (beta(1).*(X(20:end).^4)) +(beta(2).*(X(20:end).^3)) + (beta(3).*(X(20:end).^2))+ (beta(4).*(X(20:end).^1)) + beta(5);
        for i =1:101
             PolyponsetRate(1,i)=  beta(1) + ((beta(2)-beta(1)) /( beta(3)+beta(4)* exp(-beta(5)*X(i)))^(1./beta(6)));
        end
       % PolyponsetRate(1,1:20) = 0;
        PolyponsetRate=PolyponsetRate';        
    end 


    % temporary progression rate from Chen et al., 2003
    q1(1:101,1) = 0.021;
    q2(1:101,1) = 0.057;
    q3(1:101,1) = 0.063;
    q4(1:101,1) = 0; %initializing. q4 will be calculated below

    Dist_polyp=0.23; % Average of Chen Loeve and Goto 0.18 and 0.27
    % q4 =Dist_polyp*PolyponsetRate;  on Jan 9 we understand that the
    % Dist_Polyp = 0.23 is not for the polyp onset but for the Insitu onset
    % from the Chen paper, Loeve paper and Goto paper.


    polypPeople = zeros(4, t_max); %Prashant: added for steady state check 9 Jan 2016
    dt=1/15;
    for t=1:t_max
        value=zeros(101,1);
        for t1=1:1/dt

            %POP(:,5)=POP(:,5) + dt*( q3.*POP(:,4) + q4.*POP(:,1));                      %Insitu Cancer
            POP(:,4)=POP(:,4) + dt*((-mu - q3).*POP(:,4) + q2.*POP(:,3));               %Polyp stage 3 (edit to get the actual polyp size group)
            POP(:,3)=POP(:,3) + dt*((-mu - q2).*POP(:,3) + q1.*POP(:,2));               %Polyp Stage 2 (edit to get the actual polyp size group)
            POP(:,2)=POP(:,2) + dt*((-mu - q1).*POP(:,2) +(PolyponsetRate(:,1)).*POP(:,1));  %Polyp Stage 1 (edit to get the actual polyp size group)
            POP(:,1)=POP(:,1) + dt*(-mu-PolyponsetRate(:,1)).*POP(:,1);   %Disease Free

            value=value+dt*( q3.*POP(:,4));
        end
        value = value / (1 - Dist_polyp); %Making the overall insitu onset (including onset from non polyp pathway) and not just the insitu onset from adenoma polys
        q4 = (value * Dist_polyp) ./ POP(:,1) ;
        InsituOnsetRate = value./sum(POP,2); %CAN WE HAVE THIS AND ABOVE 2 LINES OUTSIDE THE FOR LOOP?

        for I=1:Stages
            POP(2:end,I)=POP(1:end-1,I);
            POP(1,I)= 0;
        end
        %POP(end,4) = POP(end,4) +lAge;
        %must decide if births are to be based on population size
        Births = probAgeA(1,1)* sum(sum(POP));%42;%InsituOnsetRate.*probAgeA' *pop0;%BirthsF*sum(POP(:));
        POP(1,1)=Births;
        %Prashant: added polypPeople to check steady state of proportion of
        %people in various stages of polyp for 500 years simulation
        polypPeople(1,t)=sum(POP(:,1))/ sum(sum(POP(:,:)));
        polypPeople(2,t)=sum(POP(:,2))/ sum(sum(POP(:,:)));
        polypPeople(3,t)=sum(POP(:,3))/ sum(sum(POP(:,:)));
        polypPeople(4,t)=sum(POP(:,4))/ sum(sum(POP(:,:)));

        %Result = InsituOnsetRate(1:101,1);%CAN WE HAVE THIS OUTSIDE THE FOR LOOP?
        Result = InsituOnsetRate(20:80,1);
    end
    if plotInd==1

        %Added by Prashant to return all the transitions from Healthy to Polyp <=5mm
        %to Polyp >=6mm to <10mm to Polyp >10mm to Insitu Cancer and from Healthy
        %to Insitu Cancer
        polypOnsetRateAgeGroup = zeros(12,5);
       % polypOnsetRateAgeGroup(:,1) = beta; 
        polypOnsetRateAgeGroup(:,2) = q1(1,1);
        polypOnsetRateAgeGroup(:,3) = q2(1,1);
        polypOnsetRateAgeGroup(:,4) = q3(1,1);

        %Calculating proportion of prevalance of Polyps as PrevG
        PrevPolys = zeros(size(AgeArrayLower,2),3);
        TotalPop=sum(POP,2);
        for i = 1 : 12        
            indB=AgeArrayLower(i):AgeArrayUpper(i);
            PrevPolys(i,1)=sum(POP(indB,2))./sum(TotalPop(indB,1));
            PrevPolys(i,2)=sum(POP(indB,3))./sum(TotalPop(indB,1));
            PrevPolys(i,3)=sum(POP(indB,4))./sum(TotalPop(indB,1));
            polypOnsetRateAgeGroup(i,1) = sum(PolyponsetRate(indB,1).*probAgeA(indB)')/sum(probAgeA(indB));
            polypOnsetRateAgeGroup(i,5) = sum(q4(indB,1).*probAgeA(indB)')/sum(probAgeA(indB));
        end
        close all

        figure
        hold on
        plot(0:100, PolyponsetRate(:,1),'r');
        plot(0:100, q4(:,1), 'b');
        title('Polyp onset rates');

        figure
        hold on
        plot(0:100, [InsituOnsetRate(:,1)]','r');
        plot(0:100, OnsetRate(1,:),'r*');
        title('Insitu rates');

        %Prashant: adding PolypOnset graph plot
        figure
        hold on
        %plot(1:500, polypPeople(1,:),'r');
        plot(1:500, polypPeople(2,:),'b*');
        plot(1:500, polypPeople(3,:),'go');
        plot(1:500, polypPeople(4,:),'r-');
        title('Polyp distribution 500 Yrs');
        
        %%
        
        [word_format_prevalence, prev_val] = get_word_format_prevalance(PrevPolys,PAR,AgeArrayLower,AgeArrayUpper);
        %%

    end
end