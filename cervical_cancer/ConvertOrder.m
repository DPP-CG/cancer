   Output2=zeros(748,38);
    Output2(:)=NaN;
    Output2=num2cell(Output2);
    [~,~,Output1]=xlsread('Temp Output.xlsx');
    Output2(1:133,1:38)=Output1(1:133,1:38);
    Output2(244:748,1:38)=Output1(134:638,1:38);
    Output2(134:242,1:38)=Output1(640:748,1:38);
    xlswrite('Output for Spectrum.xlsx',Output2);
    csvwrite('CervicalCancer.csv',Output2,748,38);