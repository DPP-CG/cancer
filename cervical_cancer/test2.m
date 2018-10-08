     Output2=zeros(749,52);
     Output2(:)=NaN;
     Output2=num2cell(Output2);
     [~,~,Output1]=xlsread('Temp Output.xlsx');
     Output2(1:133,1:52)=Output1(1:133,1:52);
     Output2(244:749,1:52)=Output1(134:639,1:52);
     Output2(134:242,1:52)=Output1(640:748,1:52);
     Output2(750:880,1:52)=Output1(750:880,1:52);
     xlswrite('Output for Spectrum.xlsx',Output2);