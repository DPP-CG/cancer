POP=zeros(80,10);
POP(:,1)=100000;
r1=0.15;
r2=0.0931; %Highrisk HPV to LSIL
r3=0.2107; %LSIL to HSIL
r4=0.1188; %LSIL to Highrisk HPV
r5=0.01715; %HSIL to Highrisk HPV
r6=0.01715; %HSIL to Normal
r7=0.1188; %LSIL to Normal
r8=0.05; 
r9=0.0568; %Lowrisk HPV to LSIL
r10=0.0921; %LSIL to HSIL
r11=0.0704; %HSIL to Lowrisk HPV
r12=0.1059; %LSIL to Lowrisk HPV
r13=0.0704; %HSIL to Normal
r14=0.1059; %LSIL to Normal
r17=0.2693; %HighRisk HPV to Normal
r18=0.2682; %LowRisk HPV to Normal
Immunity = 10;%beta(end);
r19=1./Immunity; %Immunity to Normal
r20=0.2688;


t_max=500;

r15=zeros(size(POP,1),size(POP,1));%HSIL of Highrisk HPV to CIS, varied by ages
for i=1:30
    r15(i,i)=0.0292;
end
for i=31:39
    r15(i,i)=0.0506;
end
for i=40:49
    r15(i,i)=0.1344;
end
for i=50:size(POP,1)
    r15(i,i)=0.1952;
end

r16=zeros(size(POP,1),size(POP,1));%HSIL of Lowrisk HPV to CIS, varied by ages
for i=1:30
    r16(i,i)=0.0070;
end
for i=31:39
    r16(i,i)=0.0140;
end
for i=40:49
    r16(i,i)=0.0221;
end
for i=50:size(POP,1)
    r16(i,i)=0.0445;
end

dt=1/12;
for t=1:t_max

value = 0;
  for t1=1:1/dt
         
      
        POP(:,9)=POP(:,9) + (-r19.*POP(:,9) + r7.*POP(:,3) + r6.*POP(:,4) + r14.*POP(:,6) + r13.*POP(:,7) + r20.*POP(:,10) + r17.*POP(:,2) + r18.*POP(:,5))*dt;

        POP(:,10)=POP(:,10) + (+ r4.*POP(:,3) + r5.*POP(:,4) + r12.*POP(:,6) + r11.*POP(:,7)- r20.*POP(:,10))*dt;

        POP(:,8)=POP(:,8) + (+r15*POP(:,4) + r16*POP(:,7))*dt;
        value=value + (r15*POP(:,4) + r16*POP(:,7))*dt; 
            
        
        POP(:,4)=POP(:,4) + (-r5.*POP(:,4) - r6.*POP(:,4) - r15*POP(:,4) + r3.*POP(:,3))*dt;        
        POP(:,3)=POP(:,3) + (-r7.*POP(:,3) - r4.*POP(:,3) - r3.*POP(:,3) + r2.*POP(:,2))*dt;
        POP(:,2)=POP(:,2) + (-r2.*POP(:,2) + r1'.*POP(:,1) - r17.*POP(:,2))*dt;
        
        POP(:,7)=POP(:,7) + (-r13.*POP(:,7) - r11.*POP(:,7) - r16*POP(:,7) + r10.*POP(:,6))*dt;
        POP(:,6)=POP(:,6) + (-r14.*POP(:,6) - r10.*POP(:,6) - r12.*POP(:,6) + r9.*POP(:,5))*dt;
        POP(:,5)=POP(:,5) + (-r9.*POP(:,5) + r8'.*POP(:,1) - r18.*POP(:,5))*dt;      
        
        POP(:,1)=POP(:,1) + (-r1'.*POP(:,1) - r8'.*POP(:,1) + r19.*POP(:,9))*dt;

    end
end

S=sum(POP(21,:));
disp(S);
disp(POP(:,1));