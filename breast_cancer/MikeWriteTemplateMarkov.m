function [AA,BB,CC,DD,EE,FF,GG,HH] = MikeWriteTemplateMarkov(TEMPMARKOV1,TEMPMARKOV2,TEMPMARKOV3,TEMPMARKOV4,TEMPMARKOV5,TEMPMARKOV6,TEMPMARKOV7,TEMPMARKOV8,CountryIndex,POP,OnsetRatesF,CFR,diagr,DwellTimesF,DiagsA,StageDist,AgeArrayMid,PrevA,prev,CountryName)

            AA = TEMPMARKOV1;
            BB = TEMPMARKOV2;
            CC = TEMPMARKOV3;
            DD = TEMPMARKOV4;
            EE = TEMPMARKOV5;
            FF = TEMPMARKOV6;
            GG = TEMPMARKOV7;
            HH = TEMPMARKOV8;
            
%             POP
%             OnsetRatesF
%             CountryIndex
            
CountryIndex = CountryIndex - 1;

sAA = CountryIndex .* 5;
AA{4,5+sAA} = CountryName;

for icol=5:8
    if icol ==5
        for irow=7:107
            AA{irow,icol+sAA}=POP(irow-6,1);
        end
    elseif icol==6
        for irow=7:107
            AA{irow,icol+sAA}=sum(POP(irow-6,(2:5)),2);
        end
    elseif icol==7
        for irow=7:107
            AA{irow,icol+sAA}=sum(POP(irow-6,(6:9)),2);
        end
    elseif icol==8
        for irow=7:107
            AA{irow,icol+sAA}=OnsetRatesF(irow-6,1);
        end
    end
end

sCC = CountryIndex .* 2;
CC{4,5+sCC} = CountryName;


icol=5;
   for irow=7:107
      CC{irow,icol+sCC}=OnsetRatesF(irow-6,1);
   end
   
 tmax=500;
% t=zeros(1,tmax);
% t(1,1)=1;
% 
% for T=2:tmax
%     t(1,T)=t(1,T-1)+ 1;
% end

sDD = CountryIndex .* 5;
DD{4,5+sDD} = CountryName;

   for icol=5:8
       for irow=7:107
           DD{irow,icol+sDD}=CFR(irow-6,tmax-(100*(icol-5)));
       end
   end
      
   sEE = CountryIndex .* 5;
   EE{4,5+sEE} = CountryName;

   for icol=5:8
       for irow=7:107
           EE{irow,icol+sEE}=diagr(irow-6,icol-4);
       end
   end
   
   sFF = CountryIndex .* 5;
   FF{4,5+sFF} = CountryName;

   for icol=5:8
       for irow=7:107
           FF{irow,icol+sFF}=DwellTimesF(icol-4,irow-6);
       end
   end
   
   sGG = CountryIndex .* 3;
   GG{4,5+sGG} = CountryName;

            tmax=500;
   
icol=5;
       for irow=7:10
           for t = 1:tmax
               
               
%            tmax=500;
%            t=zeros(1,tmax);
%            t(1,1)=1;
%            for T=2:tmax
%                t(1,T)=t(1,T-1)+ 1;
%            end

           GG{irow,icol+sGG}=DiagsA(irow-6,t);
           end
       end
              
       icol=6;
       for irow=7:10
           GG{irow,icol+sGG}=StageDist(irow-6);
       end
           
       sHH = CountryIndex .* 4;
       HH{4,5+sHH} = CountryName;

   for icol=5:7
       for irow= 7:(tmax+6)
           HH{irow,icol+sHH}=CFR(icol-3,irow-6);
       end
   end
       
%   AgeArrayMid
   
       sBB = CountryIndex .* 3;
       BB{4,5+sBB} = CountryName;

icol=3;
for irow=7:18
    BB{irow,icol}=AgeArrayMid(1,irow-6);
end
   
   icol=5;
for irow=7:18
    BB{irow,icol+sBB}=(prev(1,irow-6)*100);
end
   
icol=6;
for irow=7:18
    BB{irow,icol+sBB}=PrevA(1,irow-6);
end

end

