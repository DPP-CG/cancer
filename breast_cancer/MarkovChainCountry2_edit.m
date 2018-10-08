function [RES]=MarkovChainCountry2(PAR)

% [~,~,TEMPMARKOV1]=xlsread('C:\Users\Mike\Documents\Umass\Spring 2015\Cancer Research\excel\Template for OutputData Markov.xls','Total Pop');
% [~,~,TEMPMARKOV2]=xlsread('C:\Users\Mike\Documents\Umass\Spring 2015\Cancer Research\excel\Template for OutputData Markov.xls','prev & prevA');
% [~,~,TEMPMARKOV3]=xlsread('C:\Users\Mike\Documents\Umass\Spring 2015\Cancer Research\excel\Template for OutputData Markov.xls','Onset Rate');
% [~,~,TEMPMARKOV4]=xlsread('C:\Users\Mike\Documents\Umass\Spring 2015\Cancer Research\excel\Template for OutputData Markov.xls','CFR');
% [~,~,TEMPMARKOV5]=xlsread('C:\Users\Mike\Documents\Umass\Spring 2015\Cancer Research\excel\Template for OutputData Markov.xls','Time to diagnosis');
% [~,~,TEMPMARKOV6]=xlsread('C:\Users\Mike\Documents\Umass\Spring 2015\Cancer Research\excel\Template for OutputData Markov.xls','Dwell Time');
% [~,~,TEMPMARKOV7]=xlsread('C:\Users\Mike\Documents\Umass\Spring 2015\Cancer Research\excel\Template for OutputData Markov.xls','Stage at diagnosis');
% [~,~,TEMPMARKOV8]=xlsread('C:\Users\Mike\Documents\Umass\Spring 2015\Cancer Research\excel\Template for OutputData Markov.xls','DiagsA');

RES=[];

close all
nsims=10000;
%Country-specific parameters
probAgeA=PAR{1};
MortalityArray=PAR{2};
for i = 80:101
MortalityArray(i)=MortalityArray(i-1)* MortalityArray(79)/MortalityArray(78);
end
Births=PAR{3};
prevalenceActualF=PAR{4};
%prevN=prevalenceActual/sum(prevalenceActual);

%Prevalence and case fatality rate: Age Group is different from Spectrum age group so converting
AgeArrayL = PAR{7};
AgeArrayU = PAR{8};
AgeArrayM=(AgeArrayU+AgeArrayL)/2;


AgeArrayLM =[1    5   15  30	45	60	70	80];
AgeArrayUM =[4    14  29 	44	59	69	79	100];
AgeArrayMM=(AgeArrayUM+AgeArrayLM)/2;

CFR = interp1(AgeArrayMM,PAR{5},[0:100],'cubic');
%CFR(1,1:14)=0;
cancerMortalityArray =-log(1-min(.99,CFR));
%cancerMortalityArray(1,80:end) = cancerMortalityArray(1,79);
%cancerMortalityArray = max(cancerMortalityArray, MortalityArray);
PAR{5} = cancerMortalityArray;
StageDist=PAR{6};


% CountryIndex=PAR{9};


%Interpolation to single age years, flatline after 80 years
PrevalenceF = interp1(AgeArrayM,prevalenceActualF,[0:100],'cubic');
%PrevalenceF(80:end)=PrevalenceF(80);
PrevalenceF= max(0,PrevalenceF);



%Spectrum Age groups
AgeArrayLower=[1    5	10	15  20	25	30	40	50	60	70	80];
AgeArrayUpper=[4	9	14  19 	24	29	39	49	59	69	79	100];
AgeArrayMid=((AgeArrayUpper+AgeArrayLower)/2);

prevalenceActual = zeros(size(prevalenceActualF,1),size(AgeArrayLower,2));

for i = 1 : size(prevalenceActual,2)
     indA=AgeArrayLower(i):AgeArrayUpper(i);
    prevalenceActual(1,i)= dot(PrevalenceF(indA),probAgeA(indA))/sum(probAgeA(indA)) ;%PrevalenceF(AgeArrayMid(i));%
end
%prevalenceActual(1,1:3) =0;
PAR{4}=prevalenceActual;

%Progression times in undiagnoses cases
DwellTimes = zeros(size(AgeArrayLower,1),size(AgeArrayLower,2));
for i = 1 : size(AgeArrayLower,2) 
   %    DwellTimes(1,i)= 0.6841*exp( 0.0306*AgeArrayMid(i));
%    DwellTimes(2,i)= 1.2086*exp( 0.0306*AgeArrayMid(i));
%    DwellTimes(3,i)= 2.5997*exp( 0.0306*AgeArrayMid(i));
%    DwellTimes(4,i)= .3*exp(0.0306*AgeArrayMid(i));

% Below commented by Prashant on 10 July 2016. 
%   DwellTimes(1,i)= 1/(1.661*.97^AgeArrayMid(i));
%   DwellTimes(2,i)= 1/(1.661*(3/5.3)*.97^AgeArrayMid(i));
%   DwellTimes(3,i)= 1/(1.661*(3/8.4)*.97^AgeArrayMid(i));
%   DwellTimes(4,i)= 1/(1.661*(3/2.5)*.97^AgeArrayMid(i));
   
%Using below values added by Prashant on 10 July 2016 
    DwellTimes(1,i)= 1/(1.661*(3/5.22)*.97^AgeArrayMid(i));
    DwellTimes(2,i)= 1/(1.661*(3/3.045 )*.97^AgeArrayMid(i));
    DwellTimes(3,i)= 1/(1.661*(3/2.32)*.97^AgeArrayMid(i));
    DwellTimes(4,i)= 1/(1.661*(3/2)*.97^AgeArrayMid(i));

%   DwellTimes(1,i)= 1.6 + 0.2* (AgeArrayMid(i)*AgeArrayMid(i)) ./ ( 3092.3+ (-108.3)*AgeArrayMid(i)+ AgeArrayMid(i)*AgeArrayMid(i));
%   DwellTimes(2,i)= 2.8 + 0.3* (AgeArrayMid(i)*AgeArrayMid(i)) ./ ( 3092.3+ (-108.3)*AgeArrayMid(i)+ AgeArrayMid(i)*AgeArrayMid(i));
%   DwellTimes(3,i)= 6.1 + 0.7* (AgeArrayMid(i)*AgeArrayMid(i)) ./ ( 3092.3+ (-108.3)*AgeArrayMid(i)+ AgeArrayMid(i)*AgeArrayMid(i));
%   DwellTimes(4,i)= .27 + 0.03* (AgeArrayMid(i)*AgeArrayMid(i)) ./ ( 3092.3+ (-108.3)*AgeArrayMid(i)+ AgeArrayMid(i)*AgeArrayMid(i));
end

% for i = size(AgeArrayLower,2): size(AgeArrayLower,2)
%   DwellTimes(1,i)= DwellTimes(1,i-1);
%   DwellTimes(2,i)= DwellTimes(2,i-1);
%   DwellTimes(3,i)= DwellTimes(3,i-1);
%   DwellTimes(4,i)= DwellTimes(4,i-1);
%    
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Interpolation to single age years, flatline after 80 years

DwellTimes1 = interp1(AgeArrayMid,DwellTimes(1,:),[0:100],'cubic');
% DwellTimes1(80:end)=DwellTimes1(80);

DwellTimes2 = interp1(AgeArrayMid,DwellTimes(2,:),[0:100],'cubic');
% DwellTimes2(80:end)=DwellTimes2(80);

DwellTimes3 = interp1(AgeArrayMid,DwellTimes(3,:),[0:100],'cubic');
% DwellTimes3(80:end)=DwellTimes3(80);

DwellTimes4 = interp1(AgeArrayMid,DwellTimes(4,:),[0:100],'cubic');
% DwellTimes4(80:end)=DwellTimes4(80);

DwellTimesF=[DwellTimes1;DwellTimes2;DwellTimes3;DwellTimes4];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculates 10-year relative survival with and without treatment, which are
%used in Excel file to calculate  relative mortality avergaed over the 3
%countries. 
%Also calculates coverage and stage distributions, which would be good for
%verfication purposes.
[Mu_Tx, ImpactAge, BaseMortality, TxCoverage, StgDist, RelativeMortTx, RelativeMortNoTx] = Mortality(DwellTimesF, MortalityArray, PrevalenceF, StageDist, smooth(PAR{5},'lowess'),PAR);

% Calculated in Excel using 10-year relative survival  from function
% Mortality.
% RelativeMortTx = [[3.97882963	5.335307383	10.54064409	17.95052658	24.55971376	30.62048696	35.88385674	39.99670009	43.6168351	47.27567044	50.96389799	53.49259238	53.38904377	50.49120674	45.59469549	40.04857206	34.82152322	30.6945014	27.42682732	24.36682523	21.34173549	18.70471483	16.31170271	14.21898987	12.40054582	10.86314741	9.470841061	8.202383723	7.032369877	5.983429209	5.029044911	4.171832886	3.404892669	2.713401577	2.08123728	1.502830822	0.985162899	0.518299796	0.102981294	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
% ];
% [4.273832062	5.760410716	11.4442383	19.5937031	26.95484909	33.7992164	39.84731911	44.69585968	49.06728891	53.55937193	58.16989745	61.54023407	61.93782091	59.09993124	53.87715688	47.8049938	42.01812082	37.47102391	33.90342943	30.53058202	27.1345778	24.16341722	21.44166986	19.0512194	16.96940859	15.21946733	13.62408353	12.1580456	10.78717681	9.549925674	8.410837599	7.379999024	6.453530524	5.612702015	4.831136689	4.097679539	3.434085253	2.817070016	2.253701037	1.74539755	1.288945955	0.875914948	0.507041783	0.176297891	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
% ];
% [5.519525921	7.524634991	15.12911754	26.21620679	36.51643077	46.38470747	55.42713456	63.05219736	70.24353	77.86104059	85.93301791	92.45350929	94.70420939	92.04935793	85.55694502	77.47533143	69.56978274	63.45496199	58.79453168	54.2919158	49.55170487	45.38483355	41.49266029	38.05531071	35.06322093	32.60605047	30.34333627	28.23318966	26.20455937	24.3595052	22.62517856	21.04232702	19.62041642	18.32464601	17.0835151	15.8589823	14.73788301	13.6347595	12.58410012	11.60808123	10.71328395	9.860387669	9.088748894	8.379684993	7.723335506	7.102437946	6.520244729	5.912874094	5.304322855	4.724970495	4.201538026	3.71171759	3.298648023	2.945662601	2.619122863	2.303092492	2.018643718	1.751703873	1.5015148	1.275772666	1.078118612	0.901975541	0.748659533	0.616803305	0.504647072	0.410280924	0.332179603	0.269203462	0.219507352	0.181023666	0.151994722	0.130004106	0.113462997	0.10120315	0.091805193	0.084129921	0.076971885	0.070277769	0.063803768	0.057564613	0.057564613	0.057564613	0.057564613	0.057564613	0.057564613	0.057564613	0.057564613	0.057564613	0.057564613	0.057564613	0.057564613	0.057564613	0.057564613	0.057564613	0.057564613	0.057564613	0.057564613	0.057564613	0.057564613	0.057564613	0.057564613
% ];
% [55.96546778	76.02707228	152.3513649	263.1794219	365.5142194	463.0248348	551.8774169	626.3056392	696.1958921	770.1151309	848.3509932	911.1445034	931.8556088	904.4436996	839.5811565	759.4205819	681.264328	620.8703202	574.8791518	530.5706157	484.0602659	443.2503898	405.2028602	371.6595822	342.5134632	318.6310416	296.6792791	276.2419153	256.6172041	238.798781	222.0691617	206.8250941	193.1588689	180.7283123	168.8273799	157.0763103	146.3318072	135.7447232	125.6547526	116.2824994	107.6953049	99.49917924	92.09179748	85.28755988	78.98867403	73.02006289	67.41554527	61.51159553	55.54793599	49.83710804	44.6616779	39.7883485	35.68455599	32.18337686	28.92587254	25.73625023	22.84970982	20.1132479	17.52078999	15.1636106	13.08911227	11.22211402	9.581568102	8.15465105	6.923527745	5.868089794	4.973733038	4.228344176	3.614376804	3.113125128	2.711732845	2.385840216	2.122396742	1.912750531	1.740801664	1.593605833	1.456000469	1.327915662	1.204613278	1.086279894	1.086279894	1.086279894	1.086279894	1.086279894	1.086279894	1.086279894	1.086279894	1.086279894	1.086279894	1.086279894	1.086279894	1.086279894	1.086279894	1.086279894	1.086279894	1.086279894	1.086279894	1.086279894	1.086279894	1.086279894	1.086279894
% ]];
% 
% % Calculated in Excel using 10-year relative survival  from function
% % Mortality.
% RelativeMortNoTx = [[4.579092252	6.182872873	12.30911831	21.1242926	29.13876153	36.6457687	43.34270646	48.78758706	53.76450772	58.93194484	64.29709251	68.36189449	69.179435	66.40494775	60.93413101	54.45637587	48.24358009	43.39797572	39.64366395	36.07888233	32.44217361	29.26516419	26.34261288	23.77987006	21.55811702	19.71875157	18.04374723	16.50332573	15.05276467	13.74750045	12.54156207	11.45378259	10.48469914	9.611307869	8.793114097	8.009302934	7.302350178	6.627525155	6.000423239	5.429616312	4.915751365	4.438633349	4.014187938	3.631907987	3.285700238	2.96674465	2.675192598	2.384096466	2.103045144	1.843396368	1.614319003	1.405854795	1.233021556	1.088067096	0.95749116	0.834778032	0.726908133	0.628148158	0.537693045	0.457818461	0.389399973	0.329682851	0.278755596	0.235778154	0.199792376	0.169809241	0.145044997	0.124859052	0.108504657	0.095284224	0.084764839	0.076127052	0.068973931	0.063139695	0.058214245	0.053869107	0.049682708	0.04567334	0.041691937	0.037752323	0.037752323	0.037752323	0.037752323	0.037752323	0.037752323	0.037752323	0.037752323	0.037752323	0.037752323	0.037752323	0.037752323	0.037752323	0.037752323	0.037752323	0.037752323	0.037752323	0.037752323	0.037752323	0.037752323	0.037752323	0.037752323
% ];
% [4.946164883	6.709812413	13.42351839	23.14575541	32.07938102	40.54188324	48.19384561	54.53256737	60.42181337	66.60153123	73.08800331	78.17677721	79.60478571	76.90486626	71.03956835	63.9250974	57.03534506	51.68422547	47.57209325	43.63443182	39.5541083	35.97877552	32.66457512	29.74825992	27.21520298	25.12737402	23.21566993	21.44532132	19.76054573	18.23650785	16.81598746	15.52719497	14.37458344	13.33025679	12.340558	11.37727612	10.50189086	9.652178005	8.851832652	8.115321774	7.445920925	6.815124552	6.249139224	5.733942131	5.261821858	4.820215678	4.410645689	3.989342552	3.572049129	3.178563458	2.826097273	2.498919024	2.22543986	1.994038567	1.78163172	1.576911609	1.393831473	1.222511619	1.062171173	0.917987579	0.792487166	0.680798034	0.583740753	0.500221979	0.428848723	0.368096944	0.316814023	0.273984467	0.238389864	0.208872081	0.184771746	0.164653165	0.147879912	0.134144375	0.122574963	0.112473253	0.102956655	0.094050322	0.085425461	0.077099031	0.077099031	0.077099031	0.077099031	0.077099031	0.077099031	0.077099031	0.077099031	0.077099031	0.077099031	0.077099031	0.077099031	0.077099031	0.077099031	0.077099031	0.077099031	0.077099031	0.077099031	0.077099031	0.077099031	0.077099031	0.077099031
% ];
% [5.91107613	8.077642699	16.28194486	28.28520287	39.50147355	50.3124205	60.28921628	68.78242581	76.85799629	85.45866585	94.62350247	102.1451015	104.9961206	102.4216975	95.55502827	86.86625204	78.31824993	71.73506831	66.75743785	61.92581168	56.78697217	52.26815796	48.03086733	44.28718784	41.03222764	38.37821807	35.93122822	33.64380049	31.4322988	29.42017674	27.5219053	25.78878619	24.23514086	22.82082824	21.45845362	20.10003942	18.85556592	17.61688734	16.42817958	15.31922053	14.30037061	13.32061119	12.43413802	11.61760936	10.85899021	10.13514138	9.451360485	8.714337559	7.956093392	7.220509332	6.549078398	5.908788795	5.370440581	4.912057874	4.48088468	4.049886097	3.656131656	3.276494443	2.910575051	2.573911957	2.275445164	2.00307098	1.760688167	1.546791249	1.35891128	1.194056782	1.050321106	0.925741484	0.818048618	0.725143412	0.646431206	0.578529315	0.520426204	0.471840068	0.430279107	0.393730962	0.359494166	0.327689398	0.297138574	0.267880777	0.267880777	0.267880777	0.267880777	0.267880777	0.267880777	0.267880777	0.267880777	0.267880777	0.267880777	0.267880777	0.267880777	0.267880777	0.267880777	0.267880777	0.267880777	0.267880777	0.267880777	0.267880777	0.267880777	0.267880777	0.267880777
% ];
% [60.11975115	81.89466212	164.5843698	285.1344904	397.1903664	504.7047851	603.4732144	687.1146361	766.388855	850.7417704	940.5754434	1013.993129	1041.075185	1014.517103	945.683256	859.0796678	774.1057332	708.741482	659.3841783	611.5841233	560.8434284	516.2986821	474.5887235	437.7946609	405.8589813	379.8877363	355.9804479	333.6617727	312.0964021	292.5050059	274.0355451	257.1967934	242.1325712	228.4440085	215.2563717	202.0844478	190.0305749	178.0048315	166.4497318	155.6665646	145.7627967	136.2203035	127.5942286	119.6495554	112.2654338	105.2043858	98.52187863	91.2420866	83.68987893	76.321073	69.57503012	63.10490512	57.67169039	53.0521155	48.68436546	44.27484105	40.2284834	36.29618888	32.47568638	28.94147349	25.79719977	22.90904444	20.32338796	18.02594198	15.99134043	14.1879987	12.59726355	11.19827169	9.968975004	8.890330607	7.961807347	7.148721035	6.444151727	5.848898158	5.335497909	4.881749126	4.456611958	4.061873095	3.682876777	3.320072183	3.320072183	3.320072183	3.320072183	3.320072183	3.320072183	3.320072183	3.320072183	3.320072183	3.320072183	3.320072183	3.320072183	3.320072183	3.320072183	3.320072183	3.320072183	3.320072183	3.320072183	3.320072183	3.320072183	3.320072183	3.320072183
% ]];
% 

%average values into age groups
MortA=[];
MortC=[];
for J=1:length(AgeArrayLower)
    indA=AgeArrayLower(J):AgeArrayUpper(J);
    MortA(J)=dot(MortalityArray(indA),probAgeA(indA))/sum(probAgeA(indA));
    MortC(J)=dot(cancerMortalityArray(indA),probAgeA(indA))/sum(probAgeA(indA));
    for i = 1:4 % if the denominator is equal to 0, all the values are 0.
        if sum(PrevalenceF(indA))==0
            BaseMortalityC(i,J)=0;
            Impact(i,J)=0;
            TxCoverage_A(i,J)=0;
            StgDistribution(i,J)=0;
            DwellTimes(i,J)=0;
            %%%%%%%
            MORT(i,J)=0;
            BASE(i,J)=0;
        else
        BaseMortalityC(i,J) = dot(RelativeMortNoTx(i,indA),PrevalenceF(indA))/sum(PrevalenceF(indA));
        Impact(i,J) =         dot(RelativeMortTx(i,indA),PrevalenceF(indA))/sum(PrevalenceF(indA));
        TxCoverage_A(i,J) = dot(TxCoverage(indA,i)',PrevalenceF(indA))/sum(PrevalenceF(indA));
        StgDistribution(i,J) = dot(StgDist(indA,i)',PrevalenceF(indA))/sum(PrevalenceF(indA));
        DwellTimes(i,J) = dot(DwellTimesF(i, indA),PrevalenceF(indA))/sum(PrevalenceF(indA)); 
        %%%%%%%
        MORT(i,J)= dot((MortalityArray(indA)+MortalityArray(indA).*RelativeMortTx(i,indA)),PrevalenceF(indA))/sum(PrevalenceF(indA));
        BASE(i,J) = dot((MortalityArray(indA)+MortalityArray(indA).*RelativeMortNoTx(i,indA)),PrevalenceF(indA))/sum(PrevalenceF(indA));
        end
    end
  end
MortA(end)=MortA(end-1);                                                                                                              
MortC(end)=MortC(end-1);




%TimeToDiagnosis
TimeToDiagnosis = zeros(size(DwellTimes,1),size(DwellTimes,2));
for i = 1 : size(DwellTimes,2)
   
    sumTime = zeros(size(DwellTimes,1));
    for j = 1:nsims
       randNum = rand;
       Insitu = - log(randNum) * DwellTimes(1,i);
       Local = - log(randNum) * DwellTimes(2,i);
       Regional = - log(randNum) * DwellTimes(3,i);
       Distant = - log(randNum) * DwellTimes(4,i); 
       sumTime(1) = sumTime(1) + 0.5*Insitu;
       sumTime(2) = sumTime(2) + Insitu + 0.5*Local;
       sumTime(3) = sumTime(3) + Insitu + Local + 0.5*Regional;
       sumTime(4) = sumTime(4) + Insitu + Local + Regional + 0.5*Distant; 
    end
   TimeToDiagnosis(1,i) = sumTime(1) / nsims;
   TimeToDiagnosis(2,i) = sumTime(2) / nsims;
   TimeToDiagnosis(3,i) = sumTime(3) / nsims;
   TimeToDiagnosis(4,i) = sumTime(4) / nsims;
end

TimeToDiagnosis

%Modified TimeToDiagnosis (using DCO to estimate event of first occurence )
ModTimeToDiagnosis = zeros(size(DwellTimes,1),size(DwellTimes,2));

CumStageDISTwithDCO = StageDist;
DCO = 0.04;% proportion DCO
for i = 1: size (StageDist,2)
   % CumStageDISTwithDCO(i) = StageDist(i) / (sum(StageDist(i:end))+ DCO);
    CumStageDISTwithDCO(i) = sum(StageDist(1:i)) / (sum(StageDist(1:end))+ DCO);
end 

%CumStageDISTwithDCO = CumStageDISTwithDCO / (1+DCO);

for i = 1 : size(DwellTimes,2)
   ModTimeToDiagnosis(1,i) = (1 - CumStageDISTwithDCO (1))/((1/DwellTimes(1,i) + MortA(1,i)) * CumStageDISTwithDCO(1));
   ModTimeToDiagnosis(2,i) = (1 - CumStageDISTwithDCO (2))/((1/(DwellTimes(1,i) + DwellTimes(2,i)) + MortA(1,i))* CumStageDISTwithDCO(2));
   ModTimeToDiagnosis(3,i) = (1 - CumStageDISTwithDCO (3))/((1/(DwellTimes(1,i) + DwellTimes(2,i) + DwellTimes(3,i))+ MortA(1,i))* CumStageDISTwithDCO(3));
   ModTimeToDiagnosis(4,i) = (1 - CumStageDISTwithDCO (4))/((1/(DwellTimes(1,i) + DwellTimes(2,i) + DwellTimes(3,i) + DwellTimes(4,i))+ MortA(1,i))* CumStageDISTwithDCO(4));
end 


% for i = 1 : size(DwellTimes,2)
%    ModTimeToDiagnosis(1,i) = DwellTimes(1,i) * (1 - CumStageDISTwithDCO (1))/CumStageDISTwithDCO(1);
%    ModTimeToDiagnosis(2,i) = (DwellTimes(1,i) + DwellTimes(2,i)) * (1 - CumStageDISTwithDCO (2))/CumStageDISTwithDCO(2);
%    ModTimeToDiagnosis(3,i) = (DwellTimes(1,i) + DwellTimes(2,i) + DwellTimes(3,i)) * (1 - CumStageDISTwithDCO (3))/CumStageDISTwithDCO(3);
%    ModTimeToDiagnosis(4,i) = (DwellTimes(1,i) + DwellTimes(2,i) + DwellTimes(3,i) + DwellTimes(4,i)) * (1 - CumStageDISTwithDCO (4))/CumStageDISTwithDCO(4);
% end 
 %SojournTime and SojournTimePlusInsitu
 SojournTime =  zeros(1,size(DwellTimes,2));
 SojournTimePlusInsitu = zeros(1,size(DwellTimes,2));
 for j = 1: size(DwellTimes,2)
    SojournTime(1,j) = StageDist(1) * TimeToDiagnosis(1,j) + StageDist(2) * TimeToDiagnosis(2,j) + StageDist(3) * TimeToDiagnosis(3,j) +StageDist(4) * TimeToDiagnosis(4,j);
    SojournTimePlusInsitu(1,j) = SojournTime(1,j);
    
%     sumTime = 0;
%     for z = 1: nsims
%         randNum = rand;
%         sumTime = sumTime - log(randNum) * SojournTime(1,j) - log(randNum) * DwellTimes(1,j);
%     end
%     SojournTimePlusInsitu(1,j) = sumTime / nsims;

end%for j

sumprobAgeA = zeros(1, size(AgeArrayLower,2));
for j = 1: size(AgeArrayLower,2)
    for  probAgeInd = AgeArrayLower(j) : AgeArrayUpper(j)
        sumprobAgeA(1,j) = sumprobAgeA(1,j) + probAgeA(probAgeInd);
    end
end
    
%Calculating incidence by solving the differential equations 
Incidence = IncidenceFunc(probAgeA,MortalityArray,prevalenceActual,cancerMortalityArray, AgeArrayLower,AgeArrayUpper);   
%Incidence = Incidence *1000;
Incidence =PAR{4};%GLOBOCAN breast cancer data has incidence no prevalence
%Incidence = transpose(smooth(Incidence));

%Calculating incidence by solving the differential equations MArkov Chain of the Markov Process 
OnsetRates = PrevStageRates (Incidence,TimeToDiagnosis,StageDist /(1+DCO),probAgeA,MortalityArray,prevalenceActual,cancerMortalityArray, AgeArrayLower,AgeArrayUpper);
% OnsetRates = smooth(OnsetRates);
% OnsetRates = OnsetRates';
% OnsetRates(1,1:4) =0;

OnsetMatrix = zeros(1, size(AgeArrayLower,2)+ 1); 

OnsetMatrix(1,1:size(AgeArrayLower,2))= OnsetRates * 1000;
OnsetMatrix(1,size(AgeArrayLower,2)+1)= (OnsetRates * 1000) * transpose(sumprobAgeA);
 %sum((prevalenceActual) * (transpose((SojournTimeMat(i,:)) * transpose(sumprobAgeA)))) / sum(prevalenceActual)
 prevalenceTotal = prevalenceActual(1: size(AgeArrayLower,2))* transpose(sumprobAgeA) *1000;
%OnsetRates(1,1:6)= 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Interpolation to single age years, flatline after 80 years

IncidenceF = interp1(AgeArrayMid,Incidence,[0:100],'cubic');
%IncidenceF(80:end)=IncidenceF(80);

OnsetRatesF = interp1(AgeArrayMid,OnsetRates,[0:100],'cubic');
%OnsetRatesF(80:end)=OnsetRatesF(80);
%OnsetRatesF(1:20) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i = 1 : size(OnsetRates,2)
%      indA=AgeArrayLower(i):AgeArrayUpper(i);
%     onsetRateSmooth(1,i)= dot(OnsetRatesF(indA),PrevalenceF(indA))/sum(PrevalenceF(indA)) ;%PrevalenceF(AgeArrayMid(i));%
% end


% 
% beta0=[5,0.001, 5,5];
    
%% calculation of diagnostic rates starts here 
    
    for upper_start = 2:2
        % here trials for different values of starting points starts
        % these trials are essentially for finding the diagnostic rate values
        % which fits the incidence by Age and Stage
        trial_run = 1;

        % if we are fitting the diagnostic rates value by age and stage then we
        % need initialization for all the values of diagnostic rates, hence the
        % beta0 matrix will be 4x12. When we are fitting diagnostic rates values by
        % age only, in that case beta0 matrix will be 1x12
        betaRuns = ones(4,size(AgeArrayLower,2),trial_run);

        for randRuns = 1:trial_run

            % randomly choosing starting point
            beta0 = ones(4,size(AgeArrayLower,2)) * rand() * upper_start;
            % storing initial point for analyzing results
            initial_point=beta0(:,1);

            % following commented fuction (ExpSourj) is used only when we are
            % finding diagnostic rates which fit incidence/prevalence by Age
            % ONLY. As we now are fitting to Age and Stage, we don't use this
            % function (it will not be needed as curvefitting function will
            % give us each and every value. Hence, no need to generate values
            % for PC-1,2,3 from PC-4)

            %[diagr]=ExpSourj(beta0, StageDist,AgeArrayLower, AgeArrayUpper);

            %test shape by displaying
            % figure
            % hold on
            % X=0:100;
            % plot(X,diagr(:,1),'b')
            % plot(X,diagr(:,2),'r')
            % plot(X,diagr(:,3),'g')
            % plot(X,diagr(:,4),'c')
            % pause

            % following X and Y matrices are input for the curve fitting
            % function (Y can be incidence or Prevalence)
            X=AgeArrayMid;
            %Y=[prevalenceActual*100];%if fit to prevalence
            Y=StageDist(1,1:4)'*[Incidence*100];%if fit to incidence
            %[Z,diagr]=PopSim(beta0,X,PAR,IncidenceF,OnsetRatesF,DwellTimesF,AgeArrayLower,AgeArrayUpper,1)

            %Fit the prevalence model, find the diagnostic rates and plot the results
            tic
            [Z,PrevS,PrevC,diagrT,beta]=FitSoJourn(1000,beta0,X,Y,PAR,Incidence,OnsetRatesF,DwellTimesF,AgeArrayLower,AgeArrayUpper, Mu_Tx, BaseMortality);
            toc
            betaRuns(:,:,randRuns) = beta;
            initial_point_array(:,randRuns) = initial_point;
        end

        trial_detail=zeros(4,size(beta,2)+1,trial_run);
        trial_detail(:,1,:) = initial_point_array;
        trial_detail(:,2:(size(beta,2)+1),:) = betaRuns;

        for i = 1:trial_run
            trial_details(:,:,i) = transpose(trial_detail(:,:,i));
        end

        %[DiagnosticRate]=ExpSourjX(beta,AgeArrayMid, StageDist);

        rand_i = randi(trial_run);
        if trial_details(1,1,rand_i) <= 3.5
            DiagnosticRate = betaRuns(:,:,rand_i)';
        else
            while trial_details(1,1,rand_i) > 3.5
                rand_i = randi(trial_run);
                DiagnosticRate = betaRuns(:,:,rand_i)';
            end
        end

        % collecting output in tractable format
        for i_l = 1:trial_run
            for j = 1:4
                stack((j-1)*13 + 1:(j-1)*13 + 13,i_l) = trial_details(:,j,i_l);
            end
        end

        % storing the results
       RES{12 + upper_start} = stack;
   
    end
    
%% diagnostic rates calculation ends here

%save results to be mapped to associaitons file
RES{1}=PrevS;
RES{2}=PrevC;
RES{3}=Incidence;
RES{4}=OnsetRates';%onsetRateSmooth';
RES{5}=DiagnosticRate;
RES{6}=DwellTimes';
RES{7}=BaseMortalityC; %Relative Mortality without treatmet
RES{8}=MortA';% Disease free mortality
RES{9} = Impact;% %Relative Mortality with treatment
RES{10} = TxCoverage_A;% Current coverage of treatment

%pause

% [AA,BB,CC,DD,EE,FF,GG,HH] = MikeWriteTemplateMarkov(TEMPMARKOV1,TEMPMARKOV2,TEMPMARKOV3,TEMPMARKOV4,TEMPMARKOV5,TEMPMARKOV6,TEMPMARKOV7,TEMPMARKOV8,CountryIndex)
% 
%   xlswrite('Mike Markov Results.xls',AA,'Total Pop')
%   xlswrite('Mike Markov Results.xls',BB,'prev & prevA')
%   xlswrite('Mike Markov Results.xls',CC,'Onset Rate')
%   xlswrite('Mike Markov Results.xls',DD,'CFR')
%    xlswrite('Mike Markov Results.xls',EE,'Time to diagnosis')
%    xlswrite('Mike Markov Results.xls',FF,'Dwell Time')
%    xlswrite('Mike Markov Results.xls',GG,'Stage at diagnosis')
%    xlswrite('Mike Markov Results.xls',HH,'DiagsA')
%   
% 
%   'Mike Markov Results Written'
% 
% 
% pause

end%end of MarkvoChainCountry

function [PrevA,PrevS,PrevC,diagr,beta]=FitSoJourn(MaxFunEvals,beta0,X,Y,PAR,Incidence,OnsetRatesF,DwellTimesF,AgeArrayLower,AgeArrayUpper, Mu_Tx, BaseMortality)

% defining lower and upper bound
% LB and UB matrices will be 4x12 (Age-stage fitting)
LB = ones(4,size(AgeArrayLower,2))*0;
UB = ones(4,size(AgeArrayLower,2))*4;

% Strict condition on initial age grops
LB(:,1:3) = 0;
UB(:,1:3) = 0;
% strict condition on PC-1 (all age groups), as there are no cases to
% diagnose in these states, we want the rates to be strictly zero
StageDist = PAR{6};
if StageDist(1,1) == 0
    LB(1,:) = 0;
    UB(1,:) = 0;
end

%LB=[-10,-1, -10,-10];
%UB=[100,0.5, 10, 100];%if exponential curve for non-US
%UB = [2000, 0.5,10,100];%% if logistic US
%UB = [10000, 0.5,10,1000];%% if logistic non-US

%OPTIONS = optimset('largescale','off','LevenbergMarquardt','on','display','off','Diagnostics','on','MaxFunEvals',700,'TolFun',1e-8,'TolCon',1e-8);
OPTIONS = optimset('MaxFunEvals',MaxFunEvals,'Diagnostics','on','TolFun',1e-6,'TolCon',1e-6);

%% Curve fitting: fitting values by Age AND Stage
store_beta=beta0;
for age = 1:size(AgeArrayLower,2)
    % selecting current X value 
    individualX = X(1,age);
    
    for state_c = 1 : 4 % PC-1 to PC-4
        
        tic
        
        [betafit_out,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@(x,individualX )...
        PopSim(x,individualX,PAR,Incidence,OnsetRatesF,DwellTimesF,AgeArrayLower,AgeArrayUpper,Mu_Tx, BaseMortality, 0, age,store_beta,Y,state_c),...
        beta0(1,age),individualX,Y(state_c,age),LB(state_c,age),UB(state_c,age),OPTIONS);

        store_beta(state_c,age) = betafit_out;
        
        toc
    end
    
end
beta=store_beta;

%plot the results
[PrevA,diagr,PrevS,PrevC]=PopSim(beta,X,PAR,Incidence,OnsetRatesF,DwellTimesF,AgeArrayLower,AgeArrayUpper,Mu_Tx, BaseMortality,1,10^5,store_beta,Y,state_c);

end

function [diagr]=ExpSourj(b, StageDist, AgeArrayLower, AgeArrayUpper)
diagr=[];
X=0:100;
%b(3)=-1;
diagr1=ones(1,101);
diagr2=ones(1,101);
diagr3=ones(1,101);
diagr4=ones(1,101);

for i = 1:size(AgeArrayLower,2)
    value = (StageDist(1));
    if value == 0
        diagr1(AgeArrayLower(i):AgeArrayUpper(i)) = 100000000;
    else
        diagr1(AgeArrayLower(i):AgeArrayUpper(i))=(1/b(1,i))/value;
    end
    value = (StageDist(2)+StageDist(1));
    diagr2(AgeArrayLower(i):AgeArrayUpper(i))=(1/b(1,i))/value;
    value = (StageDist(3)+StageDist(2)+StageDist(1));
    diagr3(AgeArrayLower(i):AgeArrayUpper(i))=(1/b(1,i))/value;
    value = (StageDist(4)+StageDist(3)+StageDist(2)+StageDist(1));
    diagr4(AgeArrayLower(i):AgeArrayUpper(i))=(1/b(1,i))/value;
end

%     value = (StageDist(1));
%     if value == 0
%         diagr1(AgeArrayLower(i):AgeArrayUpper(i)) = 100000000;
%     else
%      diagr1(AgeArrayLower(i):AgeArrayUpper(i))=(1/b(1,i))/value;
%     end
%     value = (StageDist(2));
%     diagr2(AgeArrayLower(i):AgeArrayUpper(i))=(1/b(1,i))/value;
%     value = (StageDist(3));
%     diagr3(AgeArrayLower(i):AgeArrayUpper(i))=(1/b(1,i))/value;
%     value = (StageDist(4));
%     diagr4(AgeArrayLower(i):AgeArrayUpper(i))=(1/b(1,i))/value;
%     
% end


%%%%%%%%%%%%%%%%%
%%%LOGISTIC CURVE%%%
% %%%%%%%%%%%%%%%%%
% diagr1=b(4)./(1+(b(1)*(exp(b(2)*X))));%*(exp(-b(3)*1));
% diagr1(80:end)= diagr1(80);
% 
% value = (StageDist(1));
% %value= value^2;
% diagr2=(b(4)/value)./(1+(b(1)*(exp(-b(2)*X))));%*(exp(-b(3)*2));
% diagr2(80:end)= diagr2(80);
% 
% value = (StageDist(2));%+StageDist(1));
% %value= value^2;
% diagr3=(b(4)/value)./(1+(b(1)*(exp(-b(2)*X))));%*(exp(-b(3)*3));
% diagr3(80:end)= diagr3(80);
% 
% value = (StageDist(3));%+ StageDist(2)+ StageDist(1));
% %value= value^2;
% diagr4=(b(4)/value)./(1+(b(1)*(exp(-b(2)*X))));%*(exp(-b(3)*4));
% diagr4(80:end)= diagr4(80);
% % 

%%%%%%%%%%%%%%%%%
%%%EXPONENTIAL CURVE%%%
%%%%%%%%%%%%%%%%%
% diagr1=b(1)*(exp(b(2)*X));%*(exp(-b(3)*1));
% diagr1(80:end)= diagr1(80);
% 
% value = (StageDist(1));
% %value= value^2;
% diagr2=(b(1)/value)*(exp(b(2)*X));%*(exp(-b(3)*2));
% diagr2(80:end)= diagr2(80);
% 
% value = (StageDist(2)+StageDist(1));
% %value= value^2;
% diagr3=(b(1)/value)*(exp(b(2)*X));%*(exp(-b(3)*3));
% diagr3(80:end)= diagr3(80);
% 
% value = (StageDist(3)+ StageDist(2)+ StageDist(1));
% %value= value^2;
% diagr4=(b(1)/value)*(exp(b(2)*X));%*(exp(-b(3)*4));
% diagr4(80:end)= diagr4(80);

diagr=[diagr1;diagr2;diagr3;diagr4]';
end

function [diagr]=ExpSourjX(b,X, StageDist)
diagr=[];
%b(3)=-0.3;

%%%underdevloped countries where there is no screening diagnostic rate is likely to be highre in late stages than in early stages due to absence of screening 
for i = 1: size(X,2)
    value = (StageDist(1));
    if value == 0
        diagr1(i) = 100000000;
    else
        diagr1(i)=(1/b(1,i))/value;
    end
value = (StageDist(2)+StageDist(1));
diagr2(i)=(1/b(1,i))/value;
value = (StageDist(3)+StageDist(2)+StageDist(1));
diagr3(i)=(1/b(1,i))/value;
value = (StageDist(4)+StageDist(3)+StageDist(2)+StageDist(1));
diagr4(i)=(1/b(1,i))/value;
end

%%%US
% for i = 1: size(X,2)
%     value = (StageDist(1));
%     if value == 0
%         diagr1(i) = 100000000;
%     else
%         diagr1(i)=(1/b(1,i))/value;
%     end
%     value = (StageDist(2));
%     diagr2(i)=(1/b(1,i))/value;
%     value = (StageDist(3));
%     diagr3(i)=(1/b(1,i))/value;
%     value = (StageDist(4));
%     diagr4(i)=(1/b(1,i))/value;
% end

%%%%%%%%%%%%%%%%%
%%%LOGISTIC CURVE%%%
% %%%%%%%%%%%%%%%%%
% diagr1=b(4)./(1+(b(1)*(exp(b(2)*X))));%*(exp(-b(3)*1));
% diagr1(end)= diagr1(end-1);
% 
% value = (StageDist(1));
% %value= value^2;
% diagr2=(b(4)/value)./(1+(b(1)*(exp(-b(2)*X))));%*(exp(-b(3)*2));
% diagr2(end)= diagr2(end-1);
% 
% value = (StageDist(2));%+StageDist(1));
% %value= value^2;
% diagr3=(b(4)/value)./(1+(b(1)*(exp(-b(2)*X))));%*(exp(-b(3)*3));
% diagr3(end)= diagr3(end-1);
% 
% value = (StageDist(3));%+ StageDist(2)+ StageDist(1));
% %value= value^2;
% diagr4=(b(4)/value)./(1+(b(1)*(exp(-b(2)*X))));%*(exp(-b(3)*4));
% diagr4(end)= diagr4(end-1);

%%%%%%%%%%%%%%%%%
%%%EXPONENTIAL CURVE%%%
%%%%%%%%%%%%%%%%%
% diagr1=b(1)*(exp(b(2)*X));%*(exp(-b(3)*1));
% diagr1(end)= diagr1(end-1);
% 
% value = (StageDist(1));
% %value= value^2;
% diagr2=(b(1)/value)*(exp(b(2)*X));%*(exp(-b(3)*2));
% diagr2(end)= diagr2(end-1);
% 
% value = (StageDist(2)+StageDist(1));
% %value= value^2;
% diagr3=(b(1)/value)*(exp(b(2)*X));%*(exp(-b(3)*3));
% diagr3(end)= diagr3(end-1);
% 
% value = (StageDist(3)+ StageDist(2)+ StageDist(1));
% %value= value^2;
% diagr4=(b(1)/value)*(exp(b(2)*X));%*(exp(-b(3)*4));
% diagr4(end)= diagr4(end-1);

diagr=[diagr1;diagr2;diagr3;diagr4]';
%diagr=max(diagr,0.5);

end

function [Z,diagr,PrevS,PrevC]=PopSim(fitPAR,X,PAR,Incidence,OnsetRatesF,DwellTimesF,AgeArrayLower,AgeArrayUpper, Mu_Tx, BaseMortality, PlotModel, age, store_beta,Y,state_c)

% [~,~,TEMPMARKOV1]=xlsread('Template for OutputData Markov.xls','Total Pop');
% [~,~,TEMPMARKOV2]=xlsread('Template for OutputData Markov.xls','prev & prevA');
% [~,~,TEMPMARKOV3]=xlsread('Template for OutputData Markov.xls','Onset Rate');
% [~,~,TEMPMARKOV4]=xlsread('Template for OutputData Markov.xls','CFR');
%  [~,~,TEMPMARKOV5]=xlsread('Template for OutputData Markov.xls','Time to diagnosis');
%  [~,~,TEMPMARKOV6]=xlsread('Template for OutputData Markov.xls','Dwell Time');
%  [~,~,TEMPMARKOV7]=xlsread('Template for OutputData Markov.xls','Stage at diagnosis');
%  [~,~,TEMPMARKOV8]=xlsread('Template for OutputData Markov.xls','DiagsA');


Z=[];
PrevS=[];
PrevC=[];

AgeArrayMid=(AgeArrayUpper+AgeArrayLower)/2;

StageDist=PAR{6};

%% calculating diagnostic rates for each age ( 1 to 101 )

% PREVIOUS CALCULATIONS : In this case, we are fitting values to incidence 
% in each age and NOT fitting to incidence in each age and each state
%     if age < 13
%         if length(fitPAR)<2
%             current_beta = fitPAR;
%             fitPAR = store_beta;
%             fitPAR(1,age) = current_beta;
%         end
%     end  
%     [diagrT]=ExpSourj(fitPAR,StageDist, AgeArrayLower, AgeArrayUpper);
%     diagr = 1./diagrT

% NEW CALCULATIONS: In this case we are fitting values to incidence in each
% age and stage of cancer.
    diagr = zeros(4,101);
    % converting store_beta matrix to 4x101 size
    for ex_p = 1:length(AgeArrayLower)
        for c_s = 1:4
            diagr(c_s,AgeArrayLower(1,ex_p):AgeArrayUpper(1,ex_p)) = store_beta(c_s,ex_p);
        end
    end
    % when age is not equal to 10^5, it means we are in curve fitting
    % process
    if age ~= 10^5
        if age <= 6
            % odd age groups are of 4 years
            if mod(age,2) ~= 0
                diagr(state_c,(5*(age-1) + 1):(5*(age-1) + 1) + 3) = fitPAR;
            % even age groups are of 6 years
            elseif mod(age,2) == 0
                diagr(state_c,(5*(age-1)):(5*(age-1)) + 5) = fitPAR;
            end
        elseif age <= 11
            diagr(state_c,(5*(2*age-8)):(5*(2*age-8)) + 9) = fitPAR;
        else
            diagr(state_c,80:101) = fitPAR;
        end
    % when age = 10^5, it means that the curve fitting is done and we are
    % producing results/graphs for optimal answer
    elseif age == 10^5
        % pop simulation is using diagrT and diagr is output of the
        % function
        diagrT = fitPAR;
        for i=1:12
            if i <= 6
                if mod(i,2) ~= 0
                    for j=1:4
                        diagr(j,(5*(i-1) + 1):(5*(i-1) + 1) + 3) = diagrT(j,i);
                    end
                elseif mod(i,2) == 0
                    for j=1:4
                        diagr(j,(5*(i-1)):(5*(i-1)) + 5) = diagrT(j,i);
                    end
                end
            elseif i <= 11
                for j=1:4
                    diagr(j,(5*(2*i-8)):(5*(2*i-8)) + 9) = diagrT(j,i);
                end
            else
                for j=1:4
                    diagr(j,80:101) = diagrT(j,i);
                end
            end
        end
    end
                
    diagr = diagr';

%%        

if StageDist(1) == 0
    diagr(:,1) = 0;
    diagrT(:,1) = 0;
end

% diagr(:,1) =smooth(diagr(:,1) );
% diagr(:,2) =smooth(diagr(:,2) );
% diagr(:,3) =smooth(diagr(:,3) );
% diagr(:,4) =smooth(diagr(:,4) );
% % diagr(2:4) = StageDist(1:3)./diagrT(2:4); 
 
PrevA=zeros(size(AgeArrayLower));
%IncidenceF=IncidenceF';
OnsetRatesF=OnsetRatesF';
lambda=(1./DwellTimesF)';

probAgeA=PAR{1};probAgeA(81:end)=probAgeA(81)/length(probAgeA(81:end));
mu=PAR{2}';mu(81:end)=mu(80)+1/5;
BirthsF=PAR{3};%female pop only
prev=PAR{4};
prevN=prev/sum(prev);
muC =smooth(PAR{5}');
StageDist=PAR{6};

agemax=101;
m_a=length(probAgeA);
%age,state, duration of disease
pop0=10000;%work with pop of standard size

POP=zeros(m_a,9); % "9" are states, 1=healthy, 2=Pre-clinical0, 3=pre-clinical1 ....
POP(1:m_a,1)=probAgeA*pop0;

t_max=500;
CFR = zeros(m_a,t_max);

%% matrix for tracking deaths in preclinical

dt=1/15;
for t=1:t_max
 
    %fast time loop
    Diags=zeros(length(POP),4,t_max);
    IncidenceE = zeros(101,4);
    death_pre_clinical_4 = zeros(101,1);
    incidence_clinical_4 = zeros(101,1);
    
    for t1=1:1/dt    
        %healhy pop    
        POP(:,1)=POP(:,1) + dt*(-mu-OnsetRatesF).*POP(:,1);

        %non-diagnosed cancer
        %POP(:,2)=POP(:,2)+dt*(-mu-lambda(:,1)).*POP(:,2) + dt*(OnsetRatesF).*POP(:,1);
        POP(:,2)=POP(:,2)+dt*(-mu-lambda(:,1)-diagr(:,1)).*POP(:,2) + dt*(OnsetRatesF).*POP(:,1);
        POP(:,3)=POP(:,3)+dt*(-mu-lambda(:,2)-diagr(:,2)).*POP(:,3) + dt*(lambda(:,1)).*POP(:,2);
        POP(:,4)=POP(:,4)+dt*(-mu-lambda(:,3)-diagr(:,3)).*POP(:,4) + dt*(lambda(:,2)).*POP(:,3);
        POP(:,5)=POP(:,5)+dt*(-mu-lambda(:,4)-diagr(:,4)).*POP(:,5) + dt*(lambda(:,3)).*POP(:,4);

        % Deaths in preclinical last stage
        death_pre_clinical_4(:,1) = death_pre_clinical_4(:,1) + POP(:,5).*(lambda(:,4))*dt;

        %diagnosed cancer
        IncidenceE(:,1) = IncidenceE(:,1) + dt*(diagr(:,1)).*POP(:,2);
        POP(:,6)=POP(:,6)+dt*(-mu ).*POP(:,6) - dt*BaseMortality(1,:)'.*POP(:,6)+ dt*(diagr(:,1)).*POP(:,2);
        CFR(:,t)= CFR(:,t)+ dt*(mu ).*POP(:,6)+ dt*BaseMortality(1,:)'.*POP(:,6);

        IncidenceE(:,2) = IncidenceE(:,2) + dt*(diagr(:,2)).*POP(:,3);
        POP(:,7)=POP(:,7)+dt*(-mu ).*POP(:,7) - dt*BaseMortality(2,:)'.*POP(:,7)+ dt*(diagr(:,2)).*POP(:,3);
        CFR(:,t)= CFR(:,t)+ dt*(mu ).*POP(:,7)+ dt*BaseMortality(2,:)'.*POP(:,7);

        IncidenceE(:,3) = IncidenceE(:,3) + dt*(diagr(:,3)).*POP(:,4);
        POP(:,8)=POP(:,8)+dt*(-mu ).*POP(:,8) - dt*BaseMortality(3,:)'.*POP(:,8)+ dt*(diagr(:,3)).*POP(:,4);
        CFR(:,t)= CFR(:,t)+dt*(mu ).*POP(:,8) + dt*BaseMortality(3,:)'.*POP(:,8);

        IncidenceE(:,4) = IncidenceE(:,4) + dt*(diagr(:,4)).*POP(:,5);
        POP(:,9)=POP(:,9)+dt*(-mu ).*POP(:,9) - dt*BaseMortality(4,:)'.*POP(:,9)+ dt*(diagr(:,4)).*POP(:,5);
        CFR(:,t)= CFR(:,t)+dt*(mu ).*POP(:,9) + dt*BaseMortality(4,:)'.*POP(:,9);
        % POP(:,7)=POP(:,7) - dt*BaseMortality(2,:)'.*POP(:,7)+ dt*(diagr(:,2)).*POP(:,3);
        % POP(:,8)=POP(:,8) - dt*BaseMortality(3,:)'.*POP(:,8)+ dt*(diagr(:,3)).*POP(:,4);
        % POP(:,9)=POP(:,9) - dt*BaseMortality(4,:)'.*POP(:,9)+ dt*(diagr(:,4)).*POP(:,5);


        % incidence in clinical 4
        incidence_clinical_4 = incidence_clinical_4 + dt*(diagr(:,4)).*POP(:,5);

        %if t > t_max - 100
        Diags(:,1,t)=Diags(:,1,t)+dt*( ((diagr(:,1)).*POP(:,2)) );
        Diags(:,2,t)=Diags(:,2,t)+dt*( ((diagr(:,2)).*POP(:,3)) );
        Diags(:,3,t)=Diags(:,3,t)+dt*( ((diagr(:,3)).*POP(:,4)) );
        Diags(:,4,t)=Diags(:,4,t)+dt*( ((diagr(:,4)).*POP(:,5)) );
        %end

    end

    DiagsAll=sum(Diags,2);

    %aging
    for I=1:9
    POP(2:end,I)=POP(1:end-1,I);
    end
    POP(1,:)=0;
    %must decide if births are to be based on population size
    Births=BirthsF*10000;%BirthsF*sum(POP(:));
    %Births=10000-sum(sum(POP));
    %Births=BirthsF*sum(POP(:));

    POP(1,1)=Births;



    DiagsA(1,t) = sum(Diags(:,1,t))/ sum(sum(Diags(:,:,t)));
    DiagsA(2,t) = sum(Diags(:,2,t))/ sum(sum(Diags(:,:,t)));
    DiagsA(3,t) = sum(Diags(:,3,t))/ sum(sum(Diags(:,:,t)));
    DiagsA(4,t) = sum(Diags(:,4,t))/ sum(sum(Diags(:,:,t)));
    
end
%prevalence of cancer by age
%average mort in each age group
for J=1:length(AgeArrayLower)
    indA=AgeArrayLower(J):AgeArrayUpper(J);
    
    inc=IncidenceE(indA,1:4); 
   
    p_a=POP(indA,:);
    p_c=POP(indA,[6:9]);
    
    % death_rate in distant stage
    death_distant_age_group = death_pre_clinical_4(indA,1);
    
    %pre clinical
    p_s1=POP(indA,[2]);
    p_s2=POP(indA,[3]);
    p_s3=POP(indA,[4]);
    p_s4=POP(indA,[5]);
   
    %diagnosed/clinical
    p_c1=POP(indA,[6]);
    p_c2=POP(indA,[7]);
    p_c3=POP(indA,[8]);
    p_c4=POP(indA,[9]);
    
    if age ~= 10^5
        incidenceCalc(J) = sum(inc(:,state_c))/sum(p_a(:));
    elseif age == 10^5
        for state_c = 1:4
            incidenceCalc(state_c,J) = sum(inc(:,state_c))/sum(p_a(:));
        end
       
    end
                
    PrevA(J)=sum(p_c(:))/sum(p_a(:));
    
    %prev cancer pre clinical
    PrevS(J,1)=sum(p_s1(:))/sum(p_a(:));
    PrevS(J,2)=sum(p_s2(:))/sum(p_a(:));
    PrevS(J,3)=sum(p_s3(:))/sum(p_a(:));
    PrevS(J,4)=sum(p_s4(:))/sum(p_a(:));
    
    %prev cancer clinical
    PrevC(J,1)=sum(p_c1(:))/sum(p_a(:));
    PrevC(J,2)=sum(p_c2(:))/sum(p_a(:));
    PrevC(J,3)=sum(p_c3(:))/sum(p_a(:));
    PrevC(J,4)=sum(p_c4(:))/sum(p_a(:));
    
    
end

%Diags=Diags/sum(Diags)*100;
%Z=[PrevA*100];If fitting to prevalence
if age >1000
    Z=incidenceCalc;
   DCO =sum(death_pre_clinical_4) / (sum(death_pre_clinical_4)+sum(sum(IncidenceE)));
else
    % incidenceCalc (simulated incidence) already considers state of
    % cancer, see line ~713.
    Z=incidenceCalc(1,age)*100;
    DCO =sum(death_pre_clinical_4) / (sum(death_pre_clinical_4)+sum(sum(IncidenceE)));
    % calculating error
    err = Z-Y(state_c,age);
    

end
    
%PlotModel=1;
if(PlotModel==1)
    
close all
figure
subplot(4,1,1)
plot(POP(:,1))
title('total pop')

subplot(4,1,2)
plot(sum(POP(:,[2:5]),2) )
title('cancer not diagnosed')

subplot(4,1,3)
plot(sum(POP(:,[6:9]),2) )
title('total pop')
title('cancer diagnosed')

subplot(4,1,4)
plot(OnsetRatesF )
title('total pop')
title('Onset Rate')


%PrevA=Z(1:length(prev))
figure
hold on
plot(AgeArrayMid,prev*100,'k*')
plot(AgeArrayMid,PrevA,'r-')
title('Prevalence of Diagnosed cases by age')

% incidence fit plot 
    if age ~= 10^5
       
        % calculating actual incidence for each state of cancer
        Incidence_state = StageDist'*Incidence;

        figure
        hold on
        plot(AgeArrayMid,Incidence_state(state_c,:)*100,'k*')
        plot(AgeArrayMid,incidenceCalc*100,'r-')
        title('Incidence of Diagnosed cases by age')

    elseif age == 10^5

        % calculating actual incidence for each state of cancer
        Incidence_state = StageDist'*Incidence;

        for state_c = 1:4
            figure
            hold on
            plot(AgeArrayMid,Incidence_state(state_c,:)*100,'k*')
            plot(AgeArrayMid,incidenceCalc(state_c,:)*100,'r-')
            formatSpec = 'Incidence of Diagnosed cases by age for cancer-state %d';
            str = sprintf(formatSpec,state_c);
            title(str)
        end
    end


figure
hold on
plot(0:100,OnsetRatesF )

title('Onset Rate')

figure
hold on
plot(0:100,CFR(:,t))
plot(0:100,CFR(:,t-100))
plot(0:100,CFR(:,t-200))
plot(0:100,CFR(:,t-300))
title('CFR')



% figure
% hold on
% X=0:100;
% plot(X,diagrT(:,1),'b')
% plot(X,diagrT(:,2),'r')
% plot(X,diagrT(:,3),'g')
% plot(X,diagrT(:,4),'c')
% legend('in-situ', 'local', 'regional', 'distant');
% title('Time to diagnosis')

figure
hold on
X=0:100;
plot(X,DwellTimesF(1,:),'b')
plot(X,DwellTimesF(2,:),'r')
plot(X,DwellTimesF(3,:),'g')
plot(X,DwellTimesF(4,:),'c')
legend('in-situ', 'local', 'regional', 'distant');

title('Dwell Time')

figure
hold on
plot(1:4,DiagsA(1:4,t),'r-')
plot(1:4,StageDist(1:4),'k*')
title('Stage at diagnosis')


figure
hold on
subplot(4,1,1)
plot(1:500,DiagsA(2,:),'r-')
title('Local')

subplot(4,1,2)
plot(1:500,DiagsA(3,:),'r-')
title('Regional')

subplot(4,1,3)
plot(1:500,DiagsA(4,:),'r-')
title('Distant')




% CountryName = PAR{10};
% 
% [AA,BB,CC,DD,EE,FF,GG,HH] = MikeWriteTemplateMarkov(TEMPMARKOV1,TEMPMARKOV2,TEMPMARKOV3,TEMPMARKOV4,TEMPMARKOV5,TEMPMARKOV6,TEMPMARKOV7,TEMPMARKOV8,PAR{9},POP,OnsetRatesF,CFR,diagr,DwellTimesF,DiagsA,StageDist,AgeArrayMid,prev,PrevA,CountryName)
% 
% CountryIndex=PAR{9};
% 
%   xlswrite('Mike Markov Results.xls',AA,'Total Pop')
%   xlswrite('Mike Markov Results.xls',BB,'prev & prevA')
%   xlswrite('Mike Markov Results.xls',CC,'Onset Rate')
%   xlswrite('Mike Markov Results.xls',DD,'CFR')
%    xlswrite('Mike Markov Results.xls',EE,'Time to diagnosis')
%    xlswrite('Mike Markov Results.xls',FF,'Dwell Time')
%    xlswrite('Mike Markov Results.xls',GG,'Stage at diagnosis')
%    xlswrite('Mike Markov Results.xls',HH,'DiagsA')
%   
% 
%   'Mike Markov Results Written'
end

end


%%ESTIMATE INCIDENCE
function [incidenceArray] = IncidenceFunc(probAgeA,MortalityArray,prevalenceActual,cancerMortalityArray, AgeArrayLower,AgeArrayUpper)

 incidenceArray = zeros(1, size(AgeArrayLower,2));
 sumAgeArray= zeros(1, size(AgeArrayLower,2));
 for j = 1 : size (AgeArrayLower,2)
     
    sumAge = 0;
    mortSum = 0;
    for  i = AgeArrayLower(j) : AgeArrayUpper(j)
      sumAge = sumAge + probAgeA(i);
      mortSum = mortSum  + cancerMortalityArray(i);
    end
    
    sumAgeArray(j) = sumAge;
    mortSum = mortSum / (AgeArrayUpper(j) - AgeArrayLower(j) + 1);
    
    if j == 1
     incidenceArray(j) = (prevalenceActual(j) * sumAgeArray(j) * (mortSum + (1 / (AgeArrayUpper(j) - AgeArrayLower(j) + 1))))/ (sumAgeArray(j) - prevalenceActual(j));
    % incidenceArray(j) = max(0,(prevalenceActual(j) * sumAgeArray(j) * (mortSum + (1 / (AgeArrayUpper(j) - AgeArrayLower(j) + 1))))/ (sumAgeArray(j) - prevalenceActual(j)));
    else
     %incidenceArray(j) = max(incidenceArray(j-1),((prevalenceActual(j) * sumAgeArray(j) * (mortSum + (1 / (AgeArrayUpper(j) - AgeArrayLower(j) + 1)))) -  (prevalenceActual(j - 1) * sumAgeArray(j-1) * (1 / (AgeArrayUpper(j - 1) - AgeArrayLower(j -1) + 1))))/ (sumAgeArray(j) - prevalenceActual(j)));
     incidenceArray(j) = ((prevalenceActual(j) * sumAgeArray(j) * (mortSum + (1 / (AgeArrayUpper(j) - AgeArrayLower(j) + 1)))) -  (prevalenceActual(j - 1) * sumAgeArray(j-1) * (1 / (AgeArrayUpper(j - 1) - AgeArrayLower(j -1) + 1))))/ (sumAgeArray(j) - prevalenceActual(j));
     if incidenceArray(j) < 0
        incidenceArray(j) = incidenceArray(j-1);
     end
     
    end%if J==1
   
 end
end%IncidenceFunc

%%ESTIMATE Previous Stage 
function [OnsetRates] = PrevStageRates (CurrentStageRate,TimeToDiagnosis,StageDist,probAgeA,MortalityArray,prevalenceActual,cancerMortalityArray, AgeArrayLower,AgeArrayUpper)
numStage = size(StageDist,2);
numAgeGroups = size(AgeArrayLower,2);
minAgeCancer = 3;

%CurrentStageRate = incidence =(I_(D_a ))
%lamdaArray= lambda
%probAgeA = age distribution =(C_a in step 2 of algo)
%MortalityArray= mu
%AgeArrayLower, AgeArrayUpper =lower and upper bounds of age group

%prevalenceActual= do not need this anymore, but lets not delete it as yet 
%cancerMortalityArray = do not need this anymore, but lets not delete it as yet 

    
%%                    
                    % For analysis by age group %     
    
    % Here, all the matrices should in dimension of 1 X 12
    
    % Dimension of matrix "probAgeA" is 1 X 100, now converting this matrix
    % into 1 X 12 matrix. This has to done in accordance with predesigned
    % age groups, "AgeArrayLower" and "AgeArrayUpper", these
    % two matrices define age groups. With help of these matrices, elements
    % of probAgeA are added together (Simple "OR" operation of
    % probability). By doing this we'll get prob of a person being in that
    % group.

    a_max_group=size(AgeArrayLower,2);
    
    % defiing matrix for storing values
    C_a=zeros(1,size(AgeArrayLower,2));
    
    for i=(1:size(AgeArrayLower,2))
        
        % Reading values of lower and upper limit of "i th" age group
        l_limit=AgeArrayLower(1,i);
        u_limit=AgeArrayUpper(1,i);
        
        % with help l_limit and u_limit, we'll do summation over those
        % indexes in probAgeA
       C_a(1,i)=sum(probAgeA(1,l_limit:u_limit));
        
    end
    
    A_a=C_a;
                                    % here, we have made an assumption that
                                    % A_a = C_a, because probability of a
                                    % person being diagnosed in age "a" is
                                    % very low
    
    % Next matrix which is not in 1 X 12 dimension is MortalityArray. this
    % matrix contains rates and not prbability values, hence we cannot
    % simply add them together. We have to calculate expected value first
    % and then add them together.
    
    % calculating expected value of mortality rate
    mu_exp=probAgeA.*MortalityArray;
    
    % defining matrix for storing values
    mu=zeros(1,size(AgeArrayLower,2));
    
  
    
    % converting mu_exp in 1 X 12 matrix
    for i=(1:size(AgeArrayLower,2))
        
        % Reading values of lower and upper limit of "i th" age group
        l_limit=AgeArrayLower(1,i);
        u_limit=AgeArrayUpper(1,i);
        
        % with help l_limit and u_limit, we'll do summation over those
        % indexes in mu_exp
        mu(1,i)=sum(mu_exp(1,l_limit:u_limit))/sum(probAgeA(l_limit:u_limit));
        
    end
    
    % defining lamda
    lamda=zeros(numStage,size(AgeArrayLower,2));
    for i=(1:numStage)
        for n=(1:numAgeGroups)
            lamda(i,n)=1/TimeToDiagnosis(i,n);
        end
    end
    
    % defining matix for long term prob of being in Healthy state and being
    % in Undiagnosed state and prob of transition from H to U
  steady_h=zeros(1,size(AgeArrayLower,2));
  steady_h(1,1:minAgeCancer-1)=A_a(1,1:minAgeCancer-1);
  
   steady_u=zeros(1,size(AgeArrayLower,2));
  probability_hu=zeros(1,size(AgeArrayLower,2));
    
    % defining incidence matrix
    I_D_a=CurrentStageRate;
    I_D_a(1,1:minAgeCancer-1)=0;
    
    % creating array for storing mid age values
    mid_age_group=zeros(1,a_max_group);
    for i=(1:numAgeGroups)
        mid_age_group(1,i)=(AgeArrayLower(1,i)+AgeArrayUpper(1,i))/2;
    end
    
    % defining prob of being in one of 4 cancer states
  s_i=StageDist;
    
    % onset rate matrix
    OnsetRates=zeros(1,a_max_group);
  %% Plotting some parameters
    
    figure
    x=(1:numAgeGroups);
    plot(x,I_D_a);
    ylabel('Incidence');
    figure
    plot(x,A_a);
    ylabel('A_a/C_a');
    plot(x,mu);
    ylabel('Mortalities');
    figure
    for i=(1:numStage)
        plot(x,lamda(i,:));
        hold on
    end
    ylabel('lamda');
    figure
    for i=(1:numStage)
        plot(x,TimeToDiagnosis(i,:));
        hold on
    end
    ylabel('Time to Diagnosis');
    
    for a = minAgeCancer+1:numAgeGroups   % START ALGORITHM 
      
   %%STEP 1
        IC = I_D_a(1,a) * C_a(1,a);
        
        sumKloop =0;
        for k = minAgeCancer:a-1
                      
           sumIloopNumerator = 0;
%           t = mid_age_group(1,a)-mid_age_group(1,k);
%            for i = 1:4
%                 sumIloopNumerator = sumIloopNumerator + s_i(1,i)* (  (1- exp(-t*lamda(i,k))) - (1- exp(-(t-1)*lamda(i,k))) );
%            end

           t = (AgeArrayUpper(1,a)-mid_age_group(1,k));
           t_minus_1 = (AgeArrayLower(1,a)-mid_age_group(1,k));
            for i = 1:numStage
                sumIloopNumerator = sumIloopNumerator + s_i(1,i)* (  (1- exp(-t*lamda(i,k))) - (1- exp(-(t_minus_1)*lamda(i,k))) );
            end
        
            muProduct=1;
%             for j = k:a
%              muProduct = muProduct * exp(-(mid_age_group(j)-mid_age_group(j-1)) *mu(1,j));
%             end
            for j = round(mid_age_group(1,k)):AgeArrayUpper(1,a)
             muProduct = muProduct * exp(-MortalityArray(1,j));
            end
            sumKloop = sumKloop + steady_h(1,k) * probability_hu(1,k) * sumIloopNumerator * muProduct;
        end
        
        
        sumIloopDenominator = 0;
        t = (AgeArrayUpper(1,a)-AgeArrayLower(1,a))/2;
        for i = 1:numStage
                sumIloopDenominator = sumIloopDenominator + s_i(1,i)* (1- exp(-t*lamda(i,a)));
        end
        denominator = A_a(1,a) * sumIloopDenominator * exp(-t * mu(1,a)) - IC;
        probability_hu(1,a) = (IC - sumKloop)/ denominator;
    %% STEP 2: does not look like we use this anywhere
    
    
    %% STEP 3
     sumKloop =0;
        for k = minAgeCancer:a-1
                      
           sumIloopNumerator = 0;
%           t = mid_age_group(1,a)-mid_age_group(1,k);
             t = AgeArrayUpper(1,a)-mid_age_group(1,k);
            for i = 1:numStage
                sumIloopNumerator = sumIloopNumerator + s_i(1,i)*  exp(-t*lamda(i,k)) ;
            end
       
            muProduct=1;
%             for j = k:a
%              muProduct = muProduct * exp(-(mid_age_group(j)-mid_age_group(j-1)) *mu(1,j));
%             end
            for j = round(mid_age_group(1,k)):AgeArrayUpper(1,a)
             muProduct = muProduct * exp(-MortalityArray(1,j));
            end
            sumKloop = sumKloop + steady_h(1,k) * probability_hu(1,k) * sumIloopNumerator * muProduct;
        end
        
      steady_h(1,a) = (A_a(1,a) - sumKloop) / (1+probability_hu(1,a));
    
    
    
    end %% END ALGORITHM
   
    
 OnsetRates = -log(1- probability_hu);

 onset_ref=[0,0,0,6.75E-05,0.000171533,0.000335284,0.000579637,0.001671398,0.002955115,0.003824575,0.004989599,0.004989599];
    x=(1:101);
    figure
    plot(mid_age_group(1:12),onset_ref,'--',mid_age_group,OnsetRates);
    legend('Reference values of Onset','Estimated');
    title('Analysis for Pre-Designed Age Groups');
end