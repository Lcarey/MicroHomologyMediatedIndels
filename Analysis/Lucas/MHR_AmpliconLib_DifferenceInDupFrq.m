DATADIR = '/Users/lcarey/SynologyDrive/Projects/2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU/DataFromCluster/';
T = readtable( [ DATADIR 'PombeAmpliconSeq_E4_alltsvs.txt' ],'TreatAsEmpty','-'); T.Properties.VariableNames = {'lib' 'chr' 's1' 'e1' 's2' 'e2' 'ReadDepth' 'Ndup' 'Ncol' 'DupFrq' 'ColFrq'};
%T = readtable( [ DATADIR 'PombeAmpliconSeq_E3_alltsvs.txt' ],'TreatAsEmpty','-');T.Properties.VariableNames = {'chr' 's1' 'e1' 's2' 'e2' 'ReadDepth' 'Ndup' 'Ncol' 'DupFrq' 'ColFrq'};

T.DupFrq(isnan(T.DupFrq))=0;
T.Ndup(isnan(T.Ndup))=0;
T.Ncol = []; T.ColFrq = [] ;

T.MHlen = T.e1 - T.s1 + 1 ;
T.InterMHDistance = T.s2 - T.e1 + 1 ;
T.InterMHDistanceR = round(T.InterMHDistance./10)*10 ; 

%T=T(T.MHlen>=6,:);

%% Amplicon 4 library
WT = mean( [T.DupFrq( strcmp(T.chr,'ssp1.dup.1')) T.DupFrq( strcmp(T.chr,'ssp1.dup.2')) T.DupFrq( strcmp(T.chr,'ssp1.dup.3'))],2);  
MRC1 = mean( [T.DupFrq( strcmp(T.chr,'mrc1.ssp1.dup.2')) T.DupFrq( strcmp(T.chr,'mrc1.ssp1.dup.3')) T.DupFrq( strcmp(T.chr,'mrc1.ssp1.dup.4')) T.DupFrq( strcmp(T.chr,'mrc1.ssp1.dup.5'))],2);  
RAD52 = mean( [T.DupFrq( strcmp(T.chr,'rad52.ssp1.dup.1')) T.DupFrq( strcmp(T.chr,'rad52.ssp1.dup.2')) T.DupFrq( strcmp(T.chr,'rad52.ssp1.dup.3'))],2);  
RAD50 = mean( [T.DupFrq( strcmp(T.chr,'rad50.ssp1.dup.2')) T.DupFrq( strcmp(T.chr,'rad50.ssp1.dup.3'))],2);  
MRC1 = MRC1+0.1 ; RAD52 = RAD52+0.1 ; RAD50 = RAD50+0.1 ; WT = WT+0.1 ;
data = [ log2(MRC1./WT)  log2(RAD52./WT) log2(RAD50./WT)  ] ;
data = data( ~all(data==0,2) , :);

[~,p] = ttest(WT,MRC1) ; fprintf('MRC1 p = %0.05f\n' , p); 
[~,p] = ttest(WT,RAD52) ; fprintf('RAD52 p = %0.05f\n' , p); 
[~,p] = ttest(WT,RAD50) ; fprintf('RAD50 p = %0.05f\n' , p); 
fprintf('MRC1  > 0 = %0.01f\n' , mean(MRC1>0)*100 ); 
fprintf('RAD52 > 0 = %0.01f\n' , mean(RAD52>0)*100 ); 
fprintf('RAD50 > 0 = %0.01f\n' , mean(RAD50>0)*100 ); 
fprintf('WT    > 0 = %0.01f\n' , mean(WT>0)*100 ); 

figure; 
subplot(2,2,1)
plot( WT , MRC1 , 'ok','MarkerFaceColor','b')
set(gca,'yscale','log')
set(gca,'xscale','log')
line(xlim,xlim)
xlabel('wild-type')
ylabel('mrc1 KO')

subplot(2,2,2)
plot( WT , RAD52 , 'ok','MarkerFaceColor','r')
set(gca,'yscale','log')
set(gca,'xscale','log')
line(xlim,xlim)
xlabel('wild-type')
ylabel('rad52 KO')

subplot(2,2,3)
hold on ;
[f,x]=ecdf(data(:,1));plot(x,f,'-','DisplayName','mrc1');
[f,x]=ecdf(data(:,2));plot(x,f,'-','DisplayName','rad52');
[f,x]=ecdf(data(:,3));plot(x,f,'-','DisplayName','rad50');
xlabel('ratio (mutant/wt)')
ylabel('fraction of MHpairs')
legend('location','se')
set(gca,'ytick',0:.1:1)
grid on ;
xlim([-3 3])


subplot(2,2,4)
hold on ; 

boxplot(data ,'symbol','')
set(gca,'xtick',1:3)
set(gca,'xticklabel',{'mrc1' 'rad52' 'rad50'})
line(xlim,[0 0])
ylim([-4 4])
ylabel('mutant / WT')


%
%% %% Amplicon 3 library
WT = mean( [T.DupFrq( strcmp(T.chr,'SPCC1235.01.1')) T.DupFrq( strcmp(T.chr,'SPCC1235.01.2'))],2);  
CDS1 = mean( [T.DupFrq( strcmp(T.chr,'SPCC1235.01.cds1d.1')) T.DupFrq( strcmp(T.chr,'SPCC1235.01.cds1d.2'))],2);  
RAD52 = mean( [T.DupFrq( strcmp(T.chr,'SPCC1235.01.rad3d.1' )) T.DupFrq( strcmp(T.chr,'SPCC1235.01.rad3d.2' ))],2);  
RAD50 = mean( [T.DupFrq( strcmp(T.chr,'SPCC1235.01.rad50d.1')) T.DupFrq( strcmp(T.chr,'SPCC1235.01.rad50d.2'))],2);  
RAD27 = mean( [T.DupFrq( strcmp(T.chr,'SPCC1235.01.rad27d.1')) T.DupFrq( strcmp(T.chr,'SPCC1235.01.rad27d.2'))],2);  

CDS1 = CDS1+0.1 ; RAD52 = RAD52+0.1 ; RAD50 = RAD50+0.1 ; WT = WT+0.1 ; RAD27 = RAD27+0.1 ;
data = [ log2(CDS1./WT)  log2(RAD52./WT) log2(RAD50./WT)  log2(RAD27./WT)  ] ;
data = data( ~all(data==0,2) , :);

[~,p] = ttest(WT,CDS1) ; fprintf('CDS1 p = %0.05f\n' , p); 
[~,p] = ttest(WT,RAD52) ; fprintf('RAD52 p = %0.05f\n' , p); 
[~,p] = ttest(WT,RAD50) ; fprintf('RAD50 p = %0.05f\n' , p); 
fprintf('CDS1  > 0 = %0.01f\n' , mean(MRC1>0)*100 ); 
fprintf('RAD52 > 0 = %0.01f\n' , mean(RAD52>0)*100 ); 
fprintf('RAD50 > 0 = %0.01f\n' , mean(RAD50>0)*100 ); 
fprintf('WT    > 0 = %0.01f\n' , mean(WT>0)*100 ); 

figure; 
subtitle('SPCC1235 Amp3 lib')
subplot(2,2,1)
plot( WT , CDS1 , 'ok','MarkerFaceColor','b')
set(gca,'yscale','log')
set(gca,'xscale','log')
line(xlim,xlim)
xlabel('wild-type')
ylabel('mrc1 KO')

subplot(2,2,2)
plot( WT , RAD27 , 'ok','MarkerFaceColor','r')
set(gca,'yscale','log')
set(gca,'xscale','log')
line(xlim,xlim)
xlabel('wild-type')
ylabel('rad27 KO')

subplot(2,2,3)
hold on ;
[f,x]=ecdf(data(:,1));plot(x,f,'-','DisplayName','cds1');
[f,x]=ecdf(data(:,2));plot(x,f,'-','DisplayName','rad52');
[f,x]=ecdf(data(:,3));plot(x,f,'-','DisplayName','rad50');
[f,x]=ecdf(data(:,4));plot(x,f,'-','DisplayName','rad27');
xlabel('ratio (mutant/wt)')
ylabel('fraction of MHpairs')
legend('location','se')
set(gca,'ytick',0:.1:1)
grid on ;
xlim([-3 3])


subplot(2,2,4)
hold on ; 

boxplot(data ,'symbol','')
set(gca,'xtick',1:size(data,2))
set(gca,'xticklabel',{'cds1' 'rad52' 'rad50' 'rad27'})
line(xlim,[0 0])
ylim([-4 4])
ylabel('mutant / WT')