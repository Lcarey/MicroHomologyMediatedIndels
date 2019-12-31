%% Figure to plot the frequency of duplication events across ssp1

%% load data
DATADIR = '~/CareyLab/Projects/2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU/DataFromCluster/' ; 
FIGUREDIR = '~/CareyLab/Projects/2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU/FIGURES/' ; 

T = readtable( [ DATADIR  'PombeAmpliconSeq_E2_alltsvs.txt'] ,'TreatAsEmpty','-'); 
T.Properties.VariableNames = {'lib' 'chr' 's1' 'e1' 's2' 'e2' 'ReadDepth' 'DupCounts' 'DupFreq' 'CollapseCounts' 'CollapseFreq' };
T.DupCounts(isnan(T.DupCounts))=0 ;
T.DupFreq(isnan(T.DupFreq))=0 ;
T.HasDup = T.DupCounts>0 ;
T.MHlen = T.e1-T.s1+1 ;

idx_to_keep = strcmp( T.chr , 'ssp1') ;
T = T( idx_to_keep , :) ; 
T.midpt = mean([T.s1 T.e2],2);
T = sortrows(T ,'midpt');

U = unstack( T , {'DupCounts' 'DupFreq'} , 'lib' ,  'GroupingVariables'  , {'s1' 's2' 'e1' 'e2' 'midpt'} ) ; 

%%
idx1 = strcmp(T.lib,'lib_1');
fh = figure('units','centimeters','position',[5 5 20 9]) ;
hold on ; 
plot( T.midpt(idx1) , T.DupCounts(idx1) , '-')
plot( T.midpt(~idx1) , -1*T.DupCounts(~idx1) , '-')
axis tight ;
legend({'rep 1' 'rep 2'})
ylim([-20 20])
ylabel('# of reads supporting each duplication')
set(gca,'ytick',-20:10:20)
set(gca,'yticklabel',[20 10 0 10 20])
xlabel('Position along ssp1 (nt)')

%% correlation between the two
X = U.DupCounts_lib_1 ; Y = U.DupCounts_lib_2 ; 
%X = U.DupFreq_lib_1 ;  Y = U.DupFreq_lib_2 ; 

X(X>50)=50;
Y(Y>50)=50;
%X(X==0)=random('uniform',-0.05,0.05,sum(X==0),1);
%Y(Y==0)=random('uniform',-0.05,0.05,sum(Y==0),1);

X(X==0)=0.5;
Y(Y==0)=0.5;
fh = figure('units','centimeters','position',[5 5 9 9]) ;
hold on ; 
sh = scatter( X , Y , 'ok' , 'MarkerFaceColor',[.7 .7 .7]);
set(gca,'xscale','log');
set(gca,'yscale','log');
set(gca,'xtick',[0.5 1 2 5 10 20 50])
set(gca,'ytick',get(gca,'xtick'))
set(gca,'xticklabel',[0 1 2 5 10 20 50])
set(gca,'yticklabel',get(gca,'xticklabel'))
axis tight ;
xlabel('# of reads in rep 1')
ylabel('# of reads in rep 2')
[r,p]=corr(U.DupCounts_lib_1,U.DupCounts_lib_2,'type','Spearman')
