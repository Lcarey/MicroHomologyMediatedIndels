% We identified three MTDs in the SAGA complex histone acetyltransferase catalytic subunit gcn5, far higher than expected by 
% chance based on the gene length (p=@) or the number of MHPs (p=@). The model predicts gcn5 to be in the top @% of genes for MTD frequency, 
%  suggesting that MTDs in gcn5 should be frequent. Indeed, whole-genome sequencing of 16 suppressors of htb1G52D identified MTDs in gcn5, as well as in 
%  ubp8, where we also observed a single MTD in our high-coverage sequencing data. 
%  These results suggest that MTDs likely exist at high frequency in most genes, 
%  and are frequently the raw material on which natural selection acts. 

DATADIR = '/Users/lcarey/SynologyDrive/Projects/2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU/Sarah/MH_project/Manuscript-todo/processeddata/' ;
T = readtable([DATADIR 'MHRSumPreinGene.txt' ],'TreatAsEmpty','.');
T.Properties.VariableNames = { 'chr'	'start1'	'end2'	'gene'	'systematic_name'	'sumObs'	'sumPre'	'sumFloat'}; 
T.GeneLength = abs(T.start1 - T.end2) ./ 1000 ; 
T.SumFloat_N = T.sumFloat ./ T.GeneLength ; 
T.sumObs_N = T.sumObs ./ T.GeneLength ; 
T.sumPre_N = T.sumObs ./ T.GeneLength ; 
T = sortrows(T,'sumFloat','ascend');
T.sumFloat_Rank = (1:height(T))' ;

T = sortrows(T,'SumFloat_N','ascend');
T.sumFloat_N_Rank = (1:height(T))' ;
%%
Yo = T.sumObs ; Yo(Yo>100)=100; 
Yp = T.sumObs ; Yp(Yp>100)=100; 

fh = figure('units','centimeters','position',[5 5 8 8]);
hold on ;
histogram( Yo , -0.5:max(Yo))
set(gca,'yscale','log')
%%
Yobs = T.sumObs; 
Ypred = T.sumPre ; 
idx = find(strcmp(T.gene,'gcn5'));

pctile_pred = mean( T.sumPre>T.sumPre(idx))*100
pctile_obs = mean( T.sumObs>T.sumObs(idx))*100

pctile_pred_N = mean( T.sumPre_N>T.sumPre_N(idx))*100
pctile_obs_N = mean( T.sumObs_N>T.sumObs_N(idx))*100
%%


gcn5_lg_txt = sprintf('\\it{gcn5} , %0.0f%% obs, , %0.0f%% pred' , pctile_obs , pctile_pred) ; 

fh = figure('units','centimeters','position',[5 5 8 8]);
hold on ;
[f,x]=ecdf( T.sumObs );
plot(x,f*100,'-k','LineWidth',2,'DisplayName','Observed MTDs')

[f,x] = ecdf( T.sumPre ) ; 
plot(x,f*100,'-r','LineWidth',2,'DisplayName','Predicted MTDs')

line( [ T.sumObs(idx) T.sumObs(idx)] , ylim ,'Color','b','DisplayName', gcn5_lg_txt ,'LineWidth',3)
xlabel('# of MTDs per gene')

grid on ;
set(gca,'ytick',0:10:100)
xlim([0 5])
legend('location','se');
ylabel('% of genes')



fh = figure('units','centimeters','position',[5 5 8 8]);
hold on ;
histogram( T.sumObs , -0.5:15 ,'FaceColor','k');
%plot(x,f*100,'-k','LineWidth',2,'DisplayName','Observed MTDs')

histogram( T.sumPre , -0.5:15);
%plot(x,f*100,'-b','LineWidth',2,'DisplayName','Predicted MTDs')

line( [ T.sumObs(idx) T.sumObs(idx)] , [0 2e3] ,'Color','b','DisplayName', gcn5_lg_txt ,'LineWidth',3)
xlabel('# of MTDs per gene')

grid on ;
xlim([-1 7])
legend({'Measured MTDs' 'Predicted MTDs' gcn5_lg_txt} , 'location','ne');
set(gca,'xtick',0:100)
ylabel('# of genes')
mean( Yp>Yp(idx))*100
mean( Y>Y(idx))*100
