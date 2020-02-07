% Fig__Fraction_of_MHPs_with_MTD__downsample_amplicon_library
% if we reduce the coverage, how many MHPairs would have an MTD? 
T = readtable('~/Downloads/PombeAmpliconSeq_E2_alltsvs.txt','FileType','text','Format','%s%s%d%d%d%d%d%d%d%d%d','TreatAsEmpty','-');
T = T(ismember(T.Var2, {'ssp1' 'ssp2' 'SPCC1235.01' }) & ismember(T.Var1, {'lib_1' 'lib_2' })  ,:) ;

DIR = '~/CareyLab/Projects/2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU/DataFromCluster/';
G = readtable( [ DIR '10k_rm.sign.count.tsv'] ,'FileType','text');
GENOME_NDups = G.Var6 ; 
clear 'G' 
%% for testing
R = table();
R.xl = ([1:0.05:100 101:10:5000 6000:1e3:1e4 1e5] )' ; 


Q = T(ismember(T.Var2, {'ssp1' }) & ismember(T.Var1, {'lib_1'  })  ,:) ;
R.ndupreads_ssp1_lib1 = arrayfun( @(X)mean(round(Q.Var8 ./ X) > 0)*100, R.xl ) ; 
R.nreads_ssp1_lib1    = arrayfun( @(X)mean(round(Q.Var7 ./ X)), R.xl ) ; 

Q = T(ismember(T.Var2, {'ssp1' }) & ismember(T.Var1, {'lib_2'  })  ,:) ;
R.ndupreads_ssp1_lib2 = arrayfun( @(X)mean(round(Q.Var8 ./ X) > 0)*100, R.xl ) ; 
R.nreads_ssp1_lib2    = arrayfun( @(X)mean(round(Q.Var7 ./ X)), R.xl ) ; 

Q = T(ismember(T.Var2, { 'SPCC1235.01' }) & ismember(T.Var1, {'lib_1'  })  ,:) ;
R.ndupreads_SPCC1235_lib1 = arrayfun( @(X)mean(round(Q.Var8 ./ X) > 0)*100, R.xl ) ; 
R.nreads_SPCC1235_lib1    = arrayfun( @(X)mean(round(Q.Var7 ./ X)), R.xl ) ; 

Q = T(ismember(T.Var2, { 'SPCC1235.01' }) & ismember(T.Var1, {'lib_2'  })  ,:) ;
R.ndupreads_SPCC1235_lib2 = arrayfun( @(X)mean(round(Q.Var8 ./ X) > 0)*100, R.xl ) ; 
R.nreads_SPCC1235_lib2    = arrayfun( @(X)mean(round(Q.Var7 ./ X)), R.xl ) ; 

R.ndupreads_10k = arrayfun( @(X)mean(round(GENOME_NDups ./ X) > 0)*100, R.xl ) ; 
R.nreads_10k    = arrayfun( @(X)mean(round(1e4 ./ X)), R.xl ) ; 

%%
fh = figure('units','centimeters','position',[5 5 8 8]);
hold on ; 
plot( R.nreads_SPCC1235_lib1 , R.ndupreads_SPCC1235_lib1 ,'-','LineWidth',3,'DisplayName' , 'SPCC1235 rep 1')
plot( R.nreads_SPCC1235_lib2 , R.ndupreads_SPCC1235_lib2 ,'-','LineWidth',3,'DisplayName' , 'SPCC1235 rep 2')
plot( R.nreads_ssp1_lib1 , R.ndupreads_ssp1_lib1 ,'-','LineWidth',3,'DisplayName' , 'ssp1 rep 1')
plot( R.nreads_ssp1_lib2 , R.ndupreads_ssp1_lib2 ,'-','LineWidth',3,'DisplayName' , 'ssp1 rep 2')
plot( R.nreads_10k , R.ndupreads_10k ,'-k','LineWidth',3,'DisplayName' , '10k genome')
set(gca,'xscale','log')
xlabel('Coverage (avg. # of reads)')
ylabel('% of MHPs with an MTD')
legend('location','nw','box','off')
axis tight; 
set(gca,'xtick',logspace(1,7,7))
xlim([99 max(xlim)])
print('-dpng','~/Downloads/Fig__Fraction_of_MHPs_with_MTD__downsample_amplicon_library_lin','-r300');
set(gca,'yscale','log')
legend('off')
ylim([1e-3 max(ylim)])
set(gca,'ytick',logspace(-5,2,8))
print('-dpng','~/Downloads/Fig__Fraction_of_MHPs_with_MTD__downsample_amplicon_library_log','-r300');
close  ; 

%%
X  = R.nreads_ssp1_lib1  ; 
Y  = R.ndupreads_ssp1_lib1 ;
Y2 = Y(X>1.5*10^5);
X2 = X(X>1.5*10^5);

%X  = R.nreads_SPCC1235_lib1 ; 
%Y  = R.ndupreads_SPCC1235_lib1 ;
[xData, yData] = prepareCurveData( X2, Y2 );
ft = fittype( 'power2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
[fitresult, gof] = fit( xData, yData, ft, opts );

xl = logspace(1,9,1000) ;
Y =  feval(fitresult,xl) ; 

fh = figure('units','centimeters','position',[5 5 8 8]);
yyaxis left ; 

plot( xl , Y ,'-k','LineWidth',3,'DisplayName' , 'prediction : 10k genome')

xlim([50 1.2e8])
ylim([0 100])
ylabel('% of MHPs with an MTD')
ax.ytick = 0:10:100 ;  
set(gca,'YColor','k') ; 


yyaxis right
plot( xl , Y ,'-k','LineWidth',3,'DisplayName' , 'prediction : 10k genome')
xlim([50 1.2e8])
ylim([0 100])
set(gca,'YColor','k') ; 
% YTickLabels at 10^6 & 10^7 : (25*1e6) .* [1:-0.1:0] , assuming 25 million MTDs
ax.ytick = [0:20:100] ; 
yticklabels( {'0' '5^6' '10^6' '15^6' '20^6'  '25^6'}  ) ; 
ylabel('# of MTDs detected')
xlabel('Coverage (# of reads)')
print('-dpng','~/Downloads/Fig__Fraction_of_MHPs_with_MTD__downsample_amplicon_library_pred','-r300');
close  ; 