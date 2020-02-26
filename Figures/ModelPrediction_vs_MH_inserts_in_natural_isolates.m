% Genes with MTDs in natural isolates score higher in the model

%% load natural strain data
FIGDIR  = '~/Nutstore Files/Microhomology shared folder/Figures/Fig4/' ;
FIGNAME = [ FIGDIR 'ModelPrediction_vs_MH_inserts_in_natural_isolates' ] ; 

FN = [ FIGDIR '57 nature isolates_sequencing summary.xlsx' ] ;

T = readtable( FN , 'Sheet' , 2);
T = T( : , 1:13) ; 
T.Properties.VariableNames{6} = 'REF_N' ; 
T.Properties.VariableNames{7} = 'ALT_N' ; 
T = T(2:end,:);
T.REF_N = str2double(T.REF_N) ; 
T.ALT_N = str2double(T.ALT_N) ; 
T.HasMH = ~cellfun(@isempty,T.MH);

T.REF_len = cellfun(@length,T.REF);
T.ALT_len = cellfun(@length,T.ALT);
T.IsInsertion = T.ALT_len > T.REF_len ; 

T.Location = categorical(T.Location);

T = T( T.IsInsertion , :); 
%% load model predictions
DATADIR = '~/SynologyDrive/Projects/2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU/Sarah/MH_project/Manuscript-todo/processeddata/' ;
G = readtable([DATADIR 'MHRSumPreinGene.txt' ],'TreatAsEmpty','.');
G.Properties.VariableNames = { 'chr'	'start1'	'end2' 'dunno'	'gene'	'systematic_name'	'sumObs'	'sumPre'	'sumFloat' 'dunno2'}; 
G.GeneLengthKB = abs( G.end2 - G.start1 ) ./ 1000 ; 


D = readtable([DATADIR 'MHRinGene.txt' ],'TreatAsEmpty','.');
D.Properties.VariableNames = {'chr' 'start1' 'end1' 'type' 'systematic_name' 'N_MHPs'} ;
D.systematic_name = regexprep( regexprep( D.systematic_name , 'ID=','') , ';Name.*' , '') ;


G = innerjoin(D , G  , 'key' , 'systematic_name'); 
for I  =  1:height(G)
    if strcmp(G.gene{I},'NA')
        G.gene{I} = G.systematic_name{I};
    end
end
%% join the two
Q = outerjoin( T , G(:,{'N_MHPs' 'sumFloat' 'sumObs' 'gene' 'systematic_name' 'GeneLengthKB'}) , 'LeftKey','Gene','RightKey','gene');
Q.sumFloat_per_kb = Q.sumFloat ./ Q.GeneLengthKB ; 
Q.sumObs_per_kb = Q.sumObs ./ Q.GeneLengthKB ; 


idx_for_MH_insertion = Q.HasMH & strcmp(Q.REF,'.') & Q.IsInsertion ;


[p,~,~] = ranksum( Q.sumObs(~idx_for_MH_insertion & ~isnan(Q.sumObs) ) , Q.sumObs(idx_for_MH_insertion & ~isnan(Q.sumObs) ))
[p,~,~] = ranksum( Q.sumFloat(~idx_for_MH_insertion & ~isnan(Q.sumFloat) ) , Q.sumFloat(idx_for_MH_insertion & ~isnan(Q.sumFloat) ))

[pN1,~,~] = ranksum( Q.sumObs_per_kb(~idx_for_MH_insertion & ~isnan(Q.sumObs_per_kb) ) , Q.sumObs_per_kb(idx_for_MH_insertion & ~isnan(Q.sumObs_per_kb) ))
[pN2,~,~] = ranksum( Q.sumFloat_per_kb(~idx_for_MH_insertion & ~isnan(Q.sumFloat_per_kb) ) , Q.sumFloat_per_kb(idx_for_MH_insertion & ~isnan(Q.sumFloat_per_kb) ))

%%
fh = figure('units','centimeters','position',[5 5  6 4]) 
clrs = get(gca,'ColorOrder');
hold on ;
Y = Q.sumFloat ./ Q.Gene ; 
Y(Y>250)=250; 
%histogram(Y(~idx_for_MH_insertion), 0:25:max(Y) ,'Normalization','Probability' )
%histogram( Y(idx_for_MH_insertion) , 0:25:max(Y) , 'Normalization','Probability' )


[f,x] = ksdensity( Y(~idx_for_MH_insertion)  , 0:1:max(Y) );
plot(x,f,'-','LineWidth',3,'DisplayName','no MH insert','Color',[.4 .4 .4])


[f,x] = ksdensity( Y(idx_for_MH_insertion)  , 0:1:max(Y)  );
plot(x,f,'-','LineWidth',3,'DisplayName','MH insert found' , 'Color',clrs(1,:))

ylabel('Fraction of genes')
xlabel('MTD model score')
%legend('location','ne','box','off')
set(gca,'ytick',[0 max(get(gca,'ytick'))])
xlim([0 max(Y)])
print( '-dpng' , [ FIGNAME '_prediction_score'] , '-r300')

[p,h,stats] = ranksum( Y(~idx_for_MH_insertion & ~isnan(Y) ) , Y(idx_for_MH_insertion & ~isnan(Y) ))

%%
fh = figure('units','centimeters','position',[5 5  6 4]) ;
hold on ;

Y = Q.sumObs ; 
Y(Y>7)=7; 
%histogram(Y(~idx_for_MH_insertion), 0:25:max(Y) ,'Normalization','Probability' )
%histogram( Y(idx_for_MH_insertion) , 0:25:max(Y) , 'Normalization','Probability' )

xl = -0.5:8 ; 
histogram( Y(~idx_for_MH_insertion & ~isnan(Y) )  , xl,'Normalization','Probability','FaceColor',[.2 .2 .2]);
histogram( Y(idx_for_MH_insertion & ~isnan(Y) )  , xl ,'Normalization','Probability' , 'FaceColor',clrs(1,:) );

ylabel('Fraction of genes')
xlabel('# of MTDs observed')
%legend( {'no MH insert' 'MH insert found'} , 'location','best','box','off')
set(gca,'yscale','lin')

[p,~,stats] = ranksum( Y(~idx_for_MH_insertion & ~isnan(Y) ) , Y(idx_for_MH_insertion & ~isnan(Y) ))

print( '-dpng' , [ FIGNAME '_MTDs_observed'] , '-r300')



