T = readtable("/Users/lcarey/SynologyDrive/Projects/2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU/Sarah/MH_project/ProcessedData/loc_obs_pre1k.txt");
T.loc = T.loc ./ 1000  ;
locs_of_chr_breaks = find(diff(T.loc)<0) ; 
T.chr = NaN(height(T),1);
T.chr(1:locs_of_chr_breaks(1))=1;
T.chr(locs_of_chr_breaks(1)+1:locs_of_chr_breaks(2))=2;
T.chr(locs_of_chr_breaks(2)+1:end)=3;

T = T(1:locs_of_chr_breaks(1),:);
%T = T(locs_of_chr_breaks(1)+1:locs_of_chr_breaks(2) , : ); % chrII 
%T = T(locs_of_chr_breaks(2)+1:end , : ); % chrIII 

%T = T( T.loc > 1500 & T.loc < 3500  ,:);

clrs = cbrewer('qual','Set1',5);
fh = figure('units','centimeters','position',[5 5  70 20]) ;
t = tiledlayout(4,1,'TileSpacing','none','Padding','none')
nexttile
plot( T.loc , T.sumMHR ,'Color',clrs(1,:) )
ylabel('# of MHPs')
axis tight; 
nexttile

plot( T.loc , T.sumObs ,'Color',clrs(2,:) )
ylabel('# of observed MTDs')
axis tight; 
ylim([-0.5  6.5])
nexttile

plot( T.loc , T.sumPre ,'Color',clrs(3,:) )
ylabel('# predicted MTDs')
axis tight; 
ylim([-0.5  6.5])
nexttile

plot( T.loc , T.sumFloat ,'Color',clrs(4,:) )
xlabel('Position on chromosome (kb)')
ylabel('MTD prediction score')
axis tight; 

title(t,'MTDs across chromosome I')

print('-dpng','~/Downloads/MTD_observed_predict_across_chromosome','-r300');
close ; 