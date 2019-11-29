function MH_Event_Frequency_as_a_function_of_MHlen_and_distance( FILENAME_sign_count_tsv )
%% MH_Event_Frequency_as_a_function_of_MHlen_and_distance( FILENAME_sign_count_tsv )
%   calculate and plot duplication & collapse frequencies as a function of
%    the length of the MH repeats, and as a function of the distance between repeats
%
%    Input : the output of catch_signatures.awk
%       I	4048	4051	4060	4063	0	5050
%       I	4054	4057	4066	4069	3	1
%       I	4048	4051	4066	4069	0	2
%      chr   st      end     st      en   iCount  dCount
%    iCount = Insertions      dCount = Deletions / Collapse
%
%    output is a pair of tables with the event frequencies calculated,
%    grouped by MHlen or by inter-MH distance
% 
% October 2019, LBC
%% load data
if ~exist('FILENAME_sign_count_tsv' , 'var')
    FILENAME_sign_count_tsv = '~/Downloads/10k.sign.count.tsv' ; 
end

T = readtable( FILENAME_sign_count_tsv , 'FileType','text','Format','%s%d%d%d%d%d%d');
T.Properties.VariableNames = {'chr' 's1' 'e1' 's2' 'e2' 'DupCounts' 'CollapseCounts'};
T.HasDup = T.DupCounts > 0 ; 
T.HasColl = T.CollapseCounts > 0 ; 

T.MHlen = T.e1-T.s1+1 ;
T.TotalLen = T.e2 - T.s1 ; 
T.NTBetweenRepeats = T.s2 - T.e1 ; 
T.chr = categorical(T.chr) ; 
T.s1 = uint32(T.s1) ; T.e1 = uint32(T.e1) ; T.s2 = uint32(T.s2) ; T.e2 = uint32(T.e2) ; 

%% group by MH-length ,and by distance-between-repeats, for calculating statistics
G_MHlen = grpstats( T , 'MHlen' , 'mean'  , 'DataVars' , {'HasDup' 'HasColl'} );
G_NTBetweenRepeats = grpstats( T , 'NTBetweenRepeats' , 'mean'  , 'DataVars' , {'HasDup' 'HasColl'} );
G_TotalLen = grpstats( T , 'TotalLen' , 'mean'  , 'DataVars' , {'HasDup' 'HasColl'} );

%% Plot figure ; 
fh = figure('units','centimeters','position',[5 5 15 15 ]) ;

subplot(2,2,1); hold on; 
plot( G_MHlen.MHlen , G_MHlen.mean_HasDup ,'.-','LineWidth',2,'DisplayName','Dup')
plot( G_MHlen.MHlen , G_MHlen.mean_HasColl ,'.-','LineWidth',2,'DisplayName','Clps')
xlabel('MHlen (nt)')
ylabel('Frequency')

subplot(2,2,2); hold on; 
plot( G_NTBetweenRepeats.NTBetweenRepeats , G_NTBetweenRepeats.mean_HasDup ,'.-','LineWidth',2,'DisplayName','Dup')
plot( G_NTBetweenRepeats.NTBetweenRepeats , G_NTBetweenRepeats.mean_HasColl ,'.-','LineWidth',2,'DisplayName','Clps')
xlabel('NTBetweenRepeats (nt)')
ylabel('Frequency')

subplot(2,2,3); hold on; 
plot( G_TotalLen.TotalLen , G_TotalLen.mean_HasDup ,'.-','LineWidth',2,'DisplayName','Dup')
plot( G_TotalLen.TotalLen , G_TotalLen.mean_HasColl ,'.-','LineWidth',2,'DisplayName','Clps')
xlabel('TotalLen (nt)')
ylabel('Frequency')