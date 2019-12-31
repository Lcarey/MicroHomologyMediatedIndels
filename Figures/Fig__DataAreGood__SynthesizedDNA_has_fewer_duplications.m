%% load data
% takes as input, result of 
% Develop/MicroHomologyMediatedIndels/Analysis/Lucas/MH_Event_Frequency_as_a_function_of_MHlen_and_distance.m
%T = readtable('~/Downloads/alltsvs.txt','TreatAsEmpty','-');
%T.Properties.VariableNames = {'lib' 'chr' 's1' 'e1' 's2' 'e2' 'nreads ' 
idx = G.MHlen<=4 & regexpcmp(G.chr,'1') & G.median_ReadDepth > 20000;
Y=G.mean_DupFreq;
Y(Y>5)=5;
ym = max(Y(idx));
figure; hold on ;
histogram( Y(idx & G.ChemSyn==0),0:0.05:ym)
histogram( Y(idx & G.ChemSyn==1),0:0.05:ym)
 [~,p] = kstest2(  Y(idx & G.ChemSyn==0) ,  Y(idx & G.ChemSyn==1));
 title( sprintf('K.S. test p=%0.05f' , p )  )
xlabel('Mean Duplication Frequency')
ylabel('# of experiments')
legend({'native seq' 'chemically syn'})
%%

TotalDupEventCounts_Native = sum(T.DupCounts(strcmp(T.chr,'ssp1.dup.2')))   % 12,555
TotalDupEventCounts_SYN = sum(T.DupCounts(strcmp(T.chr,'ssp1.short.1.PCR')))   % 2

% $ samtools view -c lib_14.bam ssp1.dup.2    ->  43420306
% $ samtools view -c lib_1.bam ssp1.short.1   ->     73303
%12,555 / 43420306 = 0.000289150

% p = 1.2478e-09

% $ samtools view -c lib_2.bam ssp1.short.1.PCR  ->  2405009

% X = [#_reads_Native #_dup_events_native
%     #_reads_synth  #_dup_events_synth ] 
 
X = [ 43420306 12555 ;...
    2405009 2] ;

[~,p,festats] = fishertest( X )

% ../sequencing_raw_data/E4/lib_2.fa:1:>ssp1.short.1.PCR

 p =  4.5889e-289
 
 
 
%% SSP1 FLASH results
%$ cut -f 2 lib_15.readsAndPairs.tab |grep GAGAACGTGAAGGAAGTTCGTTAACTCACTCATGGACTTTTCAACCTGGTAAGCATAACCAGCGTCTTTATTCTGATAATTTTCAAGAGGCTCAGCGCCAGTGGAAGCGCCTGCAAGAATGGGGCGAGGTGAAGGAAACAAAA |wc -l
% 237712 events
fid=fopen('~/Downloads/ssp1_89_dup');
C = textscan(fid,'%s');
fclose(fid);
C=C{:};
C = C( ~regexpcmp( C , 'AGAGGTAGACGATAGCGAAGTTCCTCCTTCTGTCTTTCCTGAATATCCCGTCCACAAGGCCATCCAGAGAGGTAGACGATAGCGAAGTTCCTCCTTCTGTCTTTCCTGAATATCCCGTCCACAAGGCCATCCAGAAAACGTCCGATTCATTTCGTAAACGGAACTACAGCGCGGGAGATTATGTA'));


m=multialign(C) ; showalignment(m)