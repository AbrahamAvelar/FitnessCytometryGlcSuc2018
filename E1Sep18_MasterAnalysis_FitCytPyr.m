% This is the master script that joins all steps to analyze citometry data to
% obtain relative fitness in competition for this experiment
% it uses scripts from other repositories:
% ShowGatedDatam and SampleCitPlate from "Citometro"
% JitterPlot from "Utilidades_Graficas"
%%
load 20181002_MediaShort_Workspace
%MediaShort = rmfield(MediaIntegrated, 'Day')
Muts.Refs=[11 15 28 41 54 67 80 93 7 20 33 46 59 2 85 96];
mutantes = fieldnames(Muts);
for i = 1:length(mutantes)
    wells = Muts.(str2mat( mutantes(i) ));
    cepas(wells)=mutantes(i)
end


%% FIRST STEP GATE FSC SSC THIS SECTION REQUIRES THE ENTIRE DATASET, NOT
%% ONLY THE SHORT  VERSION use MediaIntegrated
Media=MediaIntegrated;
PLATO=11;
ver =1;ver2=1;gatedData=0;
x=1; y=2;
namex = 'FSC';namey = 'SSC';
samplesize=1000;
figure(333)
clf
for m=2:PLATO%:11
	%figure(200+m)
    %clf
    subplot(5,2,m-1)
    sampledData = SampleCitPlate( Media(m), 'Day', samplesize );
    ShowGatedData( sampledData, Media(m).Gate1, x, y, 'log', 'log')
    title(Media(m).name)
    xlim([10000 300000])
    ylim([5000 400000])
   % [n,c]=hist3([ (sampledData(:,x)) (sampledData(:,y)) ], [15 15]);
   % contour(c{2},c{1},n)
end
subplot(5,2,9)
% xlabel('FSC')
% ylabel('SSC')
% set(gca, 'paperposition',[1 1 10 15] )
% print('-depsc', 'GatedFSCSSC.eps')

%% SECOND STEP, GATE RFP VS CFP use MediaIntegrated
ver =1;ver2=1;gatedData=0;
x=4; y=5;
namex = 'RFP';namey = 'CFP';
samplesize=1000;
figure(334)
clf
for m=2:PLATO %11
    %figure(200+m)
    %clf
    subplot(5,2,m-1)
	sampledData = SampleCitPlate( Media(m), 'Day', samplesize );
    ShowGatedData( sampledData, Media(m).GateCFP, x, y, 'log', 'log')
    ShowGatedData( sampledData, Media(m).GateRFP, x, y, 'log', 'log')
    title(Media(m).name)
	xlim([0 300000])
    ylim([0 400000])
end
subplot(5,2,9)
xlabel('RFP')
ylabel('CFP')
print('-depsc', 'GatedRFPCFP.eps')



%% ROBUSTFIT AND PLOT OF LOG2(RFP/CFP)perGENERATIONs
Media = MediaShort;
for m = 2%2:11
    figure(m)
    clf
    subplot1(8,12,'Gap', [0.003 0.005])
    for w=1:96
        if (m==7 & w>36) 
            valid = [1 3:5 ];
        elseif (m==2 & w>58) | (m==10 & w==25)
            valid =1:4 ;
        elseif  (m==5 & find(ismember(w,[ 86 75 76 64 95 84 ])) )
            valid = 2:5;
        else
            valid = 1:5;
        end
        
        %subplot(8,12,w)
        subplot1(w)
        plot(4*valid, Media(m).RoC(valid,w), 'ok' )
        hold on
%%%%THESE TWO LINES SAVED THE RELATIVE FITNESS IT'S COMMENTED SO NOW IT
%%%%ONLY  SHOWS FIGURE
%         robfit = robustfit( valid*4, Media(m).RoC(valid,w) );%así se guardó
%         Media(m).robfit(:,w) = robfit;
        plot( [0 20], [ Media(m).robfit(1,w)  Media(m).robfit(1,w)+20*Media(m).robfit(2,w) ] )
        ylim( [-1.5 2.8] )
        xlim( [.5 20.5] )
        text(2, 2.4, strcat(cepas(w),'-Well:',num2str(w)) )
        set( gca,'xticklabel',[] )
        if find(ismember(w,Muts.Refs))
            set(gca, 'color', [ .7 .7 .7 ])
        end
    end
    
%save 20181001_MediaShort_Workspace

end
%subplot(8,12,96-11)
subplot1(90)
xlabel('Generations')
set( gca,'xtick',4:4:20,'xticklabel', 4:4:20)

subplot1(49)
ylabel('log2(RFP/CFP)')

%% PLOT REFERENCES
Muts.Refs=[11 15 28 41 54 67 80 93 7 20 33 46 59 2 85 96];
figure()
clf
mutants=fieldnames(Muts);
mapacolor=hsv;
jump=floor(64/11);
clear H
jitter=.5;

for m = 2:11 %different media
        hold on
        refs = Muts.Refs;
        
        a=JitterPlot(m-1,Media(m).robfit(2,refs), jitter, 'o',mapacolor(jump*m-1,:) );
        H(m-1)=a(1);
        
        %title( Media(m).name );
        %text(thisMutant-.5, .38, num2str(length( find(~isnan(Media(m).robfit(2,wells))) )))
        ylim( [-.05 .05] )
        xlim( [0 10.5] )
end
set(gca,'xtick',[] ,'xticklabel',MediaNames)
text( (1:10)-.1, ones(1,10)*.08, MediaNames, 'rotation', 270 )
%legend(H,MediaNames,'location','EastOutside')
ylim( [-.05 .08] )
print('-depsc', 'FitnessRefs' )

%% %% 20 subplots each mutant in 10 conditions MUTANTS & REFS
figure() 
clf
mutants=fieldnames(Muts);
mapacolor=hsv;
jump=floor(64/10)-1;
clear h
jitter=.5;
for thisMutant = 1:length(mutants)-3
    for m = 2:11 %different media
        subplot(4,5,thisMutant)
        hold on
        wells = Muts.(str2mat(mutants(thisMutant)));
        refs = Muts.Refs;
        
        h(m-1,:)=JitterPlot(m-1.2, Media(m).robfit(2,wells), jitter, 'o',mapacolor(jump*m-1,:) );
%        h(m-1)=plot( ones(1,length(wells))*m-1.2, Media(m).robfit(2,wells), '.-','color', mapacolor(jump*m-1,:) );
        plot( ones(1,length(refs))*m-.8, Media(m).robfit(2,refs), '.','color',mapacolor(jump*m-1,:) )
        title( mutants(thisMutant) )
        MediaNames(m-1)=Media(m).name;
    end
        ylim( [-.2 .2] )
    xlim([0 m])
end
legend(h(:,1),MediaNames,'location','EastOutside')

%% MEDIAN NORMALIZED 20 subplots of each mutant in 10 conditions
% Calculate Pvalue using ttest2
% Generate MedianNormFitness, a matrix with median values
figure(24) 
clf
mutants=fieldnames(Muts);
mapacolor=HSV;
jump=floor(64/10)-1;
clear h
jitter=.3;
orderMutants = [2 15 18 1 19 17 7 10 4 5 9 12 3 11 6  8 13 16 14 ] ;
con=0;
for thisMutant = orderMutants% 1:length(mutants)-3%
    con=con+1;
    for m = 2:11 %different media
        subplot(4,5,con)
        
        %figure(thisMutant)
        hold on
        wells = Muts.(str2mat(mutants(thisMutant)));
        refs = Muts.Refs;
        
        h(m-1,:)=JitterPlot(m-1.15, Media(m).robfit(2,wells)-median(Media(m).robfit(2,refs)), jitter, 'o',mapacolor(jump*m-1,:) );
        %JitterPlot( m-.8, Media(m).robfit(2,refs)-median(Media(m).robfit(2,refs)), jitter, '.', [.5 .5 .5])
        %plot( ones(1,length(refs))*(m-.8), Media(m).robfit(2,refs)-median(Media(m).robfit(2,refs)), 'o', 'color',[.5 .5 .5],'markersize',5)
        JitterPlot(m-.8, Media(m).robfit(2,refs)-median(Media(m).robfit(2,refs)), jitter-.1, '.',mapacolor(jump*m-1,:) );
        MedianNormFitness(thisMutant,m) = median(  Media(m).robfit(2,wells)-median(Media(m).robfit(2,refs)) );
        
        [ttestsH(thisMutant,m) ttestsP(thisMutant,m)] = ttest2( Media(m).robfit(2,refs)-median(Media(m).robfit(2,refs)),  Media(m).robfit(2,wells)-median(Media(m).robfit(2,refs)) );
        
        title( mutants(thisMutant) )
        MediaNames(m-1)=Media(m).name;
        set( gca, 'xtick', [] )
    end
        ylim( [-.12 .075] )
    xlim( [0 m] )
end
legend( h(:,1),MediaNames)%,'location','EastOutside' )


%% clustergram MedianNormFitness
%complete, with all data
metodo= 'euclidean';%'spearman';
clustergram(MedianNormFitness(orderMutants,2:end),'RowLabels',mutants(orderMutants),'ColumnLabels',MediaNames,'Standardize',3,'RowPDist',metodo,'ColumnPDist',metodo)

%% Clustergram significativos
metodo= 'euclidean';%'spearman';
clear significant significantMNF
significant = ttestsP<0.01 & ttestsP>0 %this selects acceptable pvalues
SignificantMNF(19,11)=0
SignificantMNF(significant) = MedianNormFitness(significant)
clustergram(SignificantMNF(orderMutants,2:end),'RowLabels',mutants(orderMutants),'ColumnLabels',MediaNames,'Standardize',3,'RowPDist',metodo,'ColumnPDist',metodo)


%%
