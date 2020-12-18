% Plot proliferation scores for TMRE sorted cells
% companion to 
%   /Volumes/CareyLab/Projects/2019__ProliferationCorrelatedExpression/ProcessedData/
%      TMRE sort doubling time figures/TMRE doubling time figures.m

%% load data
PROJDIR = '/Volumes/CareyLab/Projects/2019__ProliferationCorrelatedExpression/' ;
DATADIR = [PROJDIR '20200805 new result for proliferation/data/TMRE vs Proliferation rate/'] ;
DATAFILE = [DATADIR 'fib_esc_gsva.tsv'] ;

PROJDIR = '/Users/lcarey/Develop/Proliferation_Signature_Public/'
DATAFILE = 'Figure6/data/fib_esc_gsva__ProliferationSignatureScores_TMRE_sorted_cells.tsv'
T = readtable([PROJDIR DATAFILE],'filetype','text');

T.Properties.VariableNames = {'key' 'ProliferationScore' 'CellType' 'TMRE'};
T.CellType = categorical(T.CellType) ;
T.TMRE = categorical(T.TMRE) ;
T.replicate = categorical(regexprep(T.key,'.*_',''));
T.key = [];

T.log2_High_over_Low = NaN(height(T),1);
for I = 1:height(T)
    if (T.TMRE(I)=='high')
        idx = T.TMRE=='low' & T.replicate == T.replicate(I) & T.CellType == T.CellType(I);
        T.log2_High_over_Low(I) = log2(T.ProliferationScore(I) / T.ProliferationScore(idx));
    end
    
end
T

clrs = cbrewer('seq' , 'Reds' , 6) ;
% clrs = clrs([1 3 6],:); % low, med, high
clrs = clrs([6 1],:); %low,high

%% p-values
p_5_of_6 = myBinomTest(5,6,0.5,'one')
p_3_of_3 = myBinomTest(3,3,0.5,'one')

%p = myBinomTest(2,3,0.5,'one')
%p = myBinomTest(3,3,0.5,'one')

[~,p] = ttest(T.log2_High_over_Low(T.CellType=='FIB'))
[~,p] = ttest(T.log2_High_over_Low(T.CellType=='ESC'))

%% plot figure

lw = 1; 
v1 = [T.ProliferationScore(4),T.ProliferationScore(1)];
v2 = [T.ProliferationScore(5),T.ProliferationScore(2)];
v3 = [T.ProliferationScore(6),T.ProliferationScore(3)];


figure;
tiledlayout(2,2)
nexttile; hold on ;
gscatter([2 2 2 1 1 1],T.ProliferationScore(1:6),[1 1 1 3 3 3 ],clrs)
plot( v1 , '-pk','LineWidth',lw);
plot( v2, '-sk','LineWidth',lw);
plot( v3 , '-vk','LineWidth',lw);
xlim([0.5 2.5])
set(gca,'xtick',1:2)
set(gca,'xticklabels',{'low' 'high'})
legend('off')
ylabel('Proliferation Score')
%title('Fibroblasts','FontWeight','Normal')

v1 = [T.ProliferationScore(10),T.ProliferationScore(7)];
v2 = [T.ProliferationScore(11),T.ProliferationScore(8)];
v3 = [T.ProliferationScore(12),T.ProliferationScore(9)];
nexttile; hold on ;
gscatter([2 2 2 1 1 1],T.ProliferationScore(7:12),[1 1 1 3 3 3 ],clrs)
plot( v1 , '-pk','LineWidth',lw);
plot( v2, '-sk','LineWidth',lw);
plot( v3 , '-vk','LineWidth',lw/3);
xlim([0.5 2.5])
set(gca,'xtick',1:2)
set(gca,'xticklabels',{'low' 'high'})
legend('off')
ylabel('Proliferation Score')
%title('ESCs','FontWeight','Normal')


%% with inverse Y axis
v1 = [T.ProliferationScore(4),T.ProliferationScore(1)];
v2 = [T.ProliferationScore(5),T.ProliferationScore(2)];
v3 = [T.ProliferationScore(6),T.ProliferationScore(3)];
figure;
tiledlayout(2,2)
nexttile; hold on ;
gscatter([2 2 2 1 1 1],T.ProliferationScore(1:6).^-1,[1 1 1 3 3 3 ],clrs)
plot( v1.^-1 , '-pk','LineWidth',lw);
plot( v2.^-1, '-sk','LineWidth',lw);
plot( v3.^-1 , '-vk','LineWidth',lw);
xlim([0.5 2.5])
set(gca,'xtick',1:2)
set(gca,'xticklabels',{'low' 'high'})
legend('off')
ylabel('Proliferation Score (^{-1})')
title('Fibroblasts','FontWeight','Normal')

v1 = [T.ProliferationScore(10),T.ProliferationScore(7)];
v2 = [T.ProliferationScore(11),T.ProliferationScore(8)];
v3 = [T.ProliferationScore(12),T.ProliferationScore(9)];
nexttile; hold on ;
gscatter([2 2 2 1 1 1],T.ProliferationScore(7:12).^-1,[1 1 1 3 3 3 ],clrs)
plot( v1.^-1 , '-pk','LineWidth',lw);
plot( v2.^-1, '-sk','LineWidth',lw);
plot( v3.^-1 , '-vk','LineWidth',lw/3);
xlim([0.5 2.5])
set(gca,'xtick',1:2)
set(gca,'xticklabels',{'low' 'high'})
legend('off')
ylabel('Proliferation Score (^{-1})')
title('ESCs','FontWeight','Normal')
%% 

[p,t,stats] = anovan( T.ProliferationScore ,  {T.TMRE , T.CellType , T.replicate } ...
    , 'VarNames' ,{'TMRE' 'CellType' 'replicate'} , 'model','linear') 

%% 
ESC = T(T.CellType=='ESC',:)
FIB = T(T.CellType=='FIB',:)
figure;
tiledlayout(2,2)
nexttile;
hold on ;
plot(FIB.ProliferationScore(FIB.TMRE=='high'), FIB.ProliferationScore(FIB.TMRE=='low'),'ok');
line(xlim,xlim)
xlabel('High TMRE')
ylabel('Low TMRE')
title('FIB')

nexttile;
hold on ;
plot(ESC.ProliferationScore(ESC.TMRE=='high'), ESC.ProliferationScore(ESC.TMRE=='low'),'ok');
line(xlim,xlim)
xlabel('High TMRE')
ylabel('Low TMRE')
title('ESC')

%%

figure('Position',[500 500 125 150]); 
hold on ; 

for I = 1:3
    h = errorbar( I , G.mean_doubling_time(I) , G.std_doubling_time(I) ,'sk' ,'MarkerFaceColor',clrs(I,:) ... ,
     ,  'MarkerSize',10  );
end
set(gca,'xtick',1:3)
set(gca,'xticklabel',{'Low' 'Med' 'High'})
xlabel( 'TMRE' )
xlim([0.4 3.5])
ylabel('Doubling Time (hours)')

sprintf('ANOVA p=%0.04f' , p(1) ) 
