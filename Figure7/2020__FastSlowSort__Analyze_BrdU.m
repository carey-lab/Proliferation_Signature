%% load data
WORKDIR  = './' ;
BRDU = [ WORKDIR 'data/2020__FastSlowSort__Brdu_data.xlsx' ] ;
VIAB = [ WORKDIR 'data/2020__FastSlowSort__viability_data.xlsx' ] ;

% look at BrdU
E = readtable( BRDU , 'Sheet' , 'BrdU_ESC');
F = readtable( BRDU , 'Sheet' , 'BrdU_Fibro');
E = E( ~isnan(E.Replicate1),:);
F = F( ~isnan(F.Replicate1),:);
E.CellsInS = mean([E.Replicate1,E.Replicate2],2) ;
F.CellsInS = mean([F.Replicate1,F.Replicate2,F.Replicate3],2) ;

E.Type = regexprep(E.Type,' cells','');
F.Type = regexprep(F.Type,' cells','');

F.Speed = categorical( regexprep(F.Type,' .*','') );
E.Speed = categorical( regexprep(E.Type,' .*','') );
F.Drug = categorical( regexprep(F.Type,'.* ','') );
E.Drug = categorical( regexprep(E.Type,'.* ','') );

E.CellType = repmat('ESC',height(E),1)
F.CellType = repmat('FIB',height(F),1)

% relative to DMSO on that day
%F.Rel2DMSO = NaN(height(F),1);
for I = 1:height(F)
 %   F.Rel2DMSO(I) = log2( ...
  %      F.CellsInS(I) / F.CellsInS( F.Day==F.Day(I) & F.Drug=='DMSO' & F.Speed==F.Speed(I)) ...
  %      );
    F.Rel2DMSO_1(I) = log2( ...
        F.Replicate1(I) / F.Replicate1( F.Day==F.Day(I) & F.Drug=='DMSO' & F.Speed==F.Speed(I)) ...
        );
    F.Rel2DMSO_2(I) = log2( ...
        F.Replicate2(I) / F.Replicate2( F.Day==F.Day(I) & F.Drug=='DMSO' & F.Speed==F.Speed(I)) ...
        );
    F.Rel2DMSO_3(I) = log2( ...
        F.Replicate3(I) / F.Replicate3( F.Day==F.Day(I) & F.Drug=='DMSO' & F.Speed==F.Speed(I)) ...
        );            
end

F.Type=[];
%F.Replicate1=[];
%F.Replicate2=[];


%E.Rel2DMSO = NaN(height(E),1);
for I = 1:height(E)
  %  E.Rel2DMSO(I) = log2( ...
  %      E.CellsInS(I) / E.CellsInS( E.Day==E.Day(I) & E.Drug=='DMSO' & E.Speed==E.Speed(I)) ...
  %      );
    E.Rel2DMSO_1(I) = log2( ...
        E.Replicate1(I) / E.Replicate1( E.Day==E.Day(I) & E.Drug=='DMSO' & E.Speed==E.Speed(I)) ...
        );
    E.Rel2DMSO_2(I) = log2( ...
        E.Replicate2(I) / E.Replicate2( E.Day==E.Day(I) & E.Drug=='DMSO' & E.Speed==E.Speed(I)) ...
        );
            
end

E.Type=[];
%E.Replicate1=[];
%E.Replicate2=[];

%%
%F = F(F.Day<3,:);
%%
clrs = 'gr';
figure; 
tiledlayout(2,5)


nexttile; hold on ;
idx=F.Drug=='Antimycin';
gscatter(F.Day(idx),F.Rel2DMSO_1(idx),F.Speed(idx),clrs,'s',5);
gscatter(F.Day(idx),F.Rel2DMSO_2(idx),F.Speed(idx),clrs,'^',5);
gscatter(F.Day(idx),F.Rel2DMSO_3(idx),F.Speed(idx),clrs,'p',5);
%gscatter(F.Day(idx),F.Rel2DMSO(idx),F.Speed(idx),clrs,'+',10);
title('FIB Antimycin')
ylabel('log2(%S / DMSO)')

nexttile; hold on ;
idx=F.Drug=='Oligomycin';
gscatter(F.Day(idx),F.Rel2DMSO_1(idx),F.Speed(idx),clrs,'s',5);
gscatter(F.Day(idx),F.Rel2DMSO_2(idx),F.Speed(idx),clrs,'^',5);
gscatter(F.Day(idx),F.Rel2DMSO_3(idx),F.Speed(idx),clrs,'p',5);
%gscatter(F.Day(idx),F.Rel2DMSO(idx),F.Speed(idx),clrs,'+',10);
title('FIB Oligomycin')
ylabel('log2(%S / DMSO)')


nexttile; hold on ;
idx=F.Drug=='DMSO';
gscatter(F.Day(idx),100*F.Replicate1(idx),F.Speed(idx),clrs,'s',5);
gscatter(F.Day(idx),100*F.Replicate2(idx),F.Speed(idx),clrs,'^',5);
gscatter(F.Day(idx),100*F.Replicate3(idx),F.Speed(idx),clrs,'p',5);
gscatter(F.Day(idx),100*F.CellsInS(idx),F.Speed(idx),clrs,'+',10);
legend('off')
title('FIB DMSO')
ylabel('% of cells in S')


nexttile; hold on ;
idx=F.Drug=='Antimycin';
gscatter(vertcat(F.Day(idx),F.Day(idx),F.Day(idx)),100*vertcat(F.Replicate1(idx),F.Replicate2(idx),F.Replicate3(idx)),vertcat(F.Speed(idx),F.Speed(idx),F.Speed(idx)),clrs,'o',5);
gscatter(F.Day(idx),100*F.CellsInS(idx),F.Speed(idx),clrs,[],20);
legend('off')
title('FIB Antimycin')
ylabel('% of cells in S')

nexttile; hold on ;
idx=F.Drug=='Oligomycin';
gscatter(vertcat(F.Day(idx),F.Day(idx),F.Day(idx)),100*vertcat(F.Replicate1(idx),F.Replicate2(idx),F.Replicate3(idx)),vertcat(F.Speed(idx),F.Speed(idx),F.Speed(idx)),clrs,'o',5);
gscatter(F.Day(idx),100*F.CellsInS(idx),F.Speed(idx),clrs,[],20);
legend('off')
title('FIB Oligomycin')
ylabel('% of cells in S')

nexttile; hold on ;
idx=E.Drug=='Antimycin';
gscatter(vertcat(E.Day(idx),E.Day(idx)),vertcat(E.Rel2DMSO_1(idx),E.Rel2DMSO_2(idx)),vertcat(E.Speed(idx),E.Speed(idx)),clrs,'o',5);
%gscatter(E.Day(idx),E.Rel2DMSO(idx),E.Speed(idx),clrs,[],20);
title('ESC Antimycin')
ylabel('log2(%S / DMSO)')
legend('off')
xlabel('Day')

nexttile; hold on ;
idx=E.Drug=='Oligomycin';
gscatter(vertcat(E.Day(idx),E.Day(idx)),vertcat(E.Rel2DMSO_1(idx),E.Rel2DMSO_2(idx)),vertcat(E.Speed(idx),E.Speed(idx)),clrs,'o',5);
%gscatter(E.Day(idx),E.Rel2DMSO(idx),E.Speed(idx),clrs,[],20);
title('ESC Oligomycin')
ylabel('log2(%S / DMSO)')
legend('off')
xlabel('Day')

nexttile; hold on ;
idx=E.Drug=='DMSO';
gscatter(E.Day(idx),100*E.Replicate1(idx),E.Speed(idx),clrs,'s',5);
gscatter(E.Day(idx),100*E.Replicate2(idx),E.Speed(idx),clrs,'^',5);
gscatter(E.Day(idx),100*E.CellsInS(idx),E.Speed(idx),clrs,[],20);
title('ESC DMSO')
ylabel('% of cells in S')
legend('off')
xlabel('Day')

nexttile; hold on ;
idx=E.Drug=='Antimycin';
gscatter(vertcat(E.Day(idx),E.Day(idx)),100*vertcat(E.Replicate1(idx),E.Replicate2(idx)),vertcat(E.Speed(idx),E.Speed(idx)),clrs,'o',5);
gscatter(E.Day(idx),100*E.CellsInS(idx),E.Speed(idx),clrs,[],20);
title('ESC Antimycin')
ylabel('% of cells in S')
legend('off')
xlabel('Day')

nexttile; hold on ;
idx=E.Drug=='Oligomycin';
gscatter(vertcat(E.Day(idx),E.Day(idx)),100*vertcat(E.Replicate1(idx),E.Replicate2(idx)),vertcat(E.Speed(idx),E.Speed(idx)),clrs,'o',5);
gscatter(E.Day(idx),100*E.CellsInS(idx),E.Speed(idx),clrs,[],20);
title('ESC Oligomycin')
ylabel('% of cells in S')
legend('off')
xlabel('Day')

%% Fast vs Slow first  days ESCs
S = stack(E,{'Rel2DMSO_1' 'Rel2DMSO_2' } , 'NewDataVariableName' , 'S_rel2DMSO' ); 
%S.Drug( S.Drug=='AA+Oligomycin') = 'Oligomycin' ;
S.Drug = categorical( string(S.Drug) ) ;
S = S( S.Drug ~= 'DMSO' , :);
S = S( S.Day  <  3  ,:);

idx = S.Drug=='Antimycin' ;
[~,p] = ttest2( S.S_rel2DMSO(S.Speed=='Fast' & idx) , S.S_rel2DMSO(S.Speed=='Slow' & idx) )

figure; 
boxplot( S.S_rel2DMSO , {S.Drug S.Speed})
title('ESCs')
%% Fast vs Slow first two days FIBs
S = stack(F,{'Rel2DMSO_1' 'Rel2DMSO_2' 'Rel2DMSO_3'} , 'NewDataVariableName' , 'S_rel2DMSO' ); 
%S = stack(F,{'Replicate1' 'Replicate2' 'Replicate3'} , 'NewDataVariableName' , 'S' ); 
%S.Drug( S.Drug=='AA+Oligomycin') = 'Oligomycin' ;
S.Drug = categorical( string(S.Drug) ) ;
S = S( S.Drug ~= 'DMSO' , :);
S = S( S.Day  <  3  ,:);

idx = S.Drug=='Antimycin' ;
[~,p] = ttest2( S.S_rel2DMSO(S.Speed=='Fast' & idx) , S.S_rel2DMSO(S.Speed=='Slow' & idx) )

figure; 
boxplot( S.S_rel2DMSO , {S.Drug S.Speed})
title('FIBs')
%%

[p,tbl,stats,terms] = anovan(  S.S_rel2DMSO(idx) , {S.Drug(idx),S.Day(idx),S.Speed(idx)} ...
    ,'VarNames' , {'Drug' 'Day' 'Speed' }  ,'model','linear', 'Continuous',2)
c = multcompare( stats , 'Dimension' , [1 3])

%% publication quality figures
% show difference in % BrdU+, as expected from the sort
figname = './Sorted_Cells_BrdU_DMSO_drug_barplots' ; 

T=E(E.Drug=='DMSO',:);

figure('Position',[100 100 200 200]); 
x = [ T.Replicate1(T.Day==1 & T.Speed=='Fast') T.Replicate1(T.Day==1 & T.Speed=='Slow') ;
    T.Replicate1(T.Day==2 & T.Speed=='Fast') T.Replicate1(T.Day==2 & T.Speed=='Slow') ;
%    T.Replicate1(T.Day==3 & T.Speed=='Fast') T.Replicate1(T.Day==3 & T.Speed=='Slow') ;
    T.Replicate2(T.Day==1 & T.Speed=='Fast') T.Replicate2(T.Day==1 & T.Speed=='Slow') ;
    T.Replicate2(T.Day==2 & T.Speed=='Fast') T.Replicate2(T.Day==2 & T.Speed=='Slow') ;
%    T.Replicate1(T.Day==3 & T.Speed=='Fast') T.Replicate1(T.Day==3 & T.Speed=='Slow') ;
    ]
bar([1 2 4 5] , 100*x , 2 )
set(gca,'xticklabel',{'Day 1' 'Day2' 'Day 1' 'Day 2'})
xlabel('Rep. 1                    Rep. 2')
legend({'Fast' 'Slow'},'location','nw','box','off')
ylabel('% of cells in S   (BrdU^+)')
title('ESCs +DMSO')
ylim([0 100])
text(0.7,max(ylim)-3,'CFSE sort')
print('-dpng',[figname '_bar_DMSO_ESCs.png'] , '-r600')
close ; 

T=F(F.Drug=='DMSO',:);
figure('Position',[100 100 200 200]); 
x = [ T.Replicate1(T.Day==1 & T.Speed=='Fast') T.Replicate1(T.Day==1 & T.Speed=='Slow') ;
    T.Replicate1(T.Day==2 & T.Speed=='Fast') T.Replicate1(T.Day==2 & T.Speed=='Slow') ;
%    T.Replicate1(T.Day==3 & T.Speed=='Fast') T.Replicate1(T.Day==3 & T.Speed=='Slow') ;
    T.Replicate2(T.Day==1 & T.Speed=='Fast') T.Replicate2(T.Day==1 & T.Speed=='Slow') ;
    T.Replicate2(T.Day==2 & T.Speed=='Fast') T.Replicate2(T.Day==2 & T.Speed=='Slow') ;
%    T.Replicate1(T.Day==3 & T.Speed=='Fast') T.Replicate1(T.Day==3 & T.Speed=='Slow') ;
    ]
bar([1 2 4 5] , 100*x , 2 )
set(gca,'xticklabel',{'Day 1' 'Day2' 'Day 1' 'Day 2'})
xlabel('Rep. 1                    Rep. 2')
lh = legend({'Fast' 'Slow'},'location','nw','box','off')
ylabel('% of cells in S   (BrdU^+)')
title('FIBs +DMSO')
ylim([0 55])
text(0.7,max(ylim)-2.1,'CFSE sort')
print('-dpng',[figname '_bar_DMSO_FIBs.png'] , '-r600')
close ; 



T=F(F.Drug=='Antimycin',:);
figure('Position',[100 100 200 200]); 
x = [ T.Replicate1(T.Day==1 & T.Speed=='Fast') T.Replicate1(T.Day==1 & T.Speed=='Slow') ;
    T.Replicate1(T.Day==2 & T.Speed=='Fast') T.Replicate1(T.Day==2 & T.Speed=='Slow') ;
%    T.Replicate1(T.Day==3 & T.Speed=='Fast') T.Replicate1(T.Day==3 & T.Speed=='Slow') ;
    T.Replicate2(T.Day==1 & T.Speed=='Fast') T.Replicate2(T.Day==1 & T.Speed=='Slow') ;
    T.Replicate2(T.Day==2 & T.Speed=='Fast') T.Replicate2(T.Day==2 & T.Speed=='Slow') ;
%    T.Replicate1(T.Day==3 & T.Speed=='Fast') T.Replicate1(T.Day==3 & T.Speed=='Slow') ;
    ]
bar([1 2 4 5] , 100*x , 2 )
set(gca,'xticklabel',{'Day 1' 'Day2' 'Day 1' 'Day 2'})
xlabel('Rep. 1                    Rep. 2')
lh = legend({'Fast' 'Slow'},'location','nw','box','off')
ylabel('% of cells in S   (BrdU^+)')
title('FIBs +Antimycin')
ylim([0 55])
text(0.7,max(ylim)-2.1,'CFSE sort')
print('-dpng',[figname '_bar_Anti_FIBs.png'] , '-r600')
close ; 
%%
for drugname = {'Antimycin' 'Oligomycin'} 
    drugname = drugname{:};
T=F(F.Drug==drugname & F.Day<3,:);
figure('Position',[100 100 200 200]); 
hold on ;
clrs = get(gca,'ColorOrder');
v = vertcat(T.Rel2DMSO_1,T.Rel2DMSO_2,T.Rel2DMSO_3);
g = vertcat(T.Speed,T.Speed,T.Speed);
bh = boxplot(-1*v,g,'widths',0.8);
ylabel({'Drug effect' 'log_2(DMSO/drug)'});
[~,p] = ttest2( v(g=='Fast'),v(g=='Slow'));
[~,pF] = ttest( v(g=='Fast'));[~,pS] = ttest( v(g=='Slow'));
fprintf('FIB TTests for %s : fast= %0.05f , slow= %0.05f\n' , drugname , pF, pS);
title( sprintf('FIB +%s p=%0.04f', drugname , p ) )
line(xlim,[0 0],'LineStyle','--','Color',[.7 .7 .7])

h = findobj(gca,'Tag','Box');
patch(get(h(1),'XData'),get(h(1),'YData'),clrs(2,:),'FaceAlpha',1);
patch(get(h(2),'XData'),get(h(2),'YData'),clrs(1,:),'FaceAlpha',1);
bh = boxplot(-1*v,g,'Color','k','widths',0.8);
set(bh(6,1),'LineWidth',3,'Color','k')
set(bh(6,2),'LineWidth',3,'Color','k')
xlabel('CFSE sort')
ylim([-0.09 0.59])

print('-dpng',[figname '_box_' drugname '_FIBs.png'] , '-r600')
close; 


T=E(E.Drug==drugname & E.Day<3,:);
figure('Position',[100 100 200 200]); 
hold on ;
v = vertcat(T.Rel2DMSO_1,T.Rel2DMSO_2);
g = vertcat(T.Speed,T.Speed);
bh = boxplot(-1*v,g,'widths',0.8);
ylabel({'Drug effect' 'log_2(DMSO/drug)'});
[~,p] = ttest2( v(g=='Fast'),v(g=='Slow'));
[~,pF] = ttest( v(g=='Fast'));[~,pS] = ttest( v(g=='Slow'));
fprintf('ESC TTests for %s : fast= %0.05f , slow= %0.05f\n' , drugname , pF, pS);
title( sprintf('ESC +%s p=%0.04f', drugname , p ) )
line(xlim,[0 0],'LineStyle','--','Color',[.7 .7 .7])

h = findobj(gca,'Tag','Box');
patch(get(h(1),'XData'),get(h(1),'YData'),clrs(2,:),'FaceAlpha',1);
patch(get(h(2),'XData'),get(h(2),'YData'),clrs(1,:),'FaceAlpha',1);
bh = boxplot(-1*v,g,'Color','k','widths',0.8);
set(bh(6,1),'LineWidth',3,'Color','k')
set(bh(6,2),'LineWidth',3,'Color','k')
xlabel('CFSE sort')
ylim([-0.19 0.59])

print('-dpng',[figname '_box_' drugname '_ESCs.png'] , '-r600')
close; 

end

%%
T=F(F.Drug=='Antimycin',:);
figure('Position',[100 100 200 200]); 
x = [ T.Replicate1(T.Day==1 & T.Speed=='Fast') T.Replicate1(T.Day==1 & T.Speed=='Slow') ;
    T.Replicate1(T.Day==2 & T.Speed=='Fast') T.Replicate1(T.Day==2 & T.Speed=='Slow') ;
    T.Replicate2(T.Day==1 & T.Speed=='Fast') T.Replicate2(T.Day==1 & T.Speed=='Slow') ;
    T.Replicate2(T.Day==2 & T.Speed=='Fast') T.Replicate2(T.Day==2 & T.Speed=='Slow') ;
    ]
bar([1 2 4 5] , 100*x , 2 )
set(gca,'xticklabel',{'Day 1' 'Day2' 'Day 1' 'Day 2'})
xlabel('Rep. 1                    Rep. 2')
lh = legend({'Fast' 'Slow'},'location','nw','box','off')
ylabel('% of cells in S   (BrdU^+)')
title('FIBs +DMSO')
ylim([0 65])
text(0.7,max(ylim)-2.1,'CFSE sort')
print('-dpng',[figname '_DMSO_FIBs.png'] , '-r600')
close ; 