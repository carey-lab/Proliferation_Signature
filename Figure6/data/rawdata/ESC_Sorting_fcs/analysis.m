cd('~/Google Drive/CareyLab/Projects/2016__ProliferationCorrelatedExpression/20171215__ESCs_with_TMRE');
%% something is wrong with the 'Sort' exported .fcs files. read in the files from .csv file
ESCells = readtable('Sort/ES_001.csv','Delimiter',',','FileType','text');
param = regexp('FSC,SSC,FITC,PE,PE-TR-PI,PECY5,640-APC,640-APC ALEXA,640-APC-CY7,561-PE,561-PE-CY5,561-PE-CY7,457-CFP,457-CROMOMCYM,355-DAPI,355' ,',','split');
%% Read data from tubes
Day0 = FACSReadPlateDir( 'DAY0' ,'fcsprefix','');
Day2 = FACSReadPlateDir( 'DAY2','fcsprefix','');
Day3 = FACSReadPlateDir( 'DAY3','fcsprefix','');
%Sorting = FACSReadPlateDir( 'Sort','fcsprefix','05');
wells = [Day0 Day2 Day3 ];
%SampleMap = containers.Map( {'156989.fcs' '156991.fcs' } , {'unstained1' 'unstained2'} );
%SampleMap('156997.fcs') = 'ES P9 Day0' ;
%SampleMap('156995.fcs') = 'ES P8 Day0' ;

%% extract sorted bin from wells 
for I = 1:numel(wells)
    fn = upper( wells(I).filename );
    fn = regexprep( fn ,'^.*_P','');
    wells(I).TMRE_bin = str2double( fn(1) ) ; 
end


%% FSC A x H
 Gates_FSC = FACSDrawGateForEach(wells,'group','fcstype','max_total_cells',10000,'xvar','FSC_A','yvar','FSC_H');
 wells = FACSApplyGates(wells,'fcstype',Gates_FSC);
 wells = arrayfun(@FACSRemoveFilteredCells,wells);
%%  live cells are the negative ones for TO PRO3 (APC dye).
 Gates_Live = FACSDrawGateForEach(wells,'group','fcstype','max_total_cells',10000,'yvar','APC_A','xvar','FSC_A');
 wells = FACSApplyGates(wells,'fcstype',Gates_Live);
 wells = arrayfun(@FACSRemoveFilteredCells,wells);
  Gates_Live = FACSDrawGateForEach(wells,'group','fcstype','max_total_cells',10000,'yvar','APC_A','xvar','SSC_A');
 wells = FACSApplyGates(wells,'fcstype',Gates_Live);
 wells = arrayfun(@FACSRemoveFilteredCells,wells);

 
 
%% calculate ratio measurement and statistics
wells = arrayfun(@(x)FACSCalcNewCellValue(x,'CFSEoverFSC','x.FITC_A ./ x.FSC_A') , wells);
wells = arrayfun(@(x)FACSCalcNewCellValue(x,'CFSEoverSSC','x.FITC_A ./ x.SSC_A') , wells);
wells = arrayfun(@(x)FACSCalcNewCellValue(x,'TMREoverSSC','x.PE_A ./ x.SSC_A') , wells);
wells = arrayfun(@(x)FACSCalcNewCellValue(x,'TMREoverFSC','x.PE_A ./ x.FSC_A') , wells);

vn = wells(1).data.Properties.VarNames;
vn = vn( regexpcmp(vn,'_A'));
vn  = [vn {'CFSEoverFSC' 'CFSEoverSSC' 'TMREoverSSC' 'TMREoverFSC' } ] ;

wells = arrayfun(@(x)FACSGetStats(x,'flunames',vn),wells);

%wells = arrayfun(@(x)FACSGammaFit( x,'fluname','RoverG' ) , wells) ;
%wells = arrayfun(@(x)FACSGammaFit( x,'fluname','PE_TexasRed_A' ) , wells) ;
%wells = arrayfun(@(x)FACSGammaFit( x,'fluname','FITC_A' ) , wells) ;



DS = struct2ds(wells);
save DS.mat DS

% 
% %%
% %%
% ud = unique({wells.PlateName});
% DS = wells(1).data ;
% DS.Day =  repmat(ud(1) , nrows(DS) , 1);
% DS.Bin = repmat( wells(1).TMRE_bin , nrows(DS) , 1);
% for day = 1:numel(ud)
%     idx = find(strcmp( {wells.PlateName} , ud{day}));
%     for I = 1:numel(idx)
%         dt = wells(idx(I)).data;
%         dt.Day =  repmat(ud(day) , nrows(dt) , 1);
%         dt.Bin = repmat( wells(idx(I)).TMRE_bin , nrows(dt) , 1);
%         DS = vertcat( DS , dt ) ;
%     end
% end
% %%
% vn = wells(1).data.Properties.VarNames;
% vn = vn( regexpcmp(vn,'_A'));
% for I = 1:numel(vn)
%     fh = figure('units','centimeters','position',[5 5 10 10 ]);
%     X =  DS.FSC_A  ;
%     Y =  DS.(vn{I}) ;
%     idx = X>0 & Y>0;
%     X = X(idx);
%     Y = Y(idx);
%     G =  strcat(DS.Day,'-',string(DS.Bin)) ;
%     idx = randsample(numel(X) , 1e4 ) ;
%     gscatter(  X(idx) , Y(idx) , G(idx));
%     set(gca,'yscale','log')
%     set(gca,'xscale','log')
%     xlim([ prctile(X,0.1) prctile(X,99.9)])
%     ylim([ prctile(Y,0.1) prctile(Y,99.9)])
%     xlabel('FSC-A')
%     ylabel(regexprep( vn{I} ,'_' ,'-'));
%     print('-dpng',[vn{I} '.png'] , '-r300');
%     close;
% end
