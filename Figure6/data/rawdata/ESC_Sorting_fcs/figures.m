DS.Day = str2double( regexprep(DS.PlateName,'DAY',''));
DS.filename = upper(DS.filename);
DS = sortrows(DS , {'Day' 'TMRE_bin' 'filename'} );
%DS = DS(DS.Day>0,:);
%% metadata
DS.LowO2 = NaN(nrows(DS),1);
DS.VitC = NaN(nrows(DS),1);

for I =1:nrows(DS)
    meta = regexprep( regexprep( DS.filename{I} ,'.FCS',''),'.*_P. ','');
    switch meta
        case{'1' '2'}
            DS.LowO2(I)= false ;
            DS.VitC(I) = false ;
        case{'AA' 'AA2' 'AA 2'}
            DS.LowO2(I) = false;
            DS.VitC(I) = true;
        case{'LOW 1' 'LOW 2'}
            DS.LowO2(I) = true;
            DS.VitC(I) = false ;
        case{'LOW AA' 'LOW AA2' 'LOW AA1' 'LOW AA1'}
            DS.LowO2(I) = true ;
            DS.VitC(I) = true ;
        otherwise
            disp( DS.filename{I})
            disp( meta) 
            disp(I)
    end
    if regexpcmp(meta,'1')
        DS.replicate(I) = 1;
    else
        DS.replicate(I) = 2;
    end
end
         
%%
DS.CFSE = DS.mean_FITC_A ; 
gscatter( DS.CFSE , DS.mean_PE_A ,strcat(DS.PlateName,'-',string(DS.TMRE_bin)) , 'bgk' , 'ooo+++')
xlabel('mean CFSE')
ylabel('mean TMRE')

%%
[P,T,STATS,TERMS]=anovan( DS.CFSE , [DS.TMRE_bin  DS.Day DS.VitC DS.LowO2 DS.replicate] ,'VarNames' ...
    ,  {'TMRE bin' 'Day' 'VitC' 'Low O2' 'rep'} ,'Display','on','model',[1 0 0 0 0 ; 0 1 0 0 0 ; 0 0 1 0 0 ; 0 0 0 1 0 ; 0 0 0 0 1 ]);
[c,m,h,gnames] = multcompare(STATS , 'dimension',[4 1] )
xlabel('mean CFSE')

%%
DS = sortrows(DS , {'Day' 'TMRE_bin' 'filename'} );
DS2 = DS( DS.Day==2,:);
DS3 = DS( DS.Day==3,:);

DS2.ratio = log2(DS3.mean_FITC_A ./ DS2.mean_FITC_A) ;

fh = figure('units','centimeters','position',[5 5 5 7 ]);
boxplot(DS2.ratio , DS2.TMRE_bin)
set(gca,'xticklabel',{'low' 'med' 'high'})
xlabel('TMRE bin')
ylabel('CFSE change (Day3 / Day2)')

%% plot CFSE for low & high TMRE in normal conditions : 
idxLow  = find( DS.Day>0 & DS.LowO2==false & DS.VitC==false & DS.TMRE_bin==6) ; 
idxHigh  = find( DS.Day>0 & DS.LowO2==false & DS.VitC==false & DS.TMRE_bin==8) ; 
idx = vertcat(idxLow,idxHigh);
clrs = [1 .5 0 ; 1 .5 0 ; 1 0 0 ; 1 0 0 ; 0 .5 1 ; 0 .5 1 ; 0 0 1 ; 0 0 1];
fh = figure('units','centimeters','position',[5 5 14 7 ]);
hold on ;
for I = 1:numel(idx)
    X = DS.data_PE_A{idx(I)};
    X = log2(X(X>0));
    [f,x] = ksdensity(X,1:0.01:20);
    txt = sprintf( 'Day%d TMRE %d' , DS.Day(idx(I)) , DS.TMRE_bin(idx(I)) )
    plot(x,f,'-','LineWidth',2,'Color',clrs(I,:) , 'DisplayName',  txt);
end
xlim([3 20])
legend('location','best')
xlabel('TMRE (log_2)')
%% Day0 is the sort  
% make some histograms
idx = find(DS.Day==0);
vn = DS.Properties.VarNames;
vn = vn( regexpcmp(vn,'^data_')  )
figname = 'day0 histograms.eps'; delete(figname);
for vni = 1:numel(vn)
    X = DS.(vn{vni}) ;
    fh = figure('units','centimeters','position',[5 5 7 10 ]);
    hold on ;
    for I = 1:numel(idx)
        x = X{idx(I)};
        x = log2(x(x>0));
        [f,x] = ksdensity(x);
        plot(x,f,'-','DisplayName',regexprep( DS.filename{idx(I)} , '_' ,' ') ,'LineWidth',3 );
    end
    legend('location','NorthOutside')
    axis tight;
    xlabel( regexprep( vn{vni} ,'_' ,' ') );
    print('-dpsc2',figname,'-append');
    close;
end

%%
Y = DS.data_FITC_A{ strcmp(DS.filename,'05122017_Day 0.fcs')} ;
X = DS.data_PE_A{ strcmp(DS.filename,'05122017_Day 0.fcs')} ;
idx = X>0 & Y>0;
X = X(idx);
Y = Y(idx);

Yl = log10(Y); 
idx = Yl < ( modefit( Yl ) *1.015) & Yl > ( modefit( Yl ) * 0.985 ) ; 
X = X(idx) ;
Y = Y(idx) ; 
P6 = X<700 & X>60 ; 
P7 = X>900 & X<1050 ;
P8 = X>1100 ; 

X = log10(X) ; 
Y = log10(Y) ; 
fh = figure('units','centimeters','position',[5 5 7 10 ]);
hold on; 
[f,x]  = ksdensity(Y(P6),1:0.001:5); 
plot(x,f,'-','DisplayName','P6 (low)' ,'LineWidth',3 );
[f,x]  = ksdensity(Y(P7),1:0.001:5); 
plot(x,f,'-','DisplayName','P7 (mid)' ,'LineWidth',3 );
[f,x]  = ksdensity(Y(P8),1:0.001:5); 
plot(x,f,'-','DisplayName','P8 (high)' ,'LineWidth',3 );

xlim([3.5 3.85])