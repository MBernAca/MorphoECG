
function plot_results(geol,wc,As,geolmap,axy,labelax,TableName)
    % This function computes and plot power-law fits for width and depth
    % against drainage area, for each lithological domain.
    
    % Compute log10
    wc_log2 = log10(wc);
    As_log = log10(As);
    
    % Compute mean bin width
    Nbin =round(sqrt(mean([numel(As_log(geol==2)),numel(As_log(geol==3)),numel(As_log(geol==4)),numel(As_log(geol==5))])));
    %-------------- BASEMENT ---------------------------
    ind=find(geol==2); 
    edges = min(As_log):(max(As_log)-min(As_log))/Nbin:max(As_log);
    figure;set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    sub1=subplot(1,4,1);hold on;
    xlabel(labelax{1});ylabel(labelax{2});xlim([min(As_log) max(As_log)]); ylim([min(wc_log2) max(wc_log2)+0.5]);
    % get bin number of each data and count number of data in each bin
    [~,~,loc] = histcounts(As_log(ind),edges); SumBin=accumarray(loc(:),1,[Nbin 1]);
    % Compute the mean value in each bin
    meany = accumarray(loc(:),wc_log2(ind),[Nbin 1])./SumBin; meany(isnan(meany))=0;
    % Centered x for each bin
    xmid=0.5.*(edges(1:end-1)+edges(2:end));

    meany2=meany;meany2(meany==0)=[];xmid(meany==0)=[];SumBin(meany==0)=[];
    % Fit power law weighted by the number of count in each bin
    mdl = fitlm(xmid,meany2,'Weights',SumBin)
    %Get Coefficient (intercept and exponent)
    coef = mdl.Coefficients{:,'Estimate'};
    std_coef = mdl.Coefficients{:,'SE'};
    tStat_coef = mdl.Coefficients{:,'tStat'};
    pValue_coef = mdl.Coefficients{:,'pValue'};
    
    % Plot results
    subplot(1,4,1);scatter(As_log(ind),10.^wc_log2(ind),10,geol(ind),'fill');colormap(geolmap);caxis([0 5]);title('Basement'); ylim([0 max(10.^wc_log2)]);xlim([0 max(As_log)]); axis square;
    axis square
    plot(xmid,10.^meany2,'ko','MarkerFaceColor','w');
    % Plot fit
    [yfitted,CI] = predict(mdl,xmid');
    plot(xmid,10.^yfitted,'k-'); Rs1=mdl.Rsquared.Adjusted;
    plot(xmid,10.^yfitted+10.^CI(:,2),'k--','LineWidth',0.5);
    plot(xmid,10.^yfitted-10.^CI(:,1),'k--','LineWidth',0.5);
    
    % save coefficients
    a1=10^coef(1);b1=coef(2);CI1=coefCI(mdl);
    std1 = 10^std_coef(1); std1b = std_coef(2);
    table1 = table([a1;b1],[std1;std1b],[tStat_coef(1);tStat_coef(2)],[pValue_coef(1);pValue_coef(2)],'VariableNames',{'coefficients','SE','tStat','pValue'},'RowNames',{'1','2'});
    ax = get(gca,'XTick'); ay = get(gca,'YTick');
    text(ax(axy(1)),ay(axy(2)-1),['R^2 =',num2str(Rs1)],'Parent',sub1,'FontSize',6);
    text(ax(axy(1)),ay(axy(2)),['y = ',num2str(a1),'A^{',num2str(b1),'}'],'Parent',sub1,'FontSize',6);
    hold off;
    
    %---------------  METASEDIMENT  --------------------------------------
    ind=find(geol==4);
    sub2=subplot(1,4,2);hold on;
    xlabel(labelax{1});ylabel(labelax{2}); xlim([min(As_log) max(As_log)]); ylim([min(wc_log2) max(wc_log2)+0.5]);
    % get bin number of each data and count number of data in each bin
    [~,~,loc] = histcounts(As_log(ind),edges); SumBin=accumarray(loc(:),1,[Nbin 1]);
    % Compute the mean value in each bin
    meany = accumarray(loc(:),wc_log2(ind),[Nbin 1])./SumBin; meany(isnan(meany))=0;
    % Centered x for each bin
    xmid=0.5.*(edges(1:end-1)+edges(2:end));
    
    meany2=meany;meany2(meany==0)=[];xmid(meany==0)=[]; SumBin(meany==0)=[];
    % Fit power law weighted by the number of count in each bin
    mdl = fitlm(xmid,meany2,'Weights',SumBin)
    %Get Coefficient (intercept and exponent)
    coef = mdl.Coefficients{:,'Estimate'};
    std_coef = mdl.Coefficients{:,'SE'};
    tStat_coef = mdl.Coefficients{:,'tStat'};
    pValue_coef = mdl.Coefficients{:,'pValue'};
    % Plot results
    subplot(1,4,2);scatter(As_log(ind),10.^wc_log2(ind),10,geol(ind),'fill');colormap(geolmap);caxis([0 5]);title('Metasediment'); ylim([0 max(10.^wc_log2)]);xlim([0 max(As_log)]); axis square;
    axis square
    plot(xmid,10.^meany2,'ko','MarkerFaceColor','w');
    % Plot fit
    [yfitted,CI] = predict(mdl,xmid');
    plot(xmid,10.^yfitted,'k-'); Rs2=mdl.Rsquared.Adjusted; 
    plot(xmid,10.^yfitted+10.^CI(:,2),'k--','LineWidth',0.5);
    plot(xmid,10.^yfitted-10.^CI(:,1),'k--','LineWidth',0.5);    
    
    % save coefficients
    a2=10^coef(1);b2=coef(2);CI2=coefCI(mdl);
    std2 = 10^std_coef(1); std2b = std_coef(2);
    table2 = table([a2;b2],[std2;std2b],[tStat_coef(1);tStat_coef(2)],[pValue_coef(1);pValue_coef(2)],'VariableNames',{'coefficients','SE','tStat','pValue'},'RowNames',{'1','2'});
    ax = get(gca,'XTick'); ay = get(gca,'YTick');
    text(ax(axy(1)),ay(axy(2)-1),['R^2 =',num2str(Rs2)],'Parent',sub2,'FontSize',6);
    text(ax(axy(1)),ay(axy(2)),['y = ',num2str(a2),'A^{',num2str(b2),'}'],'Parent',sub2,'FontSize',6);
    hold off;
    
    
    % ------------ VOLCANICS ---------------------------------------------
    ind=find(geol==3); 
    sub3=subplot(1,4,3);hold on;
    xlabel(labelax{1});ylabel(labelax{2});xlim([min(As_log) max(As_log)]); ylim([min(wc_log2) max(wc_log2)+0.5]);
    % get bin number of each data and count number of data in each bin
    [~,~,loc] = histcounts(As_log(ind),edges); SumBin=accumarray(loc(:),1,[Nbin 1]);
    % Compute the mean value in each bin
    meany = accumarray(loc(:),wc_log2(ind),[Nbin 1])./SumBin; meany(isnan(meany))=0; 
    % Centered x for each bin
    xmid=0.5.*(edges(1:end-1)+edges(2:end));

    meany2=meany;meany2(meany==0)=[];xmid(meany==0)=[];SumBin(meany==0)=[];
    % Fit power law weighted by the number of count in each bin
    mdl = fitlm(xmid,meany2,'Weights',SumBin);
    %Get Coefficient (intercept and exponent)
    coef = mdl.Coefficients{:,'Estimate'};
    std_coef = mdl.Coefficients{:,'SE'};
    tStat_coef = mdl.Coefficients{:,'tStat'};
    pValue_coef = mdl.Coefficients{:,'pValue'};
    % Plot results
    subplot(1,4,3);scatter(As_log(ind),10.^wc_log2(ind),10,geol(ind),'fill');colormap(geolmap);caxis([0 5]);title('Volcanics'); ylim([0 max(10.^wc_log2)]);xlim([0 max(As_log)]); axis square;
    axis square
    plot(xmid,10.^meany2,'ko','MarkerFaceColor','w');
    % Plot fit
    [yfitted,CI] = predict(mdl,xmid');
    plot(xmid,10.^yfitted,'k-'); Rs3=mdl.Rsquared.Adjusted;
    plot(xmid,10.^yfitted+10.^CI(:,2),'k--','LineWidth',0.5);
    plot(xmid,10.^yfitted-10.^CI(:,1),'k--','LineWidth',0.5);
    
    % save coefficients
    a3=10^coef(1);b3=coef(2);CI3=coefCI(mdl);
    std3 = 10^std_coef(1); std3b = std_coef(2);
    table3 = table([a3;b3],[std3;std3b],[tStat_coef(1);tStat_coef(2)],[pValue_coef(1);pValue_coef(2)],'VariableNames',{'coefficients','SE','tStat','pValue'},'RowNames',{'1','2'});
    ax = get(gca,'XTick'); ay = get(gca,'YTick');
    text(ax(axy(1)),ay(axy(2)-1),['R^2 =',num2str(Rs3)],'Parent',sub3,'FontSize',6);
    text(ax(axy(1)),ay(axy(2)),['y = ',num2str(a3),'A^{',num2str(b3),'}'],'Parent',sub3,'FontSize',6);
    hold off;
    
    
    %----------------- SEDIMENTARY ---------------------------------------
    ind=find(geol==5);  
    sub4=subplot(1,4,4);hold on;
    xlabel(labelax{1});ylabel(labelax{2});xlim([min(As_log) max(As_log)]); ylim([min(wc_log2) max(wc_log2)+0.5]);
    % get bin number of each data and count number of data in each bin
    [~,~,loc] = histcounts(As_log(ind),edges); SumBin=accumarray(loc(:),1,[Nbin 1]);
    % Compute the mean value in each bin
    meany = accumarray(loc(:),wc_log2(ind),[Nbin 1])./SumBin; meany(isnan(meany))=0; 
    % Centered x for each bin
    xmid=0.5.*(edges(1:end-1)+edges(2:end));
    
    meany2=meany;meany2(meany==0)=[];xmid(meany==0)=[];SumBin(meany==0)=[];
    % Fit power law weighted by the number of count in each bin
    mdl = fitlm(xmid,meany2,'Weights',SumBin);
    % Get Coefficient (intercept and exponent)
    coef = mdl.Coefficients{:,'Estimate'}; 
    std_coef = mdl.Coefficients{:,'SE'};
    tStat_coef = mdl.Coefficients{:,'tStat'};
    pValue_coef = mdl.Coefficients{:,'pValue'};
    % Plot results
    subplot(1,4,4);scatter(As_log(ind),10.^wc_log2(ind),10,geol(ind),'fill');colormap(geolmap);caxis([0 5]);title('Sediments'); ylim([0 max(10.^wc_log2)]);xlim([0 max(As_log)]); axis square;xlabel('Drainage area(km^2)');ylabel('Valley width (km)');
    axis square
    plot(xmid,10.^meany2,'ko','MarkerFaceColor','w');
    % Plot fit
    [yfitted,CI] = predict(mdl,xmid');
    plot(xmid,10.^yfitted,'k-'); Rs4=mdl.Rsquared.Adjusted; 
    plot(xmid,10.^yfitted+10.^CI(:,2),'k--','LineWidth',0.5);
    plot(xmid,10.^yfitted-10.^CI(:,1),'k--','LineWidth',0.5);
    
    %save coefficients
    a4=10^coef(1);b4=coef(2);CI4=coefCI(mdl);
    std4 = 10^std_coef(1); std4b = std_coef(2);
    table4 = table([a4;b4],[std4;std4b],[tStat_coef(1);tStat_coef(2)],[pValue_coef(1);pValue_coef(2)],'VariableNames',{'coefficients','SE','tStat','pValue'},'RowNames',{'1','2'});
    ax = get(gca,'XTick'); ay = get(gca,'YTick');
    text(ax(axy(1)),ay(axy(2)-1),['R^2 =',num2str(Rs4)],'Parent',sub4,'FontSize',6);
    text(ax(axy(1)),ay(axy(2)),['y = ',num2str(a4),'A^{',num2str(b4),'}'],'Parent',sub4,'FontSize',6);
    hold off

    %Save coefficients and confidence intervals
    CIm = [10^CI1(1,1);CI1(2,1)];CIp=[10^CI1(1,2);CI1(2,2)] ;
    Rs=[Rs1;0];
    table1 = [table1,table(Rs),table(CIm),table(CIp)];
    CIm = [10^CI2(1,1);CI2(2,1)];CIp=[10^CI2(1,2);CI2(2,2)] ;
    Rs=[Rs2;0];
    table2 = [table2,table(Rs),table(CIm),table(CIp)];
    CIm = [10^CI3(1,1);CI3(2,1)];CIp=[10^CI3(1,2);CI3(2,2)] ;
    Rs=[Rs3;0];
    table3 = [table3,table(Rs),table(CIm),table(CIp)];
    CIm = [10^CI4(1,1);CI4(2,1)];CIp=[10^CI4(1,2);CI4(2,2)] ;
    Rs=[Rs4;0];
    table4 = [table4,table(Rs),table(CIm),table(CIp)];

    writetable(table1,TableName,'Sheet',1);
    writetable(table2,TableName,'Sheet',2);
    writetable(table3,TableName,'Sheet',3);
    writetable(table4,TableName,'Sheet',4);

    end
