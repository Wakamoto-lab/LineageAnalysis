function plot_cumulant()
    %% List structure
    % Row : environments and strains
    % E.coli        from R1 to R12
    % S.pombe       from R13 to R19
    % M.smegmatis   R20
    
    % Column : analyzed value
    % C1 : Ttot(hr)
    % C2 : Lambda(1/hr)
    % C3 : number of initial cells
    % C4 : number of last cells
    % C5 - C10 : cumulant derivated by n!, c_n
    % C11 - C16 : scaled cumulant, w_n
    %% data import
    disp('choose DataSummary.xlsx');
    [file,path] = uigetfile({'*.xlsx';'*.csv';'*.mat';'*.*'},'File Selector');
    cd(path);
    data = xlsread(file,'G2:AM20');
    
    %% data distribution
    T_tot = data(:,1);
    Lambda= data(:,2);
    err_Lambda = data(:,3);
    N0Ntau = data(:,4:5);
    C = [data(:,8),data(:,10),data(:,12),data(:,14),data(:,16),data(:,18)];
    err_C = [data(:,9),data(:,11),data(:,13),data(:,15),data(:,17),data(:,19)];
    W = [data(:,20),data(:,22),data(:,24),data(:,26),data(:,28),data(:,30)];
    err_W = [data(:,21),data(:,23),data(:,25),data(:,27),data(:,29),data(:,31)];
    wSD = data(:,32);
    err_wSD = data(:,33);
    

%
    
    %% sort for Lambda in each spieces
    [~,sortEc] = sort(Lambda(1:11,1));
    sortEc = flipud(sortEc);
    xEc = 0.85:0.03:1.15;
    [~,sortSp] = sort(Lambda(12:18,1));
    sortSp = flipud(sortSp) + 11;
    xSp = 1.41:0.03:1.59;
    
    %% plot (strain, Lambda)
    figure(7);clf;
    % E.coli
    errorbar(xEc,Lambda(sortEc,1),2*err_Lambda(sortEc,1),'ro',...
    'MarkerFaceColor','r','MarkerEdgeColor','k','Color','k','LineWidth', 1.5,'MarkerSize',7);hold on;
    % S.pombe
    errorbar(xSp,Lambda(sortSp,1),2*err_Lambda(sortSp,1),'bo',...
    'MarkerFaceColor','b','MarkerEdgeColor','k','Color','k','LineWidth', 1.5,'MarkerSize',7);
    % M.smegmatis
    errorbar(2,Lambda(19,1),2*err_Lambda(19,1),'go',...
    'MarkerFaceColor','g','MarkerEdgeColor','k','Color','k','LineWidth', 1.5,'MarkerSize',7);
    xticks([1 1.5 2]);
    xticklabels({'\sl E. coli','\sl S. pombe','\sl M. smegmatis'});
    xtickangle(45);
    
    % ylabel
    ylabel('$$\Lambda$$','Interpreter', 'Latex');
    
    % xlim and ylim
    xlim([0.7  2.2]);
    ylim([0 1.2]);
    
    % save figure
    savefig('strain_lambda.fig');    
    
    %% plot (strain, 1-w1)
    figure(8);clf;
    % E.coli
    errorbar(xEc,wSD(sortEc,1),2*err_wSD(sortEc,1),'ro',...
    'MarkerFaceColor','r','MarkerEdgeColor','k','Color','k','LineWidth', 1.5,'MarkerSize',7);hold on;
    % S.pombe
    errorbar(xSp,wSD(sortSp,1),2*err_wSD(sortSp,1),'bo',...
    'MarkerFaceColor','b','MarkerEdgeColor','k','Color','k','LineWidth', 1.5,'MarkerSize',7);
    % M.smegmatis
    errorbar(2,wSD(19,1),2*err_wSD(19,1),'go',...
    'MarkerFaceColor','g','MarkerEdgeColor','k','Color','k','LineWidth', 1.5,'MarkerSize',7);
    xticks([1 1.5 2]);
    xticklabels({'\sl E. coli','\sl S. pombe','\sl M. smegmatis'});
    xtickangle(45);    
    
    % ylabel
    ylabel('$$S[D]/\tau\Lambda$$','Interpreter', 'Latex');
    
    % xlim and ylim
    xlim([0.7  2.2]);
    ylim([0 0.12]);    
    
    % save figure
    savefig('strain_wSD.fig');   
  
    %% plot (strain, 1/w1)
    figure(9);clf;
    % E.coli
    errorbar(xEc,1./(1-wSD(sortEc,1)),2*err_wSD(sortEc,1)./((1-wSD(sortEc,1)).^2),'ro',...
    'MarkerFaceColor','r','MarkerEdgeColor','k','Color','k','LineWidth', 1.5,'MarkerSize',7);hold on;
    % S.pombe
    errorbar(xSp,1./(1-wSD(sortSp,1)),2*err_wSD(sortSp,1)./((1-wSD(sortSp,1)).^2),'bo',...
    'MarkerFaceColor','b','MarkerEdgeColor','k','Color','k','LineWidth', 1.5,'MarkerSize',7);
    % M.smegmatis
    errorbar(2,1/(1-wSD(19,1)),2*err_wSD(19,1)/((1-wSD(19,1))^2),'go',...
    'MarkerFaceColor','g','MarkerEdgeColor','k','Color','k','LineWidth', 1.5,'MarkerSize',7);
    xticks([1 1.5 2]);
    xticklabels({'\sl E. coli','\sl S. pombe','\sl M. smegmatis'});
    xtickangle(45);    
    
    % ylabel
    ylabel('$$ \Lambda/c_1$$','Interpreter', 'Latex');
    
    % xlim and ylim
    xlim([0.7  2.2]);
    ylim([1 1.15]);    
    
    % save figure
    savefig('strain_Lamperc1.fig');    
%     %% plot (strain, w1)
%     figure(6);clf;
%     DataArray=nan(12,3);
%     % E.coli
%     DataArray(1:11,1)=List(1:11,11);
%     % S.pombe
%     DataArray(1:7,2)=List(12:18,11);
%     % M.smegmatis
%     DataArray(1,3)=List(19,11);
%     Colors=[0.9 0 0;0 0.9 0.9; 0 0.9 0];
%     
%     % plot and xlabel
%     UnivarScatter(DataArray,'Label',{'\sl E. coli','\sl S. pombe','\sl M. smegmatis'},'MarkerFaceColor',Colors);
%     xtickangle(45);
%     
%     % ylabel
%     ylabel('$$w_1$$','Interpreter', 'Latex');
%     
%     % xlim and ylim
%     xlim([0.5  3.5]);
%     ylim([0.88 1]);
%     
%     % save figure 
%     savefig('strain_w1.fig');
%     saveas(gcf,'strain_w1.eps','epsc');
end

function randx = make_randX(nx)
    randx = zeros(1,nx);
    randx(1,1) = randn(1,1)/10;
    for ix = 2:nx
        min_dist = 0.001;
        r = 0;
        while min_dist < 0.015
            r = randn(1,1)/10;
            min_dist = abs(r-randx(1,1));
            for jx = 1:(ix-1)
                min_dist = min(min_dist, abs(r-randx(1,jx)));
            end
        end
        randx(1,ix) = r;
    end
end