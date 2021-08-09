function [outFLx,outM,outC, outW]=onedim_calculus(data,motherdir,t_tot,i)
    cd(motherdir);
    %% data import
    ndiv=data(1,:);
    x=data(i,:);
    char = phenotypename(i);
    
    %% determine bin center
    if i==1
        xc = min(x):1:max(x);
        nbinx = max(x) - min(x) + 1;
    else
        alpha = 0.4; % binwidth = alpha * iqr(x);
        wbinx=alpha*iqr(x);
        hx=histogram(x); % generate the histogram (not normalized)
        hx.BinWidth=wbinx; % set binwidth with wbinx
        nbinx=hx.NumBins; % correct the bin number
        xedges=hx.BinEdges;% generate the bin edges
        %xcenter
        xc=zeros(1,nbinx);
        for ibinx=1:nbinx
            xc(1,ibinx)=(xedges(ibinx)+xedges(ibinx+1))/2;
        end
    end
    
    %% calculate cl/rs distributions, fitness lanscape
    [cl,rs,N0,nlineage]=onedim_clrs_distr_given_bincenter(x,ndiv,xc,i);
    Lambda = log(nlineage/N0);
    SD = onedim_selection_strength_D(ndiv);
    FLx = onedim_fitness_landscape(cl,rs,Lambda,t_tot);
    Sx = onedim_selection_strength(cl,rs,xc);
    SxperSD = [Sx(1)/SD(1), Sx(2)/SD(2)];
    M = [Lambda, Sx, SxperSD];
    
    %% bootstarp errorbar estimation
    [errFLx,meanFLx, errM,meanM, errW, meanW] = onedim_calculus_bootstrap(xc,x,ndiv,t_tot,i);
    %onedim_calculus_shuffled_bootstrap(xc,x,ndiv,t_tot,i, M);
    %meanFLx = zeros(1,nbinx);
    %errFLx = zeros(1,nbinx);
    %meanM = zeros(1,5);
    %errM = zeros(1,5);
    outFLx = [xc; cl; FLx; meanFLx; errFLx];
    outM = [M; meanM; errM];

    %% H-function
    [outC, outW] = Hfunction(FLx*t_tot,cl,xc,i);

    outW = [outW; meanW; errW];
end

function [cldis_x,rsdis_x,N0,nlineage]=onedim_clrs_distr_given_bincenter(x,ndiv,xc,i)
    %% bin setting
    [~,nlineage]=size(x); % get the numer of lineages
    hx=histogram(x); % generate the histogram (not normalized)
    [~,nbinx]=size(xc);
    xedges=zeros(1,nbinx+1);
    wbinx=xc(2)-xc(1);
    xedges(1)=xc(1)-wbinx/2;
    for ibinx=1:nbinx
        xedges(ibinx+1)=xedges(ibinx)+wbinx;
    end
    hx.BinEdges=xedges;
    
    %% calculate the retrospective distribution
    rsdis_x=hx.BinCounts/nlineage; % normalize uniformaly
    rsdis_x=rsdis_x/wbinx; % convert frequency to density

    %% calculate the choronological distribution    
    cldis_x=zeros(1,nbinx);
    for ibinx=1:(nbinx-1)
        for l=1:nlineage
            if(x(l)>=xedges(ibinx) && x(l)<xedges(ibinx+1))
                cldis_x(ibinx)= cldis_x(ibinx)+2^(-ndiv(l));
            end
        end
    end
    for l=1:nlineage
        if (x(l)>=xedges(nbinx) && x(l)<=xedges(nbinx+1))
            cldis_x(nbinx)= cldis_x(nbinx)+2^(-ndiv(l));
        end
    end
    N0 = sum(cldis_x);
    cldis_x=cldis_x/sum(cldis_x); % normalize uniformaly
    cldis_x=cldis_x/wbinx; % convert frequency to density
end

function SD = onedim_selection_strength_D(ndiv)
    [~,nlineage]= size(ndiv);
    SD = [0,0];
    ch = 2.^(-ndiv);
    ch = ch/sum(ch);
    rs = ones(1,nlineage)/nlineage;

    for i=1:nlineage
        if rs(i)*ch(i) ~=0
            SD(1) = SD(1) + ch(i)*log(ch(i)/rs(i));
            SD(2) = SD(2) + rs(i)*log(rs(i)/ch(i));
        end
    end
end

function Sx = onedim_selection_strength(ch,rs,xc)
    [~,nbinx]= size(ch);
    Sx = [0,0];
    wbinx=xc(2)-xc(1);
    for i=1:nbinx
        if rs(i)*ch(i) ~=0
            Sx(1) = Sx(1) + ch(i)*log(ch(i)/rs(i))*wbinx;
            Sx(2) = Sx(2) + rs(i)*log(rs(i)/ch(i))*wbinx;
        end
    end
end

function FLx = onedim_fitness_landscape(cl,rs,Lambda,t_tot)
    %% get the number of bins
    [~,nbinx]=size(cl);
    %% calculate the fitnelss landscape
    % intialize
    FLx=zeros(1,nbinx);
    for i=1:nbinx
        if cl(i) ~= 0
            FLx(i) = Lambda + log(rs(i)/cl(i)); 
        end
    end
    FLx = FLx/t_tot;
end

function [errFLx, meanFLx, errM, meanM, errW, meanW] = onedim_calculus_bootstrap(xc, x, ndiv,t_tot,i)
    nBoot = 20000;
    nLineage = size(ndiv,2);
    nBin = size(xc,2);
    %% randomized (ndiv, x)s
    % calculate cumulative chronoloigial distribution
    minD = min(ndiv);
    Pcl = 2.^(-(ndiv - minD));
    Pcl = Pcl/sum(Pcl);
    cumPcl = cumsum(Pcl,2);
    %% boot arrays
    FLx_list = zeros(nBoot, nBin);
    M_list = zeros(nBoot,5);
    W_list = zeros(nBoot,6);

    %% bootstrap
    alert = 10;
    for iBoot = 1:nBoot
        % alert 
        if (iBoot/nBoot)*100 >= alert
            disp(horzcat('bootstrap: ',num2str(alert),'%....'));
            alert = alert + 10;
        end
        % create random sample dataset
        ndiv_boot = zeros(1,nLineage);
        x_boot = zeros(1,nLineage);
        for iLineage = 1:nLineage
            index = randi(nLineage);
            ndiv_boot(1,iLineage) = ndiv(1,index);
            x_boot(1,iLineage) = x(1,index);
        end
        [cl,rs,N0,nlineage]=onedim_clrs_distr_given_bincenter(x_boot,ndiv_boot,xc,i);
        Lambda = log(nlineage/N0);
        FLx = onedim_fitness_landscape(cl,rs,Lambda,t_tot);
        FLx_list(iBoot,:) = FLx;
        Sx = onedim_selection_strength(cl,rs,xc);
        SD = onedim_selection_strength_D(ndiv_boot);
        M_list(iBoot,:) = [Lambda, Sx, Sx(1)/SD(1), Sx(2)/SD(2)];
        W = cumulant_weight(FLx*t_tot,cl,xc,Lambda);
        W_list(iBoot,:) = W;
    end
    figure(1);
    hist(M_list(:,1));
    savefig( horzcat('hist_LambdaBootList_',phenotypename(i),'.fig'));
    hist(M_list(:,2));
    savefig(horzcat('hist_S1XBootList_',phenotypename(i),'.fig'));
    hist(M_list(:,3));
    savefig(horzcat('hist_S2XBootList_',phenotypename(i),'.fig'));
    hist(M_list(:,4));
    savefig(horzcat('hist_S1XperS1DBootList_',phenotypename(i),'.fig'));
    hist(M_list(:,5));
    savefig(horzcat('hist_S2XperS2DBootList_',phenotypename(i),'.fig'));
    clf;
    %% calculate 2SD for bootstrap samples
    % fitness lanscape
    meanFLx = mean(FLx_list);
    errFLx = zeros(1,nBin);
    for iBin = 1:nBin
        index = find(FLx_list(:,iBin) ~=0);
        if size(index,1)/nBoot >=0.95
            errFLx(1,iBin) = 2*std(FLx_list(index,iBin));
        end
    end
    % selection strength and population growth rate
    errM = 2*std(M_list,0,1);
    meanM = mean(M_list,1);
    
    errW = 2*std(W_list,0,1);
    meanW = mean(W_list,1);
end

function onedim_calculus_shuffled_bootstrap(xc, x, ndiv,t_tot,i, M)
    nBoot = 20000;
    nLineage = size(ndiv,2);
    nBin = size(xc,2);

    %% boot arrays
    M_list = zeros(nBoot,5);

    %% bootstrap
    alert = 10;
    for iBoot = 1:nBoot
        % alert 
        if (iBoot/nBoot)*100 >= alert
            disp(horzcat('bootstrap: ',num2str(alert),'%....'));
            alert = alert + 10;
        end
        % create random sample dataset
        x_boot = x(randperm(nLineage));

        % calculate cl and rs for permutated pairs of (D, X)
        [cl,rs,N0,nlineage]=onedim_clrs_distr_given_bincenter(x_boot,ndiv,xc,i);
        Lambda = log(nlineage/N0);
        Sx = onedim_selection_strength(cl,rs,xc);
        SD = onedim_selection_strength_D(ndiv);
        M_list(iBoot,:) = [Lambda, Sx, Sx(1)/SD(1), Sx(2)/SD(2)];
    end
    %% output figure
    figure(1);
    % Lambda
    hist(M_list(:,1));
    hold on;
    xline(M(1));
    savefig( horzcat('hist_Lambda_shuffled_BootList_',phenotypename(i),'.fig'));
    clf;
    % S1
    hist(M_list(:,2));
    hold on;
    xline(M(2));
    savefig(horzcat('hist_S1X_shuffled_BootList_',phenotypename(i),'.fig'));
    clf;
    % S2
    hist(M_list(:,3));
    hold on;
    xline(M(3));
    savefig(horzcat('hist_S2X_shuffled_BootList_',phenotypename(i),'.fig'));
    clf;
    % S1rel
    hist(M_list(:,4));
    hold on;
    xline(M(4));
    savefig(horzcat('hist_S1relX_shuffled_BootList_',phenotypename(i),'.fig'));
    clf;
    % S2rel
    hist(M_list(:,5));
    hold on;
    xline(M(5));
    savefig(horzcat('hist_S1relX_shuffled_BootList_',phenotypename(i),'.fig'));
    clf;
end

function [outC, outW] =  Hfunction(FLx,cl,xc,ipheno)
    figure(10);
    clf;
    beta = 0:0.01:1;
    H = ones(1,101);
    nbinx = size(FLx,2);
    wbinx = xc(2) - xc(1);
    for i = 0:100
        j = i + 1;
        b = beta(j);
        expH = 0;
        for ibinx = 1:nbinx
            expH = expH + exp(b*FLx(1,ibinx)) * cl(ibinx) * wbinx;
        end
        H(j) = log(expH);
    end
    plot(beta,H);
    hold on;

    % approx line by interpolation by beta = 0 and 1
    linear = H(1,1) * ones(1,101) + (0:0.01:1)*(H(1,101)- H(1,1));
    plot(beta,linear,'k--');
    grid on;
    ylim([0,ceil(H(1,101))+1]);
    savefig(horzcat('Hfunction_',phenotypename(ipheno),'.fig'));
    
    figure(12);
    clf;
    plot(beta, H-linear);
    savefig(horzcat('diff_Hfunction_',phenotypename(ipheno),'.fig'));

    %% calculate the 1st to 6th cumulants
    m=zeros(1,6);   % array for momments
    c=zeros(1,6);   % array for cumulants
    
    % calculate moments
    for order = 1:6
        for ibinx = 1:nbinx
            m(order) = m(order) + ((FLx(ibinx))^order) * cl(ibinx)*wbinx;
        end
    end
    
    % calculate cumulants
    c(1) = m(1);
    c(2) = (m(2) - m(1)^2)/2;
    c(3) = (m(3) - 3*m(2)*m(1) + 2*m(1)^3)/6;
    c(4) = (m(4) - 4*m(3)*m(1) -3*m(2)^2 + 12*m(2)*m(1)^2 - 6*m(1)^4)/24;
    c(5) = (m(5)-5*m(4)*m(1)-5*m(3)*m(2)+15*m(3)*m(1)^2+15*m(2)^2*m(1)-35*m(2)*m(1)^3+14*m(1)^5)/120;
    c(6) = (m(6)-6*m(5)*m(1)-15*m(4)*m(2)+30*m(4)*m(1)^2-10*m(3)^2+120*m(3)*m(2)*m(1)-120*m(3)*m(1)^3+30*m(2)^3-270*m(2)^2*m(1)^2+360*m(2)*m(1)^4-120*m(1)^6)/720;

    outC = cumsum(c);
    outW = outC/H(101);

    figure(13);
    clf;
    plot(1:6,outW);
    savefig(horzcat('Weight_',phenotypename(ipheno),'.fig'));
end

function W = cumulant_weight(FLx,cl,xc,Lambda)

    %% calculate the 1st to 6th cumulants
    m=zeros(1,6);   % array for momments
    c=zeros(1,6);   % array for cumulants

    nbinx = size(FLx,2);
    wbinx = xc(2) - xc(1);
    
    % calculate moments
    for order = 1:6
        for ibinx = 1:nbinx
            m(order) = m(order) + ((FLx(ibinx))^order) * cl(ibinx)*wbinx;
        end
    end

    % calculate cumulants
    c(1) = m(1);
    c(2) = (m(2) - m(1)^2)/2;
    c(3) = (m(3) - 3*m(2)*m(1) + 2*m(1)^3)/6;
    c(4) = (m(4) - 4*m(3)*m(1) -3*m(2)^2 + 12*m(2)*m(1)^2 - 6*m(1)^4)/24;
    c(5) = (m(5)-5*m(4)*m(1)-5*m(3)*m(2)+15*m(3)*m(1)^2+15*m(2)^2*m(1)-35*m(2)*m(1)^3+14*m(1)^5)/120;
    c(6) = (m(6)-6*m(5)*m(1)-15*m(4)*m(2)+30*m(4)*m(1)^2-10*m(3)^2+120*m(3)*m(2)*m(1)-120*m(3)*m(1)^3+30*m(2)^3-270*m(2)^2*m(1)^2+360*m(2)*m(1)^4-120*m(1)^6)/720;

    outC = cumsum(c,2);
    W = outC/Lambda;
end