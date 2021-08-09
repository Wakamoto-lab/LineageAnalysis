function evaluate_pvalues_for_conditions(motherdir)
    dir_pvalue = horzcat(motherdir,'/data/pvalue');
    % general notation of phenotype
    % p1 = D, p2 = c_GFP, p3 = c_mCherry, p4 = p_GFP, p5 = p_mCherry, p6 = lambda 

    %% 1.onedim selection strength S[X]
    % Row: S[p1]/S[p2]/S[p3]/S[p4]/S[p5]/S[p6]
    % Column: 1- glc37 vs. gly37, 2- gly37 vs. glc30, 3- glc30 vs. glc37 
    pS=zeros(3,6);
    % data import %
    S=zeros(9,6);
    for iconds = 1:3
        switch iconds
            case 1
                condname='/glc37';
            case 2
                condname='/gly37';
            case 3
                condname='/glc30';
        end
        dir=horzcat(motherdir,'/data',condname,'/onedim');
        cd(dir);
        s=csvread(horzcat(dir,'/pvalue_onedimSS.csv'));
        switch iconds
            case 1
                S(1:3,:) = s(1:3,:);
            case 2
                S(4:6,:) = s(1:3,:);
            case 3
                S(7:9,:) = s(1:3,:);
        end
    end
    
    % calculate p value %
    for ipheno = 1:6
        s1 = S(1:3,:); %glc37
        s2 = S(4:6,:); %gly37
        s3 = S(7:9,:); %glc30
        for iconds = 1:3
            switch iconds 
                case 1
                    disp(horzcat('mean = ', num2str(mean(s1(:,ipheno))), ', std = ', num2str(std(s1(:,ipheno)))));
                    [~, pS(iconds,ipheno)] = ttest2(s1(:,ipheno),s2(:,ipheno));
                case 2
                    [~, pS(iconds,ipheno)] = ttest2(s2(:,ipheno),s3(:,ipheno));
                case 3
                    [~, pS(iconds,ipheno)] = ttest2(s3(:,ipheno),s1(:,ipheno));    
            end
        end
    end
    
    % output pvalue%
    cd(dir_pvalue);
    savename=horzcat('pvalue_onedimSS.csv');
    csvwrite(savename,pS);
    clear s1 s2 s3;
    
    %% 2.twodim selection strength S[X,Y]
    % Row:
    % S[p2,p3]/S[p2,p4]/S[p2,p5]/S[p2,p6]/S[p3,p4]/S[p3,p5]/S[p3,p6]/S[p4,p5]/S[p4,p6]/S[p5,p6]
    % Column: 1- glc37 vs. gly37, 2- gly37 vs. glc30, 3- glc30 vs. glc37 
    pSXY=zeros(3,10);
    % data import %
    SXY=zeros(9,10);
    for iconds = 1:3
        switch iconds
            case 1
                condname='/glc37';
            case 2
                condname='/gly37';
            case 3
                condname='/glc30';
        end
        dir=horzcat(motherdir,'/data',condname,'/twodim/list');
        cd(dir);
        s=csvread(horzcat(dir,'/pvalue_twodimSS.csv'));
        switch iconds
            case 1
                SXY(1:3,:) = s(1:3,:);
            case 2
                SXY(4:6,:) = s(1:3,:);
            case 3
                SXY(7:9,:) = s(1:3,:);
        end
    end
    
    % calculate p value %
    for ipheno = 1:10
        s1 = SXY(1:3,:); %glc37
        s2 = SXY(4:6,:); %gly37
        s3 = SXY(7:9,:); %glc30
        for iconds = 1:3
            switch iconds 
                case 1
                    [~, pSXY(iconds,ipheno)] = ttest2(s1(:,ipheno),s2(:,ipheno));
                case 2
                    [~, pSXY(iconds,ipheno)] = ttest2(s2(:,ipheno),s3(:,ipheno));
                case 3
                    [~, pSXY(iconds,ipheno)] = ttest2(s3(:,ipheno),s1(:,ipheno));    
            end
        end
    end
    
    % output pvalue%
    cd(dir_pvalue);
    savename=horzcat('pvalue_twodimSS.csv');
    csvwrite(savename,pSXY);
    clear s1 s2 s3;
    
    %% 3.ij conditional selection strength S[X|Y]
    % Row:
    % S[p2|p3]/S[p2|p4]/S[p2|p5]/S[p2|p6]/S[p3|p4]/S[p3|p5]/S[p3|p6]/S[p4|p5]/S[p4|p6]/S[p5|p6]
    % Column: 1- glc37 vs. gly37, 2- gly37 vs. glc30, 3- glc30 vs. glc37 
    pcSij=zeros(3,10);
    % data import %
    cSij=zeros(9,10);
    for iconds = 1:3
        switch iconds
            case 1
                condname='/glc37';
            case 2
                condname='/gly37';
            case 3
                condname='/glc30';
        end
        dir=horzcat(motherdir,'/data',condname,'/twodim/list');
        cd(dir);
        s=csvread(horzcat(dir,'/pvalue_ij_conditional_selection_strength.csv'));      
        switch iconds
            case 1
                cSij(1:3,:) = s(1:3,:);
            case 2
                cSij(4:6,:) = s(1:3,:);
            case 3
                cSij(7:9,:) = s(1:3,:);
        end
    end
    
    % calculate p value %
    for ipheno = 1:10
        s1 = cSij(1:3,:); %glc37
        s2 = cSij(4:6,:); %gly37
        s3 = cSij(7:9,:); %glc30
        for iconds = 1:3
            switch iconds 
                case 1
                    [~, pcSij(iconds,ipheno)] = ttest2(s1(:,ipheno),s2(:,ipheno));
                case 2
                    [~, pcSij(iconds,ipheno)] = ttest2(s2(:,ipheno),s3(:,ipheno));
                case 3
                    [~, pcSij(iconds,ipheno)] = ttest2(s3(:,ipheno),s1(:,ipheno));    
            end
        end
    end
    
    % output pvalue%
    cd(dir_pvalue);
    savename=horzcat('pvalue_ij_conditional_SS.csv');
    csvwrite(savename,pcSij);
    clear s1 s2 s3;
    
    %% 4.ji conditional selection strength S[X|Y]
    % Row:
    % S[p3|p2]/S[p4|p2]/S[p5|p2]/S[p6|p2]/S[p4|p3]/S[p5|p3]/S[p6|p3]/S[p5|p4]/S[p6|p4]/S[p6|p5]
    % Column: 1- glc37 vs. gly37, 2- gly37 vs. glc30, 3- glc30 vs. glc37 
    pcSji=zeros(3,10);
    % data import %
    cSji=zeros(9,10);
    for iconds = 1:3
        switch iconds
            case 1
                condname='/glc37';
            case 2
                condname='/gly37';
            case 3
                condname='/glc30';
        end
        dir=horzcat(motherdir,'/data',condname,'/twodim/list');
        cd(dir);
        s=csvread(horzcat(dir,'/pvalue_ji_conditional_selection_strength.csv'));
        switch iconds
            case 1
                cSji(1:3,:) = s(1:3,:);
            case 2
                cSji(4:6,:) = s(1:3,:);
            case 3
                cSji(7:9,:) = s(1:3,:);
        end
    end
    
    % calculate p value %
    for ipheno = 1:10
        s1 = cSji(1:3,:); %glc37
        s2 = cSji(4:6,:); %gly37
        s3 = cSji(7:9,:); %glc30
        for iconds = 1:3
            switch iconds 
                case 1
                    [~, pcSji(iconds,ipheno)] = ttest2(s1(:,ipheno),s2(:,ipheno));
                case 2
                    [~, pcSji(iconds,ipheno)] = ttest2(s2(:,ipheno),s3(:,ipheno));
                case 3
                    [~, pcSji(iconds,ipheno)] = ttest2(s3(:,ipheno),s1(:,ipheno));    
            end
        end
    end
    
    % output pvalue%
    cd(dir_pvalue);
    savename=horzcat('pvalue_ji_conditional_SS.csv');
    csvwrite(savename,pcSji);
    clear s1 s2 s3;
    
    %% 5.PECO gamma[X,Y]
    % Row:
    % gamma[p2,p3]/gamma[p2,p4]/gamma[p2,p5]/gamma[p2,p6]/gamma[p3,p4]/gamma[p3,p5]/gamma[p3,p6]/gamma[p4,p5]/gamma[p4,p6]/gamma[p5,p6]
    % Column: 1- glc37 vs. gly37, 2- gly37 vs. glc30, 3- glc30 vs. glc37
    pgam=zeros(3,10);p2gam=zeros(3,10);
    % data import %
    gam=zeros(9,10);
    for iconds = 1:3
        switch iconds
            case 1
                condname='/glc37';
            case 2
                condname='/gly37';
            case 3
                condname='/glc30';
        end
        dir=horzcat(motherdir,'/data',condname,'/twodim/list');
        cd(dir);
        s=csvread(horzcat(dir,'/pvalue_gammaXY.csv'));
        switch iconds
            case 1
                gam(1:3,:) = s(1:3,:);
            case 2
                gam(4:6,:) = s(1:3,:);
            case 3
                gam(7:9,:) = s(1:3,:);
        end
    end
    
    % calculate p value %
    for ipheno = 1:10
        s1 = gam(1:3,:); %glc37
        s2 = gam(4:6,:); %gly37
        s3 = gam(7:9,:); %glc30
        for iconds = 1:3
            switch iconds 
                case 1
                    [~, pgam(iconds,ipheno)] = ttest2(s1(:,ipheno),s2(:,ipheno));
                case 2
                    [~, pgam(iconds,ipheno)] = ttest2(s2(:,ipheno),s3(:,ipheno));
                case 3
                    [~, pgam(iconds,ipheno)] = ttest2(s3(:,ipheno),s1(:,ipheno));    
            end
        end
    end
    
    % output pvalue%
    cd(dir_pvalue);
    savename=horzcat('pvalue_gammaXY_env.csv');
    csvwrite(savename,pgam);
    clear s1 s2 s3;
    
    % calculate p value %
    for ipheno = 1:10
        s1 = gam(1:3,:); %glc37
        s2 = gam(4:6,:); %gly37
        s3 = gam(7:9,:); %glc30
        for iconds = 1:3
            switch iconds 
                case 1
                    [~, p2gam(iconds,ipheno)] = ttest(s1(:,ipheno));
                case 2
                    [~, p2gam(iconds,ipheno)] = ttest(s2(:,ipheno));
                case 3
                    [~, p2gam(iconds,ipheno)] = ttest(s3(:,ipheno));   
            end
        end
    end
    
    % output pvalue%
    cd(dir_pvalue);
    savename=horzcat('pvalue_gammaXY_to0.csv');
    csvwrite(savename,p2gam);
    clear s1 s2 s3;
end