function auto_correlation_on_lineages(motherdir,tint)
    %% window width
    nslice=input('enter the total number of slices: ');
    delslice=nslice/2;
    lineagedir=horzcat(motherdir,'/lineage');
    %% output
    PhiG =zeros(10,delslice);%[Phicl_day1;Phicl_day2;Phicl_day3;Phirs_day1;Phirs_day2;Phirs_day3;Phicl_mean;Phicl_std;Phirs_mean;Phirs_std;]
    PhiR =zeros(10,delslice);%[Phicl_day1;Phicl_day2;Phicl_day3;Phirs_day1;Phirs_day2;Phirs_day3;Phicl_mean;Phicl_std;Phirs_mean;Phirs_std;]
    for day=1:3
        d=horzcat('day',num2str(day));
        disp(d);
        %% import num_dpl file
        dplprint=sprintf('num_dpl_day%01d.csv',day);
        readfile=horzcat(motherdir, '/phenotype/',dplprint);
        num_dpl=csvread(readfile);
        %% data import
        disp('data importing...');
        cd(lineagedir);
        totlin=0;
        [~,npos]=size(num_dpl);
        Nlin=sum(num_dpl(2,:));
        ch1=zeros(nslice,Nlin);ch2=zeros(nslice,Nlin);
        ROI=zeros(nslice,Nlin);division=zeros(nslice,Nlin);
        for ipos=1:npos
            pos=num_dpl(1,ipos);nlin=num_dpl(2,ipos);
            for ilin = 1:nlin
                lineageid=sprintf('lineage_day%1d_pos%04d_%04d.csv',day,pos,ilin);
                readfile=horzcat(motherdir, '/lineage/',lineageid);
                ldata=csvread(readfile);
                ch1(:,totlin+ilin)=ldata(:,3);
                ch2(:,totlin+ilin)=ldata(:,4);
                ROI(:,totlin+ilin)=ldata(:,7);
                division(:,totlin+ilin)=ldata(:,8);
            end
            totlin = totlin+nlin;
        end
        disp('ending data import...');
        %% calculate the auto correlation fucntion
        disp('calculate the auto correlation funciton');
        for m = 1:delslice
            d=horzcat('now auto correlation window is ', num2str(m*tint), ' min');
            disp(d);
            for ista=1:2 % 1sta=1(cl), ista=2(rs) 
                % initialize the parts of auto correlation
                EGcc=0;EGc_s=0;EGc_e=0;EGc2_s=0;EGc2_e=0; % GFP
                ERcc=0;ERc_s=0;ERc_e=0;ERc2_s=0;ERc2_e=0; % GFP
                W=0; % [cl;rs]
                for endslice=(m+1):nslice
                    %% searching the endcell at endslice
                    totlin=0;
                    for ipos = 1:npos
                        endROI=0;nlin=num_dpl(2,ipos);
                        for ilin=1:nlin
                            if ismember(endROI,ROI(endslice,totlin+ilin)) == 0
                                endROI = [endROI, ROI(endslice,totlin+ilin)];
                                % weight
                                if ista==1
                                    w=2^(-sum(division((endslice-m):endslice,totlin+ilin)));
                                else
                                    w=1;
                                end
                                W = W + w;
                                % E and V
                                % GFP
                                Gc_s=ch1((endslice-m),totlin+ilin);Gc_e=ch1(endslice,totlin+ilin);
                                EGc_s = EGc_s + w*Gc_s;
                                EGc_e = EGc_e + w*Gc_e;
                                EGcc = EGcc + w*Gc_s*Gc_e; 
                                EGc2_s = EGc2_s + w*Gc_s*Gc_s;
                                EGc2_e = EGc2_e + w*Gc_e*Gc_e;
                                % RFP
                                Rc_s=ch1((endslice-m),totlin+ilin);Rc_e=ch1(endslice,totlin+ilin);
                                ERc_s = ERc_s + w*Rc_s;
                                ERc_e = ERc_e + w*Rc_e;
                                ERcc = ERcc + w*Rc_s*Rc_e; 
                                ERc2_s = ERc2_s + w*Rc_s*Rc_s;
                                ERc2_e = ERc2_e + w*Rc_e*Rc_e;
                            end
                        end
                        clear endROI;
                    end
                end
                %GFP
                VG_s=(EGc2_s/W) - (EGc_s/W)^2;
                VG_e=(EGc2_e/W) - (EGc_e/W)^2;
                PhimG=((EGcc/W)-(EGc_s/W)*(EGc_e/W))/sqrt(VG_s*VG_e);
                %RFP
                VR_s=(ERc2_s/W) - (ERc_s/W)^2;
                VR_e=(ERc2_e/W) - (ERc_e/W)^2;
                PhimR=((EGcc/W)-(ERc_s/W)*(ERc_e/W))/sqrt(VR_s*VR_e);
                % auto correlation
                PhiG(2*(day-1)+ista,m)=PhimG;
                PhiR(2*(day-1)+ista,m)=PhimR;
            end
        end
    end
    %% mean and std of auto correltaion
    for m = 1:delslice
        PhiG(7,m)=mean(PhiG(1:3,m));PhiG(8,m)=std(PhiG(1:3,m));PhiG(9,m)=mean(PhiG(4:6,m));PhiG(10,m)=std(PhiG(4:6,m));
        PhiR(7,m)=mean(PhiR(1:3,m));PhiR(8,m)=std(PhiR(1:3,m));PhiR(9,m)=mean(PhiR(4:6,m));PhiR(10,m)=std(PhiR(4:6,m));
    end
    %% plot auto correlations
    t =tint*(1:delslice);
    subplot(2,1,1);
    errorbar(t,PhiG(7,:),PhiG(8,:),'b');hold on;
    errorbar(t,PhiG(9,:),PhiG(10,:),'r');
    subplot(2,1,2);
    errorbar(t,PhiR(7,:),PhiR(8,:),'b');hold on;
    errorbar(t,PhiR(9,:),PhiR(10,:),'r');
    
    %% save auto correlation
    phenotypedir=horzcat(motherdir,'/phenotype');
    cd(phenotypedir);
    savename = sprintf('auto_correlation_GFP.csv');
    csvwrite(savename,PhiG);
    savename = sprintf('auto_correlation_RFP.csv');
    csvwrite(savename,PhiR);
end