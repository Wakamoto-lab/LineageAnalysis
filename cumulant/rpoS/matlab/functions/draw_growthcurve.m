function draw_growthcurve()
    disp('choose the direcory to save growth curve');
    figdir=uigetdir('choose the direcory to save growth curve');
    for icond = 1:3
        nslice=input('enter the total number of slices: ');
        tint=3;
        t=tint:tint:tint*nslice;
        t=t/60;
        %% import loop (loop for day)
        loop_index=1;day=0;
        conddir=uigetdir('choose the condition directory');
        growthcurve=zeros(3,nslice);
        while(loop_index>0)    
            cd(conddir);
            %% day information
            day=day+1;
            growthcurve(day,nslice)=0;
            %% data import
            daydir=uigetdir('choose the day direcotry to be analysed');
            datadir=horzcat(daydir,'/YFP');
            %import filenames of combined data 
            [dname,ndata]=getFilename(datadir);
            %% analysis on each position
            for pos=1:ndata            
                pos_show=['now analysing the ', num2str(pos), ' th position...'];
                disp(pos_show);
                % import idata th combined data
                Data=data_autoimport(datadir,dname,pos);
                % separate to phenotypic data
                lastindex=Data(:,6);prevROI=Data(:,7);

                %linking and divisoin information
                nlineage=sum(lastindex);
                [nROI,~]=size(prevROI);
                is_linked=zeros(1,nROI);

                %% separate the data per lineages
                %arrays for save
                allROI=zeros(nslice,nlineage);
                nlin=0;
                % separate into lineage
                for iROI=1:nROI
                    lineagedata=zeros(nslice,1); %ROI
                    if lastindex(iROI)==1
                        nlin=nlin+1;
                        %data separation to lieagedata
                        ROI=iROI;
                        for islice=1:nslice
                            is_linked(1,ROI) = 1;
                            jslice = nslice+1-islice;
                            lineagedata(jslice,1) = ROI;
                            ROI=prevROI(ROI);
                        end
                        %contain the phenotype data
                        allROI(:,nlin) = lineagedata(:,1);
                    end
                    clear lineagedata;
                end
                for islice = 1:nslice
                    [~,ncell] = size(unique(allROI(islice,:)));
                    growthcurve(day,islice) = ncell;
                end
                figure(1);
                growthcurve(day,:) = growthcurve(day,:)/growthcurve(day,1);
                plot(t,growthcurve(day,:));hold on;
            end
            loop_index=input('any other day?: yes=1, no=0: ');
        end
        clf;
        %% draw mean-std growth curve
        % separately plot for each day
        figure(2);
        % mean and std plot
        mN=mean(growthcurve);stN=std(growthcurve);
        switch icond
            case 1
                errorbar(t,mN,stN,'b'); hold on;
            case 2
                errorbar(t,mN,stN,'r'); hold on;
            case 3
                errorbar(t,mN,stN,'y');
        end
        clear t mN stN growthcurve;
    end
    cd(figdir);
    savefig('growthcurve_meanstd');
    legend('glc37','gly37','glc30');
end