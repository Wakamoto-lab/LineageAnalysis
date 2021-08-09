function data_preparation_dualrep()
    micron_per_pixel = 0.064;
    tint = 3;
    disp('choose the data folder "mmdd" ');
    motherdir = uigetdir('choose the data folder "mmdd" '); %% motherdir which contains data folders "mmdd"
    cd(motherdir);
    %% save directories
    mkdir lineage; 
    lineagedir=horzcat(motherdir,'/lineage');
    mkdir phenotype; 
    phenotypedir=horzcat(motherdir,'/phenotype');
    %make new directories for save the phenotype data
    cd(phenotypedir);
    mkdir list; 
    listdir=horzcat(phenotypedir, '/list');
    mkdir com; 
    comdir=horzcat(phenotypedir,'/com');
    cd(comdir);
    mkdir list; 
    com_listdir=horzcat(comdir,'/list');
    mkdir all; 
    com_alldir=horzcat(comdir,'/all');

    %% set total number of slices
    nslice = input('enter the total number of slices: ');
    %% set data from the start slice to the final slice
    s_slice = input('enter the start slice: ');
    e_slice = input('enter the final slice: ');
    trimed_slice = e_slice - s_slice + 1;
    t_slices = s_slice:e_slice;

    seg_index = input('segmented by YFP or RFP? YFP=1, RFP=2: ');
    
    %% import loop (loop for day)
    loop_index=1;day=0;
    cd(motherdir);

    while(loop_index>0)    
        totpos=0;
        num_dpl=zeros(2,1);
        %% day information
        day=day+1;
        
        %% data import
        %combine data
        datadir = combine_dualrep(motherdir,seg_index);
        %import filenames of combined data 
        [dname,ndata]=getFilename(datadir);
        
        %% analysis on each position
        for pos=1:ndata            
            pos_show = ['now analysing the ', num2str(pos), ' th position...'];
            disp(pos_show);
            totpos = totpos+1;
            % import idata th combined data
            filename = dname(pos,:);
            if contains(filename, '.csv');
                Data = data_autoimport(datadir,filename);
                % separate to phenotypic data
                area = Data(:,1)*power(micron_per_pixel,2);
                ch1mean = Data(:,2) - Data(:,8);
                xm = Data(:,3);
                ym = Data(:,4);
                slice = Data(:,5);
                lastindex = Data(:,6);
                prevROI = Data(:,7);
                ch2mean = Data(:,9) - Data(:,10);

                %linking and divisoin information
                nlineage = sum(lastindex);
                [nROI,~] = size(prevROI);
                is_linked = zeros(1,nROI);

                %% separate the data per lineages
                %arrays for save
                allarea = zeros(nslice,nlineage);
                allch1 = zeros(nslice,nlineage);
                allch2 = zeros(nslice,nlineage);
                allxm = zeros(nslice,nlineage);
                allym = zeros(nslice,nlineage);
                allROI = zeros(nslice,nlineage);
                alldivision = zeros(nslice,nlineage);
                nlin=0;
                cd(lineagedir);
                % separate into lineage
                for iROI = 1:nROI
                    lineagedata = zeros(nslice,7); %time/area/ch1/ch2/xm/ym/ROI
                    if lastindex(iROI) == 1
                        nlin = nlin+1;
                        %data separation to lieagedata
                        ROI = iROI;
                        for islice = 1:nslice
                            is_linked(1,ROI) = 1;
                            jslice = nslice+1-islice;
                            lineagedata(jslice,:) = [tint*slice(ROI),area(ROI),ch1mean(ROI),ch2mean(ROI),xm(ROI),ym(ROI),ROI];
                            ROI=prevROI(ROI);
                        end
                        %contain the phenotype data
                        allarea(:,nlin) = lineagedata(:,2);
                        allch1(:,nlin) = lineagedata(:,3);
                        allch2(:,nlin) = lineagedata(:,4);
                        allxm(:,nlin) = lineagedata(:,5);
                        allym(:,nlin) = lineagedata(:,6);
                        allROI(:,nlin) = lineagedata(:,7);
                    end
                    clear lineagedata;
                end

                %% division count per lineage
                is_divide = divisionROI(prevROI,is_linked);
                for ilin = 1:nlin
                    alldivision(:,ilin) = division_boolean(allROI(:,ilin),is_divide);
                end

                %% trim
                allarea = allarea(t_slices,:);
                allch1 = allch1(t_slices,:);
                allch2 = allch2(t_slices,:);;
                allxm = allxm(t_slices,:);;
                allym = allym(t_slices,:);;
                allROI = allROI(t_slices,:);;
                alldivision = alldivision(t_slices,:);

                %% data output per lineage
                lineage = zeros(trimed_slice,8);
                for ilin = 1:nlin
                    lineage(:,1) = (t_slices*tint)';
                    lineage(:,2) = allarea(:,ilin);
                    lineage(:,3) = allch1(:,ilin);
                    lineage(:,4) = allch2(:,ilin);
                    lineage(:,5) = allxm(:,ilin);
                    lineage(:,6) = allym(:,ilin);
                    lineage(:,7) = allROI(:,ilin);
                    lineage(:,8) = alldivision(:,ilin);
                    savename=sprintf('lineage_day%1d_pos%04d_%04d.csv',day,pos,ilin);
                    csvwrite(savename,lineage);
                end
                %% make phenotype list 
                pheno_list = zeros(3,nlin); %% D, ave_x and ave_y for nlineages 
                for ilin = 1:nlin
                    pheno_list(1,ilin) = sum(alldivision(:,ilin));
                    pheno_list(2,ilin) = mean(allch1(:,ilin));
                    pheno_list(3,ilin) = mean(allch2(:,ilin));
                end 

                %% dataoutput per phenotype
                cd(listdir);
                savename = sprintf('phenotypelist_day%1d_pos%04d.csv',day,pos);
                csvwrite(savename,pheno_list);
                
                %% number of lineages of position in the day
                num_dpl(:,totpos) = [pos;nlin];
            end
        end
        %% combine the dataset for position and save to com directories
        disp('combining for phenotype list...');
        combine_dataset_fordays(3,listdir,com_listdir,day);
        clear dname ndata;
        loop_index = input('any other data sets?: yes=1, no=0: ');
    end
    %% combine the dataset for several days and save to "all" directory
    disp('combining for phenotype list for all days...');
    combine_dataset(3, com_listdir,com_alldir);
end

function Data=data_autoimport(datadir,filename)
    readfile=horzcat(datadir, '/', filename);
    if contains(readfile,'.csv');
        Data=readmatrix(readfile);
    end
end

function div=division_boolean(ROI,is_divide)
    %% get the number of slices
    [nslice,~] = size(ROI);
    
    %% output 
    div = zeros(nslice,1);
    
    %% main
    for islice = 1:nslice
        if ROI(islice) == 0
            disp('there is empty roi in a lineage');
            disp(ROI);
        else
            div(islice,1)=is_divide(ROI(islice));
        end
    end
end

function is_div=divisionROI(prevROI,is_linked)
    %% initialzing
    %number of ROIs
    [nROI,~]=size(prevROI);
    
    %% output
    %initially 0 and is_div(*,1) = 1 if there exists a ROI with the same previous ROI
    is_div=zeros(1,nROI);
    
    %% main 
    for ROI=1:nROI
        % for the linked ROIs
        if is_linked(ROI) == 1
            % not for the initial cells 
            if prevROI(ROI) ~= 0
                % detect the previous ROI
                pROI = prevROI(ROI);
                for iROI = 1:nROI
                    %search the linked iROI whose pROI(iROI) = pROI
                    if prevROI(iROI) == pROI && iROI ~= ROI && is_linked(iROI) ==1
                        is_div(1,ROI) = is_div(1,ROI) + 1;
                    end
                end
            end
        end
    end
end

function savedir=combine_dualrep(resdir,seg_index)
    cd(resdir);
    disp('choose the directory to be analysed');
    dir=uigetdir('choose the directory to be analysed');
    cd(dir);
    mkdir com; 
    savedir=horzcat(dir,'/com');
    %% directory setting
    if seg_index==1
        segdir=horzcat(dir,'/YFP');
        intdir=horzcat(dir,'/RFP');
    else
        intdir=horzcat(dir,'/YFP');
        segdir=horzcat(dir,'/RFP');
    end
    
    %% import the filename 
    [dname,ndata]=getFilename(segdir);
    
    %% combine and save
    for idata=1:ndata
        filename = dname(idata,:);
        if contains(filename,'.csv')
            Data1=data_autoimport(segdir,filename);
            Data2=data_autoimport(intdir,filename);
            cd(savedir);
            Data=[Data1,Data2(:,2),Data2(:,8)];
            csvwrite(filename,Data);
        end
    end
end

function combine_dataset_fordays(nslice,datadir,comdir,nday)
    cd(datadir);
    [dname,ndata]=getFilename(datadir);
    savedata=zeros(nslice,1);
    cumlineage=0;
    for idata=1:ndata
        if contains(dname(idata,:),horzcat('day',num2str(nday))) == 1
            %data import for i th position
            data=data_autoimport(datadir,dname(idata,:));
            %count the number of lineages in i th position
            [~,nlineage]=size(data);
            %expand the savedata configuration
            savedata(nslice,cumlineage+nlineage)=0;
            slin=cumlineage+1;elin=cumlineage+nlineage;
            %transport data to savedata
            savedata(:,slin:elin)=data;
            %update the cumulative number of lineages
            cumlineage=cumlineage+nlineage;
        end
    end
    cd(comdir);
    savename=extractBefore(dname(1,:),'_day');
    savename=horzcat(savename, '_day',num2str(nday),'.csv');
    csvwrite(savename,savedata);
end

function combine_dataset(nslice,datadir,comdir)
    cd(datadir);
    disp(datadir);
    [dname,ndata]=getFilename(datadir);
    savedata=zeros(nslice,1);
    cumlineage=0;
    for idata=1:ndata
        %data import for i th position
        data=data_autoimport(datadir,dname(idata,:));
        %count the number of lineages in i th position
        [~,nlineage]=size(data);
        %expand the savedata configuration
        savedata(nslice,cumlineage+nlineage)=0;
        slin=cumlineage+1;elin=cumlineage+nlineage;
        %transport data to savedata
        savedata(:,slin:elin)=data;
        %update the cumulative number of lineages
        cumlineage=cumlineage+nlineage;
    end
    cd(comdir);
    savename=extractBefore(dname(1,:),'_day');
    savename=horzcat(savename,'.csv');
    csvwrite(savename,savedata);
end