function make_cropped_ndivfile(start_slice,end_slice,ndivdir,ROIdir)
    %% make croppedROI directory under the ROI directory
    cd(ROIdir);
    mkdir croppedROI; cROIdir = horzcat(ROIdir,'/croppedROI');
    
    %% get filename of ROI
    [dname,ndata]=getFilename(ROIdir);
    
    %% crop ROI from start_slice to end_slice
    pos=0; % the current position
    for idata=1:ndata
        if contains(dname(idata,:),'ROI_')==1
            cd(ROIdir);
            pos = pos+1;% count up the current position
            
            
            % data import
            allROI=csvread(dname(idata,:));
            [nslice,nlineage]=size(allROI); % get the number of slices and lineages
            nROI =max(max(allROI));         % the maximum number of ROI 
            prevROI=zeros(nROI,1); % initialize the array of previous ROI
            
            % resortate the prevROI
            for ilineage=1:nlineage
                for islice=2:nslice
                    prevROI(allROI(islice,ilineage))= allROI(islice-1,ilineage);
                end
            end
            
            % detect the linked ROIs: 
            is_linked=zeros(1,nROI); % initialize is_linked array
            for iROI = 1:nROI
                % If iROI exists in allROI, is_linked(1,iROI) = 1,
                % else is_linked(1,iROI) = 0.
                if sum(sum(allROI==iROI)) >0
                    is_linked(1,iROI)=1;
                end
            end
            
            % crop ROI data
            croppedROI = allROI(start_slice:end_slice,:); % restore ROI information for the cropped time-window
            [cnslice,~] = size(croppedROI);     % the total number of slices of cropped data

            % remove overlapping lineages
            is_only = ones(1,nlineage); % is_only(1,ilineage) = 1 if the lineage does not have any overlapping lineages
            for ilineage = 1:(nlineage-1)
                if is_only(1,ilineage) == 1
                    for jlineage = (ilineage+1):nlineage
                        if isequal(croppedROI(:,ilineage),croppedROI(:,jlineage)) ==1
                            is_only(1,jlineage)=0;
                        end
                    end
                end
            end
            
            % count the un-overlapped lineages
            cnlineage=sum(is_only);
            
            % restore the cropped ROI sequences 
            cROI=zeros(cnslice,cnlineage);
            cnlin = 0; % current numbr of lineage in cropped data
            for clin = 1:nlineage
                if is_only(1,clin) == 1
                    % count up the current number of lineage
                    cnlin = cnlin + 1;
                    % copy the cropped ROI sequence
                    cROI(:,cnlin) = croppedROI(:,clin);
                end
            end
            
            %% division count per lineage
            
            is_divide=divisionROI(prevROI,is_linked);
            
            %% make the list of the number of divisions for each lineage 
            alldivision=zeros(cnslice,cnlineage); % initialize alldivision array
            ndiv = zeros(1,cnlineage); % initialize ndiv array
            for ilin=1:cnlineage
                alldivision(:,ilin) = division_boolean(cROI(:,ilin),is_divide); % If cell in ilin-th lineage divides at islice, alldivision(islice,ilin) = 1
                ndiv(1,ilin) = sum(alldivision(:,ilin)); % the number of divisions along ilin-th lineage
            end
            
            %% get initical-cell index and separate ndiv data per a lineage tree
            [init_ID,maxID] = get_initalcell_index(cROI,nROI,cnlineage);
            
            %% data output per single lineage tree
            % ROI
            cd(cROIdir);
            savename=sprintf('croppedROI_pos%04d.csv',pos);
            csvwrite(savename,cROI);            
            % ndiv list per single lineage tree
            cd(ndivdir);
            for iID = 1:maxID
                ndiv_singleLT = separate_ndiv(ndiv,init_ID,iID); % extract the number of divisions for lineages whose ancestor is iID cell
                if size(ndiv_singleLT,2) > 1
                    savename=sprintf('ndiv_pos%04d_%01d.csv',pos,iID);
                    csvwrite(savename,ndiv_singleLT); 
                end
            end
            clear allROI prevROI is_linked croppedROI cnslice is_only cnlineage cROI alldivision is_divide alldivision pheno_list ndiv_singleLT;
        end
    end
end