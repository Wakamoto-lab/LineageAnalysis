function [allROI,alldivision,nROI,nlineage,nslice, isFullTree]=separate_ROIs(prevROI,lastindex,slice)
    %% initialization
    nlineage = sum(lastindex);        % total number of lineages
    nslice = max(slice);              % total number of slices
    [nROI,~] = size(prevROI);         % total number of ROIs
    allROI = zeros(nslice,nlineage);  % save array of ROIs
    nlin = 0;                         % current number of lineages
    isFullTree = 1;
    %% separate into lineages
    for iROI = 1:nROI
        % array to contain ROIs
        lineagedata = zeros(nslice,1);
        
        % backward tracking if the ROI is at the end of the lineage
        if lastindex(iROI) == 1
            % count up the current number of lineages
            nlin = nlin+1;
            
            % backward tracking
            ROI=iROI;
            while(ROI>0)
                lineagedata(slice(ROI,1),:) = ROI;
                % update ROI to the previous ROI
                ROI = prevROI(ROI);
            end
            
            % fulfill the unfilled element by the end ROI number  
            % if the length of the sequence of ROIs is not equal to the number of slice
            if slice(iROI) ~= nslice
                for jslice = (slice(iROI)+1):nslice
                    lineagedata(jslice,:) = lineagedata(slice(iROI),:);
                end
                isFullTree = 0;
            end
            
            % contain the phenotype data
            allROI(:,nlin) = lineagedata;
        end
        clear lineagedata;
    end
    
    %% complete the defocused slices and remove lineages which dose not reach to the initical slice
    for ilineage = 1:nlineage
        % counting down the number of lineages
        jlin = nlineage+1-ilineage;
        for islice = 2:(nslice-1)
            % fulfill element by the next fulfilled element for unfulfilled element
            if allROI(islice,jlin) == 0
                i = 0;
                while(allROI(islice,jlin) == 0)
                    i = i + 1;
                    if allROI(islice + i,jlin) ~= 0
                        allROI(islice,jlin) = allROI(islice + i,jlin);
                    end
                end
            end
        end
        
        % remove lineages with unfulfilled elements
        C = cumprod(allROI(:,jlin) > 0);
        if C(nslice,1) == 0
            allROI(:,jlin) =[];
        end
    end
    
    %% detect the linked ROIs
    % update the number of lineages
    [~,nlineage] = size(allROI);
    
    % array of linked ROI: linked = 1, unlinked = 0
    is_linked = zeros(1,nROI);
    for iROI = 1:nROI
        % is_linekd(1,iROI) = 1 if there are iROI in allROI
        if sum(sum(allROI == iROI)) >0
            is_linked(1,iROI) = 1;
        end
    end

    %% detect the newborn cells
    is_divide = divisionROI(prevROI,is_linked); % is_divide(1,ROI) = 1 if the ROI is a newborn cells
    alldivision = zeros(nslice,nlineage);       % array to save is_divide per lineages
    for ilin = 1:nlineage
        for islice = 1:nslice
            alldivision(islice,ilin) = is_divide(allROI(islice,ilin));
        end
    end
end