function  [initial_index,totID] = get_initalcell_index(allROI,nROI,nlineage)
    %% initialize arrays
    cellID = zeros(1,nROI);     % cellID on each ROI
    totID = 0;                  % current number of cellID
    initial_index = zeros(1,nlineage); % initial_index(1,ilineage) = initial ID of the ancestor of ilineage-th lineage
    
    %% distribution of  ancestor cellID
    for ilineage = 1:nlineage
        if cellID(allROI(1,ilineage)) == 0 % If the ancestior ROI is not distributed its cellID
            totID = totID + 1; % count up cellID
            cellID(allROI(1,ilineage)) = totID;
        end
        % contain the ancestor cellID to initial_index
        initial_index(1,ilineage) = cellID(allROI(1,ilineage));
    end
end