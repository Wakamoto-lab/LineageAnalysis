function is_div=divisionROI(prevROI,is_linked)
    %% initialzing
    %number of ROIs
    [nROI,~]=size(prevROI);
    
    %% output
    %initialize is_div array: is_div(*,1) = 1 if there exists a ROI with the same previous ROI
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
                        is_div(1,ROI) =  1;
                    end
                end
            end
        end
    end
end