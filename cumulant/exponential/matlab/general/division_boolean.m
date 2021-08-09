function div=division_boolean(ROI,is_divide)
    %% get the number of slices
    [nslice,~] = size(ROI);
    
    %% initialize output array: div(islice,1) = 1 if cell in the lineage divides at islice 
    div = zeros(nslice,1);
    
    %% main
    for islice = 1:nslice
        if ROI(islice) == 0 % if there is a slice for which any ROI is distributed 
            disp('there is empty roi in a lineage');
            disp(ROI);
        else
            div(islice,1)=is_divide(ROI(islice));
        end
    end
end