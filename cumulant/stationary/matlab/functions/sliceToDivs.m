function [divs,cell_counter,is_completeTree] = siiceToDivs(slice)
    % 
    ndata = size(slice,1);
    currentD = 0;
    currentLineage = 1;

    % output
    divs = [];
    cell_counter = [];

    % count the number of division from slice data
    for idata = 1:ndata
        if  (idata == ndata) | (slice(idata) <= slice(idata+1)) 
            if currentD > 0
                cell_counter(1:slice(idata), currentLineage) = 2^(-currentD) .* ones(slice(idata),1);
            else
                cell_counter(1:51, currentLineage) = 2^(-currentD) .* ones(51,1);
            end
            divs(1,currentLineage) = currentD; 
            currentLineage = currentLineage + 1;
            currentD = 0;
        else
            if currentD == 0
                cell_counter(slice(idata+1):51,currentLineage) = 2^(-currentD) .* ones(51-slice(idata+1)+1,1);
            else
                cell_counter(slice(idata+1):slice(idata),currentLineage) = 2^(-currentD) .* ones(slice(idata)-slice(idata+1)+1,1);
            end
            currentD = currentD + 1;
        end
    end
    

    nLineage = size(divs,2);
    %disp(divs);
    % check the conssitency of the tree by the number of initial cells
    N0 = 0;
    is_completeTree = 1;
    for iLineage = 1:nLineage
        N0 = N0 + 2^(-divs(1,iLineage));
    end
    if abs(N0 - round(N0)) > 0.001
        is_completeTree = 0;
    end
end