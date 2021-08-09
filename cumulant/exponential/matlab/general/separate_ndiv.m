function ndiv_singleLT = separate_ndiv(ndiv,init_ID,iID)
    %% initialization
    ndiv_singleLT = zeros(1,1); % list of the number of divisions for lineages in the iID-th tree
    nlineage = size(ndiv,2);
    nlin = 0; % counter of lineages in iID-th tree

    for ilineage = 1:nlineage
        if init_ID(1,ilineage) == iID % if ilineage is in the iID=th lineage
            nlin = nlin + 1; % count up the number of lineage in the iID-th lineage
            ndiv_singleLT(1,nlin) = ndiv(ilineage);
        end
    end
end