function char = phenotypename(i)
    %% filename according to phenotype   
    if i==1
        char = 'ndiv';
    elseif i==2
        char = 'gfp_mean';
    elseif i==3
        char = 'rfp_mean';
    elseif i==4
        char = 'gfp_prod';
    elseif i==5
        char = 'rfp_prod';
    else 
        char = 'elongation';
    end
end