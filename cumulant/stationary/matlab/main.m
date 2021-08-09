function divisionCounter()
    %% initialize
    ndivList = [];
    cellCounterList = [];
    
    %% choose directory containting Results.csv
    results_dir = uigetdir();

    %% count the number of Results file
    [dname,ndata] = getFilename(results_dir);    % typename : list of name of cellular strains
    
    for idata = 1:ndata
        if contains(dname(idata,:),'Results') == 1
            Data = dlmread(strtrim(dname(idata,:)),',',1,1);
            slice = Data(:,3);
            [divs, cell_counter, complete_tree] = sliceToDivs(slice);
            if complete_tree == 1
                ndivList = [ndivList, divs];
                cellCounterList = [cellCounterList, cell_counter];
            else
               disp(horzcat('There are some errors on ', dname(idata,:))); 
            end
            divs = [];
        end
    end
    %% save the ndiv list
    mkdir analysis;
    ana_dir = horzcat(results_dir,'/analysis');
    cd(ana_dir);

    % total time (/hr)
    T_tot = 3*(51-1)/60; 

    % calculate cumulants, weights and the population growth rate
    [C, W, lambda, N0, Ntau, S1, S2, s1, s2, D, Qcl] = SScumulants(ndivList,T_tot,0); 

    % save TNL
    TNL(1,:) = [T_tot, lambda, N0, Ntau, 1,51, S1, S2, s1, s2];

    %% plot beta vs. <Dln2>_beta
    figure(1);
    %Hfunction(D,Qcl);

    %% plot growth curve
    figure(2);
    grwothCurve = sum(cellCounterList,2);
    plot(1:51, grwothCurve);
    
    % bootstrap for estimating uncertainty of variables
    [Cerr, Werr, LamwSDerr] = bootstrap_ndiv_Lineage(20000, ndivList, T_tot);
    
    %% save the cumulants, weights and parameters in main/Ecoli or main/Msmegmatis
    cd(ana_dir);
    csvwrite('ndivList.csv',ndivList);
    csvwrite('cellCounterList.csv',cellCounterList);
    csvwrite('growthCurve.csv',grwothCurve);
    csvwrite('Ttot_lambda.csv',TNL);
    csvwrite('Qcl.csv',[D;Qcl]);
    csvwrite('cumulants.csv',C);
    csvwrite('weights.csv',W);
    csvwrite('err_Ttot.csv',LamwSDerr);
    csvwrite('err_cumulants.csv',Cerr);
    csvwrite('err_weights.csv',Werr);
end