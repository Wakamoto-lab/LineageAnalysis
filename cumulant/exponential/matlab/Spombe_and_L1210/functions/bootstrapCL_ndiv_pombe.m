function [cerr, werr, LwSDerr] = bootstrapCL_ndiv_pombe(ndiv, nInitial, nBootstrap, T_tot)
    %% initialize output arrays
    list_c = zeros(nBootstrap, 6); % list of bootstrap samples of 1st - 6th cumulants
    list_w = zeros(nBootstrap, 6); % list of bootstrap samples of 1st - 6th cumulant weights
    list_LwSD = zeros(nBootstrap, 5); % list of bootstrap samples of population growth rate and contribution of selection strength for population growth rate
    N0 = size(ndiv,2);
    %% main: Bootstrap method
    for iboot = 1:nBootstrap
        ndiv_boot = zeros(1,nInitial); % initialize bootstrap array
        for iInitial = 1:nInitial
            index = randi([1,N0]); % choose lineage randomly
            ndiv_boot(1,iInitial) = ndiv(1,index);
        end
        % calculate cumulants, weights and growth rate
        [C,W,lambda,~,S1,S2,s1,s2] = cumulant_to_6th_pombe(ndiv_boot,T_tot, 1);
        lwsd_boot = [lambda,S1,S2,s1,s2];

        % save the results for randomly generated list to list of bootstrap samples 
        list_c(iboot,:) = C(1,:);
        list_w(iboot,:) = W(1,:);
        list_LwSD(iboot,:) = lwsd_boot;
    end
    list_C = cumsum(list_c,2);
    list_W = cumsum(list_w,2);

    %% calcualte error bars
    cerr = zeros(2,6);
    werr = zeros(2,6);
    LwSDerr = zeros(1,5);
    for i=1:6
        cerr(1,i) = 2*std(list_c(:,i));
        cerr(2,i) = 2*std(list_C(:,i));
        werr(1,i) = 2*std(list_w(:,i));
        werr(2,i) = 2*std(list_W(:,i));
        if i < 6
           LwSDerr(1,i) = 2*std(list_LwSD(:,i)); 
        end
    end
end