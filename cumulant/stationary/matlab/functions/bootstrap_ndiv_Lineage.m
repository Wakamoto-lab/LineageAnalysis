function [Cerr, Werr, LwSDerr] = bootstrap_ndiv_Lineage(nBootstrap, ndivList, T_tot)
    %% initialize output arrays
    list_c = zeros(nBootstrap, 6); % list of bootstrap samples of 1st - 6th cumulants 
    list_C = zeros(nBootstrap, 6); % list of bootstrap samples of 1st - 6th cumulative cumulants
    list_w = zeros(nBootstrap, 6); % list of bootstrap samples of 1st - 6th cumulant weights
    list_W = zeros(nBootstrap, 6); % list of bootstrap samples of 1st - 6th cumulative weights
    list_LwSD = zeros(nBootstrap, 5); % list of bootstrap samples of [tauLambda, S1, S2, s1, s2]
    
    nLineage = size(ndivList,2); % number of lineage

    % calculate random sampling function from ndivList
    lineageWeightCL = cumsum(2.^(-ndivList),2)/sum(2.^(-ndivList));

    %% main: Bootstrap method
    alart  = 10; % show alart per "alart"% progress
    
    for iboot = 1:nBootstrap
        % make randomly generated list of number of divisions
        ndiv_boot = zeros(1,nLineage); % array for randomly generated list of number of divisions
        for ilin = 1:nLineage
            index = randi(nLineage);
            ndiv_boot(1,ilin) = ndivList(1,index); 
        end
        % calculate cumulants, cumulant weights and population growth rate
        % for randomly generated list of number of divisions
        [C,W,lambda,~,~,S1, S2, s1, s2, ~, ~] = SScumulants(ndiv_boot,T_tot, 1);
        lwsd_boot = [lambda, S1, S2, s1, s2];
        % save the results for randomly generated list to list of bootstrap
        % samples
        list_c(iboot,:) = C(1,:);
        list_C(iboot,:) = C(2,:);
        list_w(iboot,:) = W(1,:);
        list_W(iboot,:) = W(2,:);
        list_LwSD(iboot,:) = lwsd_boot;
        % show the progress alart 
        if alart/100 <= iboot/nBootstrap
           disp(horzcat(num2str(alart), '% of bootstrap has done...'));
           alart = alart + 10;
        end
    end
    
    %% calcualte error bars
    % initialize output arrays
    Cerr = zeros(2,6); % list of standard deviations of 1st - 6th cumulants
    Werr = zeros(2,6); % list of standard deviaitons of 1st - 6th cumulant weights
    LwSDerr = zeros(1,5); % standard deviations of population growth rate and contribution of selection strength for population growth rate
    % save the standard deviation to the output arrays
    for i=1:6
       Cerr(1,i) = 2*std(list_c(:,i));
       Cerr(2,i) = 2*std(list_C(:,i));
       Werr(1,i) = 2*std(list_w(:,i));
       Werr(2,i) = 2*std(list_W(:,i));
       if i < 6
          LwSDerr(1,i) = 2*std(list_LwSD(:,i)); 
       end
    end
end