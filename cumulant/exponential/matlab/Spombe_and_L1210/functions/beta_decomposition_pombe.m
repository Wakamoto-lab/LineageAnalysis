function beta_decomposition_pombe(D,Qcl,ave_div,figdir)
    % get the number of bin for D
    [~,nbinD] = size(D);
    
    % initialization 
    beta = 0:0.01:1;
    meanD_beta = zeros(1,101);

    %% calculate <Dln2>_beta
    for ibeta=1:101
        % calculate Q_beta 
        Q_beta = zeros(1,nbinD);
        for ibinD = 1:nbinD
            Q_beta(1,ibinD) = exp(beta(1,ibeta)*D(ibinD)*log(2))*Qcl(1,ibinD);
        end
        Q_beta = Q_beta/sum(Q_beta); % normalization
        meanD_beta(1,ibeta) = log(2)*D*Q_beta'; % average of Dln2 with respect to Q_beta
    end

    %% plot 
    cd(figdir);
    figure(5);clf;
    plot(beta,meanD_beta);
    hold on;
    % approx line by interpolation by beta = 0 and 1
    linear = meanD_beta(1,1) * ones(1,101) + (0:0.01:1)*(meanD_beta(1,101)- meanD_beta(1,1));
    plot(beta,linear,'k--');
    % approx line by 3rd cumulant
    second = meanD_beta(1,1) * ones(1,101) +  (0:0.01:1)*(D.*D*Qcl' - (D*Qcl')^2)/2;
    plot(beta,second,'r:');
    
    grid on;
    ylim([0,ceil(meanD_beta(1,101))+1]);
    savefig(horzcat('meanDln2_beta_',num2str(ave_div)));

    figure(6); clf;
    r = 1 - 2.^(beta-1);
    plot(r,meanD_beta);
    grid on;
    ylim([0,ceil(meanD_beta(1,101))+1]);
    savefig(horzcat('meanDln2_r_',num2str(ave_div)));
end