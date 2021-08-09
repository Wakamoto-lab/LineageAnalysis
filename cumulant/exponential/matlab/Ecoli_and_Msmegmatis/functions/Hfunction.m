function Hfunction(D,Qcl)
    % get max and min D
    [~,nbinD] = size(D);
    
    %% plot <Dln2>_beta
    meanD_beta = zeros(1,101);
    beta = 0:0.01:1; % beta in [0,1] with bin = 0.01
    for ibeta=1:101
        Q_beta = zeros(1,nbinD);
        % calculate Q_beta(D)
        for ibinD = 1:nbinD
            Q_beta(1,ibinD) = exp(beta(1,ibeta)*D(1,ibinD)*log(2))*Qcl(1,ibinD);
        end
        Q_beta = Q_beta/sum(Q_beta);
        % calculate <Dln2>_beta
        meanD_beta(1,ibeta) = log(2)*D*Q_beta';
    end
    % plot beta vs. <Dln2>_beta
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
    savefig('Hfunction');
end