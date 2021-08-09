function [C,W,lambda,N0,S1,S2,s1,s2]=cumulant_to_6th_pombe(ndiv,T_tot, isboot,ave_div,figdir)
    %% initialize the parameters
    [~,N0]=size(ndiv);    % the total number of lineages
    minD=min(ndiv);             % maximum number of divisions along lineages
    maxD=max(ndiv);             % minimum number of divisions along lineages
    D=minD:1:maxD;              % array of number of divisions       
    [~,nbinD]=size(D);          % total bin number of the histogram
    Qcl=zeros(1,nbinD);         % chronological frequency of the number of divisions
    Qrs=zeros(1,nbinD);         % retrospective frequency of the number of divisions

    %% calculate the retrospective distribution and population growth rate
    exp_lambda = 0;
    for ibinD=1:nbinD
        for ilin = 1:N0
            if ndiv(1,ilin)==(ibinD-1+minD)
                Qcl(1,ibinD) = Qcl(1,ibinD) + 1;
                exp_lambda = exp_lambda + 2^ndiv(1,ilin)/N0;
            end
        end
    end
    lambda = log(exp_lambda)/T_tot;
    
    %% normalize the chronological frequency
    Qcl = Qcl/sum(Qcl);
    
    %% calculate the restospective distribution
    for ibinD=1:nbinD
        Qrs(1,ibinD) = 2^ibinD * Qcl(1,ibinD);
    end
    %% normalize the retrospective frequency
    Qrs = Qrs/sum(Qrs);

    %% plot the chronological frequency
    if isboot == 0
        cd(figdir);
        figure(4);clf;
        plot(D,Qcl,'o');
        savefig(horzcat('QclD_',num2str(ave_div)));
        % plot <Dln2>_beta
        beta_decomposition_pombe(D,Qcl,ave_div,figdir)
    end
    
    %% calculate the 1st to 6th cumulants
    m=zeros(1,6);   % array for momments
    c=zeros(1,6);   % array for cumulants
    %output
    C=zeros(2,6);   % array for cumulants
    W=zeros(2,6);   % array for cumulants
    
    % calculate moments
    for order = 1:6
        for ibinD=1:nbinD
            m(order) = m(order) + ((D(ibinD)*log(2))^order) * Qcl(ibinD);
        end
    end
    
    % calculate cumulants
    c(1) = m(1);
    c(2) = (m(2) - m(1)^2)/2;
    c(3) = (m(3) - 3*m(2)*m(1) + 2*m(1)^3)/6;
    c(4) = (m(4) - 4*m(3)*m(1) -3*m(2)^2 + 12*m(2)*m(1)^2 - 6*m(1)^4)/24;
    c(5) = (m(5)-5*m(4)*m(1)-5*m(3)*m(2)+15*m(3)*m(1)^2+15*m(2)^2*m(1)-35*m(2)*m(1)^3+14*m(1)^5)/120;
    c(6) = (m(6)-6*m(5)*m(1)-15*m(4)*m(2)+30*m(4)*m(1)^2-10*m(3)^2+120*m(3)*m(2)*m(1)-120*m(3)*m(1)^3+30*m(2)^3-270*m(2)^2*m(1)^2+360*m(2)*m(1)^4-120*m(1)^6)/720;
    c=c/T_tot;
    
    C(1,:) = c;
    C(2,:) = cumsum(c,2);
    %% calculate the 1st to 6th cumulant weight
    w=c/lambda;
    W(1,:) = w;
    W(2,:) = cumsum(w,2);

    %% calculate selection strength
    S1 = 0;
    S2 = 0;
    for ibinD = 1:nbinD
        if Qcl(ibinD) > 0
            S1 = S1 + Qcl(ibinD) * log(Qcl(ibinD)/Qrs(ibinD));
            S2 = S2 + Qrs(ibinD) * log(Qrs(ibinD)/Qcl(ibinD));
        end
    end
    s1 = S1/(T_tot*lambda);
    s2 = S2/(T_tot*lambda);

end