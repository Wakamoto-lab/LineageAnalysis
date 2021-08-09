function [C,W,lambda,N0,Ntau,S1, S2, s1, s2, D, Qcl]=SScumulants(ndiv,T_tot, isboot)
    %% initialize the parameters
    [~,nlineage]=size(ndiv);    % the total number of lineages
    minD=min(ndiv);             % maximum number of divisions along lineages
    maxD=max(ndiv);             % minimum number of divisions along lineages
    D=minD:1:maxD;              % array of number of divisions from minD to maxD     
    [~,nbinD]=size(D);          % total bin number of the histogram
    Qcl=zeros(1,nbinD);         % chronological distribution of the number of divisions
    Qrs=zeros(1,nbinD);         % retrospective distribution of the number of divisions
    
    %% calculate the chronological distribution Qcl(D)
    % retrospective frequency Q'rs(D)
    for ibinD=1:nbinD
        for ilin = 1:nlineage
            if ndiv(1,ilin)==(ibinD-1+minD)
                Qrs(1,ibinD) = Qrs(1,ibinD) + 1;
            end
        end
    end
    
    % normalize the retrospective frequency
    Qrs = Qrs/sum(Qrs);
    
    % convert Qrs to Qcl
    for ibinD=1:nbinD
        Qcl(1,ibinD) = Qrs(1,ibinD) * exp(-D(1,ibinD)*log(2));
    end
    
    % calculate the population growth rate
    lambda = log(1/sum(Qcl))/T_tot;
    
    % normalize the chronological frequency
    Qcl = Qcl/sum(Qcl);

    %% plot the chronological frequency
    if isboot == 0
        figure(5);clf;
        plot(D,Qcl,'o');
        savefig('QclD');
    end
    
    %% calculate the 1st to 6th cumulants
    m=zeros(1,6);   % array for momments
    c=zeros(1,6);   % array for cumulants
    
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
    
    c = c/T_tot;
    C(1,:) = c;
    C(2,:) = cumsum(c,2);

    % calculate the 1st to 6th cumulant weight and cumulative wight
    w=c/lambda;
    W(1,:) = w;
    W(2,:) = cumsum(w,2);
    
    %% calculate the number of initial and final cells
    N0 = 0;
    for ilineage = 1:nlineage
        N0 = N0 + 2^(-ndiv(ilineage));
    end
    Ntau = nlineage; % Ntau is equal to nlineages
   
    % calculate selection strength
    S1 = 0;
    S2 = 0;
    for ibinD=1:nbinD
        if Qcl(1,ibinD) > 0
            S1 = S1 + Qcl(1,ibinD) *log(Qcl(1,ibinD)/Qrs(1,ibinD));
            S2 = S2 + Qrs(1,ibinD) *log(Qrs(1,ibinD)/Qcl(1,ibinD));
        end
    end
    s1 = S1/(T_tot*lambda);
    s2 = S2/(T_tot*lambda);
end