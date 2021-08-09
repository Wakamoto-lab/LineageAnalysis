function out = calculateSelectionStrength(ndiv)
    % get the number of lineage
    nLineage = size(ndiv,2);

    %calculate N0
    N0 = 0;
    for iLineage = 1:nLineage
        N0 = N0 + 2^(-ndiv(1,iLineage));
    end

    % calculate Pcl(sigma) and Prs(sigma)
    Pcl = zeros(1,nLineage);
    Prs = zeros(1,nLineage);
    for iLineage = 1:nLineage
        Pcl(1,iLineage) = 2^(-ndiv(1,iLineage))/N0;
        Pcl(1,iLineage) = 1/nLineage; 
    end


    % calculate <Dln2>cl
    Dln2_cl = 0;
    Dln2_rs=0;
    for iLineage = 1:nLineage
        Dln2_cl = Dln2_cl + ndiv(1,iLineage)*log(2) * Pcl(1,iLineage);
        Dln2_rs = Dln2_rs + ndiv(1,iLineage)*log(2) * Prs(1,iLineage);
    end


    % calculate tau*Lambda
    Lambda = 0;
    for iLineage = 1:nLineage;
        Lambda = Lambda + 2^(ndiv(1,iLineage)) * Pcl(1,iLineage);
    end

   % calculate S[D]
   S1 = Lambda - Dln2_cl;
   S2 = Dln2_rs - Lambda;
   
   % output
   out = [Lambda,S1];
end