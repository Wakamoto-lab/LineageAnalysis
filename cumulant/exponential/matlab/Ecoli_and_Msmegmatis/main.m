function cumulant_analysis()
     
    %% initialize the directory
    disp('choose the directory containing results files for colony growth data: ');
    maindir=uigetdir();
    cd(maindir);

    %% initialize arrays of cumulants and weights and parameters
    C = zeros(2,6);     % cumulants and cumulative sums of cumulants up to 6th cumulant
    W = zeros(2,6);     % weights and cumulative weights up to 6th weight
    TNL = zeros(1,10);   % array for [T_tot, Lambda, N_0, N_tau, s_slice, e_slice, S1, S2, s1, s2] 
                            % T_tot: width of time windows (hr),
                            % Lambda: population growth rate (/hr),
                            % N_0 : initial number of cells
                            % N_tau : final number of cells
                            % s_slice : slice number of first slice of cropped data
                            % e_slice : slice number of last slice of cropped data
                            % S1 : S_{KL}^{(1)}[D]
                            % S2 : S_{KL}^{(2)}[D]
                            % S1/tauLambda : normalized selection strength 1, s1 = S1/tauLambda
                            % S2/tauLambda : normalized selection strength 2, s2 = S2/tauLambda
    % bootstrap error estimatations
    Cerr = zeros(2,6); % 2sd of bootstrap estimations for array c
    Werr = zeros(2,6); % 2sd of bootstrap estimations for array w
    LamwSDerr = zeros(1,5); % 2sd of bootstrap estimations for array [tauLambda, S1, S2, s1, s2]
    ndivList = [];

    %% creat lists of the numebr of divisions along lineages for each positions
    % creat directories to save the list of the number of divisions
    % and ROI information
    mkdir ndiv;
    mkdir ROI;
    
    % get time interval of time-laps measurement
    str_tint = extractAfter(extractBefore(maindir,'min'),'Tint');
    tint = str2num(str_tint);
    
    % creat list of the number of divisions
    [ndiv,nslice, s_slice, e_slice] = create_ndivfile(maindir,tint);    % ndiv: list of the numebr of divisions,                            
    ndivList = [ndivList, ndiv]; 

    % total time (/hr)
    T_tot = tint*(nslice-1)/60; 

    % calculate cumulants, weights and the population growth rate
    [C, W, lambda, N0, Ntau, S1, S2, s1, s2, D, Qcl] = SScumulants(ndivList,T_tot,0); 

    % save TNL
    LamwSD = [T_tot, N0, Ntau, s_slice, e_slice, lambda, S1, S2, s1, s2];

    %% plot beta vs. <Dln2>_beta
    Hfunction(D,Qcl);
    
    % bootstrap for estimating uncertainty of variables
    [cerr, werr, Cerr, Werr, LamwSDerr] = bootstrap_ndiv_Lineage(20000, ndivList, T_tot);
    LamwSDerr = [0,0,0,0,0, LamwSDerr(1,:); 0,0,0,0,0, LamwSDerr(2,:)];
    %% save the cumulants, weights and parameters in main/Ecoli or main/Msmegmatis
    cd(maindir);
    csvwrite('ndivList.csv',ndivList);
    csvwrite('Ttot_lambda.csv', [LamwSD;LamwSDerr]);
    csvwrite('cumulants.csv',[C(1,:);cerr]);
    csvwrite('cumulative_cumulants.csv',[C(2,:);Cerr]);
    csvwrite('weights.csv',[W(1,:);werr]);
    csvwrite('cumulative_weights.csv',[W(2,:);Werr]);
end