function [M,C,W]=cumulant_for_L1210()
    %% initialization of output arrays
    ncond = 1;
    ave_div = 1210;
    C = zeros(ncond,6);
    cumC = zeros(ncond,6);
    W = zeros(ncond,6);
    cumW = zeros(ncond,6);
    M = zeros(ncond,7);

    Cerr = zeros(ncond,6);
    Werr = zeros(ncond,6);
    cumCerr = zeros(ncond,6);
    cumWerr = zeros(ncond,6);
    Merr = zeros(ncond,5);

    % time parameters 
    tint = 10; % time interval of the time-laps experiment
    nframe = 361; % number of frame to analyze
    T = (nframe-1)*tint/60; % time window of cumulant analysis, # of frames = 600;
    disp('choose the directory which contains data');
    dir = uigetdir();
    cd(dir);
    mkdir data;
    datadir = horzcat(dir,'/data');
    mkdir figure;
    figdir = horzcat(dir, '/figure');
    cd(datadir);
    [dname,ndata] = getFilename(datadir);
    
    % get ndivList
    nlin = 0;
    ndiv = zeros(1,1);
    for idata = 1:ndata
        data = csvread(strtrim(dname(idata,:)));
        if max(data(:,1)) > nframe
            nlin = nlin + 1;
            ngen = size(data,1);
            for igen = 1:ngen
                if data(igen,1) > nframe
                    div = igen;
                    break;
                end
            end
            ndiv(1,nlin)=div;
        end
        clear data;                                                                                                                                                                                               
    end
    
    %% calculate the extimated population growth and cumulants 
    icond = 1;
    [c,w,lambda,N0,S1,S2,s1,s2]=cumulant_to_6th_pombe(ndiv,T,0,1210,figdir);
    m = [T,N0,lambda,S1,S2,s1,s2];
    % output array
    W = w(1,:);
    cumW = w(2,:);
    C=c(1,:);
    cumC = c(2,:);
    M=m;

    %% bootstrap for error estimation
    % chronological sampling
    [cerr, werr, Cerr, Werr, Merr] = bootstrap_ndiv_pombe(ndiv, N0, 20000, T);
    Merr = [0,0,Merr(1,:); 0,0, Merr(2,:)];


    clear Data w c m;
 
    %% output results
    cd(dir);
    csvwrite(horzcat('Ttot_lambda_',num2str(ave_div),'.csv'),[M;Merr]);
    csvwrite(horzcat('cumulants_',num2str(ave_div),'.csv'),[C;cerr]);
    csvwrite(horzcat('weights_',num2str(ave_div),'.csv'),[W;werr]);
    csvwrite(horzcat('cumulative_cumulants_',num2str(ave_div),'.csv'),[cumC;Cerr]);
    csvwrite(horzcat('cumulative_weights_',num2str(ave_div),'.csv'),[cumW;Werr]);
end