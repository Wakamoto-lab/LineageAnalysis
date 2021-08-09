function cumulant_for_pombe()
    set_default();
    %% get directory and filename
    disp('choose S.pombe directory');
    dir=uigetdir();
    cd(dir);
    [condname,ncond] = getFilename(dir);
    
    %% initialization of output arrays
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

    ave_div = 40;

    %% main
    for icond =1:ncond
        %% data import
        conddir = horzcat(dir,'/',strtrim(condname(icond,:)));
        cd(conddir);
        mkdir result; 
        result_dir = horzcat(conddir,'/result');
        disp(conddir);
        [dname,ndata] = getFilename(conddir);

        %disp(dname);
        for idata = 1:ndata
            cd(conddir);
            if contains(dname(idata,:),'Results')==1
                Data=dlmread(strtrim(dname(idata,:)),'\t',1,2);
                nROI = size(Data,1);

                % initialize paramters
                tint = 3; %time interval of time-laps
                prev_ndiv=0;

                %% find
                itotslice = 1;
                while itotslice > 0
                    totslice = 10 * itotslice + 1; 
                    ndiv=zeros(1,1);
                    nlineage=0;
                    %% make chronological ndiv list
                    for iROI = 1:nROI
                        if Data(iROI,5)==1 % for initial cell (slice = 1)
                            % tracking initializing
                            nslice=1;
                            div=0;
                            age = 0;
                            cROI = iROI; % current ROI
                            is_complete = 1;
                            while nslice < totslice
                                if Data(cROI,8) == 0
                                    if Data(cROI,10) == 0
                                        age = age +tint;
                                    elseif Data(cROI,10) == 1
                                        if (age < 20)
                                            %disp('segmentation error!');
                                            is_complete = 0;
                                            break;
                                        else
                                            div = div + 1; % add 1 when cell divides
                                            age = 0;
                                        end
                                    end
                                    cROI = Data(cROI,6); % update ROI
                                    nslice = Data(cROI,5); % count up the slice number
                                else
                                    is_complete = 0;
                                    break;
                                end
                            end
                            % for lineages which reach to the "totslice"
                            if is_complete == 1
                                if div > ave_div/2
                                    nlineage = nlineage + 1;
                                    ndiv(1,nlineage) = 0;
                                    ndiv(1,nlineage) = div;
                                end
                            end
                        end
                    end
                    prev_totslice = totslice;
                    prev_ndiv = ndiv;
                    itotslice = itotslice + 1;
                    %disp(horzcat('at ',num2str(itotslice - 1), 'hr, ','remained N0 = ', num2str(size(prev_ndiv,2))));
                    if mean(ndiv) > ave_div
                        break;
                    end
                end
                
                ndiv = prev_ndiv;
                totslice = prev_totslice;
                T_tot = tint*(totslice-1)/60; %(hr)

                %% calculate the extimated population growth and cumulants 
                mkdir figure; 
                figdir = horzcat(conddir,'/figure');
                [c,w,lambda,N0,S1,S2,s1,s2]=cumulant_to_6th_pombe(ndiv,T_tot,0,ave_div,figdir);
                m = [T_tot,N0,lambda,S1,S2,s1,s2];
                % output array
                W = w(1,:);
                cumW = w(2,:);
                C=c(1,:);
                cumC = c(2,:);
                M=m;

                %% bootstrap for error estimation
                % chronological sampling
                [cerr, werr, Cerr, Werr, Merr] = bootstrap_ndiv_pombe(ndiv, N0, 20000, T_tot);
                Merr = [0,0,Merr(1,:); 0,0, Merr(2,:)];


                clear Data w c m;

                %% output results
                cd(conddir);
                csvwrite(horzcat('Ttot_lambda_',num2str(ave_div),'.csv'),[M;Merr]);
                csvwrite(horzcat('cumulants_',num2str(ave_div),'.csv'),[C;cerr]);
                csvwrite(horzcat('weights_',num2str(ave_div),'.csv'),[W;werr]);
                csvwrite(horzcat('cumulative_cumulants_',num2str(ave_div),'.csv'),[cumC;Cerr]);
                csvwrite(horzcat('cumulative_weights_',num2str(ave_div),'.csv'),[cumW;Werr]);
            end
        end
        
    end
end