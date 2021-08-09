function [min_nslice, N]=growth_curve(tint,dname,ndata,dir)
    %% directories
    datadir = horzcat(dir,'/data');
    ROIdir = horzcat(dir,'/ROI');
    
    %% get nslice
    % initialize parameters
    min_nslice=100000;  % minimal number of total slices among Results files
    max_nslice=0;       % maximal number of total slices among Resutls files
    
    %% get maximum and minimum number of slice
    for idata=1:ndata
        if contains(dname(idata,:),'.csv') == 1
            % get filename
            filename = strtrim(dname(idata,:));
            
            if contains(filename,'Results') == 1
                cd(datadir);
                % data import
                Data = csvread(filename);
                
                % get se sequence of slice number
                slice = Data(:,5);
                
                % update the max_nslice and min_nslice
                max_nslice = max(max_nslice,max(slice));
                min_nslice = min(min_nslice,max(slice));
            end
        end
        clear Data slice filename;
    end
    
    %% initializing arrays for plotting growth curves
    N = ones(1,max_nslice);           % sequence of number of cells
    t = tint:tint:tint*max_nslice;  % sequence of time (min);
    npos = 0;                       % current number of positions
    
    %% main
    for idata=1:ndata
        if contains(dname(idata,:),'.csv') == 1
            filename = strtrim(dname(idata,:));
            if contains(filename,'Results') == 1
                % count up the number of positions
                cd(datadir);
                npos = npos+1;
                
                % data import
                Data = csvread(filename);
                
                % separate data into slice, lastindex, prevROI
                slice = Data(:,5);
                lastindex = Data(:,6);
                prevROI = Data(:,7);
                
                %% separate ROI information per lineages
                [allROI,alldivision,nROI,nlineage,Nslice,isFullTree]=separate_ROIs(prevROI,lastindex,slice);
                isFullTree = 1;
                if isFullTree == 1
                    %% save ROI information
                    cd(ROIdir);
                    savename=sprintf('ROI_pos%04d.csv',npos);
                    csvwrite(savename,allROI);

                    %% sequence of the number of cells
                    for islice = 1:Nslice
                        [~,ni] = size(unique(allROI(islice,:)));
                        N(1,islice) = N(1,islice) + ni;
                    end
                end
            end
        end
    end
    
    %% plot growth curve
    disp(horzcat('N0 = ', num2str(N(1))));
    figure(1);clf;
    plot(tint*(1:min_nslice),log2(N(1:min_nslice)));
    hold on;
    grid on;
end