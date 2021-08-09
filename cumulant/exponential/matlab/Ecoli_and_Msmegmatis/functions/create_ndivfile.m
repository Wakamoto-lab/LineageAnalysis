function [all_ndiv,nslice, start_slice, end_slice]=create_ndivfile(dir,tint)
    %% direvtories to save list of number of divisions (ndiv) and ROI information (ROI) for each lineage
    ndivdir=horzcat(dir,'/ndiv');
    ROIdir=horzcat(dir,'/ROI');
    
    %% get filename of Results files
    datadir = horzcat(dir,'/data');
    cd(datadir);
    [dname,ndata]=getFilename(datadir);
    
    %% check the growth curve
    [min_nslice, N]=growth_curve(tint,dname,ndata,dir);
    
    %% input the slices to crop data during exponential growth
    disp(horzcat('nslice = ', num2str(min_nslice), ', tint = ', num2str(tint)));
    start_slice = input('enter the start slice for cropping: ');
    end_slice =input('enter the end slice for cropping: ');
    
    %% save the figure of growth curve
    % plot start_slice line (vertical line)
    plot([tint*start_slice tint*start_slice], [log2(N(1)) log2(max(N))],'k');
    % plot end_slice line (vertical line)
    plot([tint*end_slice tint*end_slice], [log2(N(1)) log2(max(N))],'k');
    % save figure
    cd(dir);
    savefig('growthcurve');
    
    %% make ndiv list for cropped time window
    make_cropped_ndivfile(start_slice,end_slice,ndivdir,ROIdir);
    
    %% combine the ndiv files
    % get ndiv file name
    cd(ndivdir);
    [Dname,Ndata]=getFilename(ndivdir);
    
    % array to save all ndiv
    all_ndiv=zeros(1,1);
    
    % combine ndiv files
    cumlineage=0; % cumulative number of  lineages
    for idata=1:Ndata
        if contains(Dname(idata,:),'.csv') == 1
            % data import
            data=csvread(strtrim(Dname(idata,:)));
            [~,nlineage]=size(data);
            
            % check the initial number of cells
            N0 = 0;
            for ilineage = 1:nlineage
                N0 = N0 + 2^(-data(1,ilineage));
            end
            N0int = round(N0);
            
            % combine ndiv if the estimated number of cells is integer
            if abs(N0int - N0) < 0.00001
                % set the start number of array and the end number
                slin=cumlineage+1;
                elin=cumlineage+nlineage;
                
                % expand the save array
                all_ndiv(1,elin)=0;
                
                %transport data to all_ndiv
                all_ndiv(1,slin:elin)=data;
                
                %update the cumulative number of lineages
                cumlineage=cumlineage+nlineage;
            end
        end
    end
    
    %% nslice
    nslice = end_slice - start_slice+ 1;
end