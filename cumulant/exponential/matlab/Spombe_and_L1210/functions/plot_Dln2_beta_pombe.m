function plot_Dln2_beta(dir)
    % go to "ndiv"
    ndivdir=horzcat(dir,'/ndiv');
    cd(ndivdir);
    
    % get ndiv file name
    [Dname,Ndata]=getFilename(ndivdir);
    
    % array to save all ndiv
    all_ndiv=zeros(1,1);
    
    % combine ndiv files
    cumlineage=0; % cumulative number of  lineages
    for idata=1:Ndata
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
    cd(dir);
    beta_decomposition(all_ndiv);
end
