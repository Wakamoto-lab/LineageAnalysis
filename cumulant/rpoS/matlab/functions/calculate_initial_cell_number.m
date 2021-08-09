function [N0s,dname] = calculate_initial_cell_number()
    dir=uigetdir();
    [dname,ndata]=getFilename(dir);
    N0s=zeros(1,ndata);
    for idata=1:ndata
        dfname=horzcat(dir,'/',dname(idata,:));
        division=csvread(dfname);
        ndiv=sum(division);
        [~,nlineage]=size(ndiv);
        for i=1:nlineage
            N0s(idata) = N0s(idata) + 2^(-ndiv(i));
        end
        clear division ndiv nlineage 
    end
end