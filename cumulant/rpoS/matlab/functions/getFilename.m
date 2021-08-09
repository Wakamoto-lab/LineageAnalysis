function [dname,ndata]=getFilename(datadir)
    cd(datadir);
    D=dir;
    [ndata,~]=size(D);
    dname=strings(ndata,1);
    j=0;
    for i=1:ndata
        Name=D(i).name;
        if Name(1)~= '.'
            j=j+1;
            dname(j,:)=Name;
        end
    end
    dname=dname(1:j,:);
    dname=char(dname);
    ndata=j;
end