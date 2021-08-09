function Data=data_autoimport(datadir,dname,idata)
    filename=dname(idata,:);
    readfile=horzcat(datadir, '/', filename);
    Data=csvread(readfile);
end