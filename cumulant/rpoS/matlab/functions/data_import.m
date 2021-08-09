function Data=data_import()
    disp('choose the directory');
    dir=uigetdir;
    disp('choose the file');
    [filename, ~,~]=uigetfile({'*.csv','CSV files'}, 'choose the file', dir);
    readfile=horzcat(dir, '/', filename);
    Data=csvread(readfile);
end