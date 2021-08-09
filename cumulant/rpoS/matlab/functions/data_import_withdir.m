function Data=data_import_withdir(dir)
    cd(dir);
    [filename, ~,~]=uigetfile({'*.csv','CSV files'}, 'choose the file', dir);
    readfile=horzcat(dir, '/', filename);
    Data=csvread(readfile);
end