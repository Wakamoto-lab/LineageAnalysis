function onedim_statistics_fitness_selection(nslice)
    tint = 3;
    motherdir = uigetdir('choose the folder contains data folders "mmdd" '); %% motherdir which contains data folders "mmdd"

    t_tot = tint*(nslice-1)/60; %% total measurement time
    %% data import
    alldir = horzcat(motherdir,'/phenotype/com/all');
    cd(alldir);
    disp('choose the "phenotype list" file');
    data = data_import_withdir(alldir);
    %% onedim analysis
    for i=1:3
        char_phenotype = phenotypename(i);
        disp(horzcat('now analying for ', char_phenotype));
        [FLx,M,C,W] = onedim_calculus(data,motherdir,t_tot,i);
        savename=horzcat('QclFLx_', char_phenotype ,'.csv');
        csvwrite(savename,FLx);
        savename=horzcat('M_', char_phenotype ,'.csv');
        csvwrite(savename,M);
        savename=horzcat('cumulant_', char_phenotype ,'.csv');
        csvwrite(savename,C);
        savename=horzcat('weight_', char_phenotype ,'.csv');
        csvwrite(savename,W);
    end
end