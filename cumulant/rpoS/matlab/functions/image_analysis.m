function image_analysis()
    micron_per_pixel=0.064;
    tint=input('enter the time interval of time lapse(min): ');
    disp_micpix=horzcat('micron per pixel is ',num2str(micron_per_pixel));
    disp(disp_micpix);
    disp('choose the mother directory');
    motherdir=uigetdir();
    choose_ana=disp_analysis_menu();
    while(choose_ana>=0)
        switch choose_ana
            case 0
                data_preparation_singlerep(motherdir,micron_per_pixel,tint);
            case 1
                data_preparation_dualrep(motherdir,micron_per_pixel,tint);
            case 2
                plot_dynamics_signle_phenotype(motherdir,tint);
            case 3
                plot_dynamics_double_phenotype(motherdir,tint);
            case 4
                onedim_statistics_fitness_selection(motherdir,tint);
            case 5
                twodim_statistics_fitness_selection(motherdir,tint);
            case 6
                combine_FL_for_conditions(motherdir);
            case 7
                combine_onedimSS(motherdir);
            case 8
                combine_twodimSS(motherdir);
            case 9
                auto_correlation_on_lineages(motherdir,tint);
            case 10
                combine_AC(motherdir);
            case 11
                evaluate_pvalues_for_conditions(motherdir);
        end
        choose_ana=disp_analysis_menu();
    end
end

function choose_ana=disp_analysis_menu()
    disp('0. data preparation for single reporter');
    disp('1. data preparation for dual reporter');
    disp('2. dynamics of single phenotype on lineages');
    disp('3. dynamics of a pair of phenotypes on lineages');
    disp('4. calculate the "onedim" selection strength and fitness landscape');
    disp('5. calculate the "twodim" selection strength and fitness landscape');
    disp('6. combine the fitness landscape for condtions');
    disp('7. combine the oendim selection strength');
    disp('8. combine the twodim selection strength');
    disp('9. calculate the auto correlation function for GFP and RFP');
    disp('10. combine the autocorrelation');
    disp('11. combine pvalue');
    choose_ana=input('choose analysis(to exit->-1): ');
end