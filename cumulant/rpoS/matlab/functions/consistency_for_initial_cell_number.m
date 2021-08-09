function consistency_for_initial_cell_number(nslice)
    disp('choose the result file');
    Data=data_import();
    area=Data(:,1)*power(0.1,2);ch1mean=Data(:,2)-Data(:,8);xm=Data(:,3);ym=Data(:,4);slice=Data(:,5);lastindex=Data(:,6);prevROI=Data(:,7);
    %linking and divisoin information
    nlineage=sum(lastindex);
    [nROI,~]=size(prevROI);
    is_linked=zeros(1,nROI);
    tint=input('enter the time interval of timelapse: ');
    t=tint:tint:tint*nslice;
    %% separate the data per lineages
    %arrays for save
    allROI=zeros(nslice,nlineage);allch1=zeros(nslice,nlineage);alldivision=zeros(nslice,nlineage);allsize=zeros(nslice,nlineage);
    nlin=0;
    % separate into lineage
    for iROI=1:nROI
        lineagedata=zeros(nslice,6); %time/area/ch1/xm/ym/ROI
        if lastindex(iROI)==1
            nlin=nlin+1;
            % data separation per lieagedata
            ROI=iROI;
            for islice=1:nslice
                is_linked(1,ROI)=1;
                jslice=nslice+1-islice;
                lineagedata(jslice,:)=[3*(slice(ROI)),area(ROI),ch1mean(ROI),xm(ROI),ym(ROI),ROI];
                ROI=prevROI(ROI);
            end
            allROI(:,nlin)=lineagedata(:,6);
            allch1(:,nlin)=lineagedata(:,3);
            allsize(:,nlin)=lineagedata(:,2);
        end
        clear lineagedata;
    end
    
    Ntau=sum(lastindex);
    %% division count per lineage
    is_divide=divisionROI(prevROI,is_linked);
    for ilin=1:nlin
        alldivision(:,ilin)=division_boolean(allROI(:,ilin),is_divide);
    end
    
    %% write lineage tree
    nfig=0;nfig=nfig+1;
    figure(nfig);clf();
    write_lineagetree(allROI,alldivision,prevROI,nROI,nlineage,nslice,t);
    xlim([tint, tint*nslice]);
    %% calculate the initial number of cells
    ndiv=sum(alldivision);
    [~,nlineage]=size(ndiv);
    N0=0;
    for i=1:nlineage
        N0 = N0 + 2^(-ndiv(i));
    end
    consis=['Initial cell number and last cell number is ', num2str(N0), ', ', num2str(Ntau)];
    disp(consis);
    nplot=0;
    iplot=0;
    figure(10);
    plot(t,allch1);
    xlim([tint, tint*nslice]);
    hold on;
    figure(11);
    plot(t,allsize);
    xlim([tint, tint*nslice]);
    hold on;
    while(nplot<nlin)
        nplot=nplot+25;
        slin=nplot-24;ppp=min(nplot,nlin);
        i=0;nfig=nfig+1;
        figure(nfig);
        for ilin=slin:ppp
            iplot=iplot+1;
            i=i+1;
            subplot(5,5,i);
            yyaxis left;
            plot(t,allsize(:,ilin));
            yyaxis right;
            plot(t,allch1(:,ilin));
            title_lin=['lineage ',num2str(ilin)];
            title(title_lin);
            xlim([tint, tint*nslice]);
%             %color index
%             r = fix(iplot/100);
%             rplot =iplot -r*100;
%             g = fix(rplot/10);
%             b = rplot - g*10;
%             color = [0.1*r 0.1*g 0.1*b];
        end
    end
end