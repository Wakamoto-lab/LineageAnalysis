function xc=generate_bincenter(xedges)
    [~,nedges]=size(xedges);
    nbinx=nedges-1;
    wbin=xedges(2)-xedges(1);
    xc=zeros(1,nbinx);
    for i=1:(nbinx)
        xc(i)=xedges(1)+(i-0.5)*wbin;
    end
end