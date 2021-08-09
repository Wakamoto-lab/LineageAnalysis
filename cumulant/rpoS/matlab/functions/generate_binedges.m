function edges=generate_binedges(xc)
    [~,nbin]=size(xc);
    edges=zeros(1,nbin+1);
    wbin=xc(2)-xc(1);
    edges(1)=xc(1)-wbin/2;
    for i=1:nbin
        edges(i+1)=edges(i)+wbin;
    end
end