site=siteo{1}
matdir=['./out/mat/' site '/'];
datadir=['~/CR_work/ARM/DATA/' site '/'];
flo=dir([matdir 'MMCR*.mat']);
j=1;
for i=1:length(flo)
    if flo(i).bytes>10000
        fl(j,:)=flo(i).name;
        j=j+1;
    end
end
clear flo
fn=length(fl);
j=0;
nanthres=10;
ldrthres=-10.5;
load([matdir fl(1,:)],'h');
load([matdir fl(1,:)],'t');
nh=length(h);
nt=length(t);
maxtlen=nt;
nanthres=10;
ref1mask=zeros(nh,nt);   
ref2mask=zeros(nh,nt);
ldrthres=-10.5;
k=1;
rainl=0;
rainli=1;
rains15=zeros(maxtlen,1);
rains00=zeros(maxtlen,1);
rains10=zeros(maxtlen,1);
rains15m=zeros(maxtlen,12);
rains00m=zeros(maxtlen,12);
rains10m=zeros(maxtlen,12);
ntotal=zeros(maxtlen,1);
ntotalm=zeros(maxtlen,12);
i=1;
while h(i)<300
    i=i+1;
end
m300=i;
for i=1:fn
    matfname=strcat([matdir,fl(i,:)]);
    load(matfname)
    rainmask=NaN(nh,maxtlen);
    rainflag=NaN(1,maxtlen);
    rainflag00=NaN(1,maxtlen);
    rainflag10=NaN(1,maxtlen);
    rainmask(~isnan(ref))=0;
    mfnl=length(matfname);
    ymd=matfname(mfnl-11:mfnl-4);
    k=0;
    rainflag(nanmask>10)=0;
    fday=str2num(ymd(7:8));
    fyear=str2num(ymd(1:4));
    nansatart=0;
    nanlength=zeros(maxtlen,1);
    m=str2num(ymd(5:6));
    for ti=1:maxtlen
        if nanmask(ti)>10
            ntotal(ti)=ntotal(ti)+1;
            ntotalm(ti,m)=ntotalm(ti,m)+1;
            if ref(m300,ti)>-15
                rainflag(ti)=1;
                rains15(ti)=rains15(ti)+1;
                rains15m(ti,m)=rains15m(ti,m)+1;
                if ref(m300,ti)>0
                    rainflag00(ti)=1;
                    rains00(ti)=rains00(ti)+1;
                    rains00m(ti,m)=rains00m(ti,m)+1;
                    if ref(m300,ti)>10
                        rainflag10(ti)=1;
                        rains10(ti)=rains10(ti)+1;
                        rains10m(ti,m)=rains10m(ti,m)+1;
                    end
                end
                rainmask(1:m300,ti)=1;
            end
        end
    end
rainmask(isnan(ref))=NaN;
    matfname
    save(matfname,'rainmask','rainflag','-append');
    save(matfname,'rainmask','rainflag00','-append');
    save(matfname,'rainmask','rainflag10','-append');
end

save(['MMCRrainechos_' site],'rains*','ntotal*')
ARM_maskcomp_MMCR
% save('rainlist','rainl','rainls');
