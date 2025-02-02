site=siteo{1}
matdir=['./out/mat/' site '/'];
datadir=['~/CR_work/ARM/DATA/' site '/'];
BBdir=['~/CR_work/ARM/mat/' site '/'];
if isMMCR
    flo=dir([matdir 'MMCR*.mat']);
else
    flo=dir([matdir 'day*.mat']);
end
j=1;
clear fl
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
%nh=load([matdir fl(1,:)],'nh');
%nt=load([matdir fl(1,:)],'nt');
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
    i=i+1
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
    rainflag(~validmask)=0;
    fday=str2num(ymd(7:8));
    fyear=str2num(ymd(1:4));
    nansatart=0;
    nanlength=zeros(maxtlen,1);
    m=str2num(ymd(5:6));
    for ti=1:maxtlen
        if validmask(ti)
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
for ti=1:maxtlen
    if rainflag(ti)>0
        hi=1;
        while hi<nh&&isnan(refmask(hi,ti))
            hi=hi+1
        end
        if hi<nh ~isnan(refmask(hi,ti))
            rainmask(refmask==refmask(hi,ti))=refmask(hi,ti);
        end
    end
end
    matfname
    save(matfname,'rainmask','rainflag','-append');
    save(matfname,'rainmask','rainflag00','-append');
    save(matfname,'rainmask','rainflag10','-append');
    try
    BBfd=dir([BBdir,ymd,'.mat']);
    catch
    end

   if ~isempty(BBfd)
       load([BBfd.folder '/' BBfd.name])
        save(matfname,'rainmask','BB*','-append');
   end 
end
if isMMCR
save(['MMCRrainechos' site],'rains*','ntotal*')
else
save(['KAZRrainechos' site],'rains*','ntotal*')
end
ARM_maskcomp
% save('rainlist','rainl','rainls');
