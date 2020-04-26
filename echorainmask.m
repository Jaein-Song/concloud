clear
fl1=ls('./nqcnc/*.cfradial');
fl2=ls('./cilnc/*.cfradial');
fn=length(fl1);
j=0;
maxtlen=720;
nanthres=10;
ref1mask=zeros(1000,720);   
ref2mask=zeros(1000,720);
ldrthres=-10.5;
k=1;
rainl=0;
rainli=1;
rains15=zeros(720,1);
rains00=zeros(720,1);
rains10=zeros(720,1);
rains15m=zeros(720,12);
rains00m=zeros(720,12);
rains10m=zeros(720,12);
ntotal=zeros(720,1);
ntotalm=zeros(720,12);
for i=1:fn
    fname1=strcat('./nqcnc/',fl1(i,:));
    matfname=strcat('mat/day',fname1(33:45));
    fname2=strcat('./cilnc/',fl2(i,:));
    ref2inst=ncread(fname2,'reflectivity_h');
    rainmask=NaN(1000,720);
    rainflag=NaN(1,720);
    rainflag00=NaN(1,720);
    rainflag10=NaN(1,720);
    rainmask(~isnan(ref2inst))=0;
%     vel=ncread(fname2,'mean_doppler_velocity_h');
    load(matfname)
    k=0;
    nanmask=ncread(fname1,'nyquist_velocity');
    rainflag(nanmask>10)=0;
    fday=str2num(fl1(i,26:28));
    fyear=str2num(fl1(i,30:33));
    nansatart=0;
    nanlength=zeros(720,1);
    m=str2num(matfname(17:18));
    for ti=1:maxtlen
        if nanmask(ti)>10
            ntotal(ti)=ntotal(ti)+1;
            ntotalm(ti,m)=ntotalm(ti,m)+1;
            if ref2inst(20,ti)>-15
                rainflag(ti)=1;
                rains15(ti)=rains15(ti)+1;
                rains15m(ti,m)=rains15m(ti,m)+1;
                if ref2inst(20,ti)>0
                    rainflag00(ti)=1;
                    rains00(ti)=rains00(ti)+1;
                    rains00m(ti,m)=rains00m(ti,m)+1;
                    if ref2inst(20,ti)>10
                        rainflag10(ti)=1;
                        rains10(ti)=rains10(ti)+1;
                        rains10m(ti,m)=rains10m(ti,m)+1;
                    end
                end
                rainmask(1:20,ti)=1;
            end
        end
    end
rainmask(isnan(ref2inst))=NaN;

    save(matfname,'rainmask','rainflag','-append');
    save(matfname,'rainmask','rainflag00','-append');
    save(matfname,'rainmask','rainflag10','-append');
end

save('rainechos','rains*','ntotal*')
% save('rainlist','rainl','rainls');
