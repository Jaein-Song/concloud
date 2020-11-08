refmask=NaN(maxhlen,maxtlen);
velmask=NaN(maxhlen,maxtlen);
ref=ncread(fname,radar.refname{radarindex});
vel=ncread(fname,radar.Vname{radarindex});
if isMMCR
    ref=double(ref)/100;
    vel=double(vel)/1000*-1;
end

ref(ref<-100)=NaN;
vel(vel<-100)=NaN;
try
    ldr=ncread(fname,radar.LDRname{radarindex});
    ldr(vel<-100)=NaN;
catch
    ldr=NaN(size(ref));
end
k=0;
iaf=ncread(fname,radar.IAFname{radarindex}); 
if radarindex == 1
    mid = iaf;
    iaf=mid(1,:)*1;
    clear mid
    validmask=zeros(1,length(iaf));
    validmask=mod(iaf,10)==0;
else
    iaf(iaf<-100)=0;
    validmask=zeros(1,length(iaf));
    validmask=mod(iaf,2);
end
clear iaf 
nansatart=0;
nanlength=zeros(maxtlen,1);
%%INITIALIZE END    

%% Check if the start of file is not valid
if ~validmask(1)
    nanstart=1;
    nanlength(1)=1;
end
%% Count non-valid columns
for ti=2:maxtlen
    if ~validmask(ti-1)&&~validmask(ti)
        nanlength(nanstart:ti)=nanlength(nanstart)+1;
    elseif validmask(ti-1)&&~validmask(ti)
        nanstart=ti;
        nanlength(ti)=1;
    end
end
%% Fill nonvalid cells 
for ti=1:maxtlen
    if ti==1&&~validmask(ti)&&nanlength(ti)<nanthres
        ref(:,ti)=ref(:,nanlength(ti)+1)
    elseif ti>1&&~validmask(ti)&&nanlength(ti)<nanthres
        ref(:,ti)=ref(:,ti-1);
    end
end