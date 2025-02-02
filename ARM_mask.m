%%Consecutive Cloud Detector CONCLUDE
site=siteo{1}
matdir=['./out/mat/' site '/'];
if isMMCR
    radarindex = 1;
else
    radarindex = 2;
end
radar.tlen=[8640,21600];
radar.hlen=[512,596];
radar.datadir={['~/ARM_CRML/MMCR/' site '/'],['~/CR_work/ARM/DATA/' site '/']};
radar.refname={['ReflectivityBestEstimate'],['reflectivity_best_estimate']};
radar.Vname={['MeanDopplerVelocity'],['mean_doppler_velocity']};
radar.LDRname={[],['linear_depolarization_ratio']};
radar.IAFname={['ModeId'],['instrument_availability_flag']};
radar.tname={['time_offset'],['time']};
radar.hname={['Heights'],['height']};
radar.minByte=[10000000,20000000];
radar.expansion={['*.cdf'],['*.nc']}
flo=dir([radar.datadir{radarindex} radar.expansion{radarindex}])
j=1;
floSize=size(flo);
for i=1:floSize(1)
    if flo(i).bytes > radar.minByte(radarindex)
        fl(j,:)=flo(i).name;
        j=j+1;
    end
end
clear flo
fn=length(fl);
j=0;
nanthres=10;
k=1;
for i=1:fn
    %%INITIALIZE
    %% Get vars
    fname=strcat(radar.datadir{radarindex} ,fl(i,:))
    try
        t=ncread(fname,radar.tname{radarindex});
        h=ncread(fname,radar.hname{radarindex});
        hi = 1
        while h(hi)<300
            hi=hi+1;
        end
        m300=hi;
        td=t(2)-t(1);
        nanthres=1200/td;
        maxtlen=length(t);
        maxhlen=length(h);
        if maxtlen==radar.tlen(radarindex)&&maxhlen==radar.hlen(radarindex)
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
            %% From left to right Propagation
            for ti=1:maxtlen
                ref_nn=find(~isnan(ref(:,ti)));
                if ~isempty(ref_nn)&&length(ref_nn)<maxhlen
                    cbhi=ref_nn(1);
                    cthi=0;
                    cth=ref_nn(diff(ref_nn)>1);
                    if isempty(cth)
                        cbh=ref_nn(1);
                        cth=ref_nn(length(ref_nn));
                        layer=1;
                    else 
                        cbh(1)=ref_nn(1);
                        layeri=1;
                        for layerj=1:length(ref_nn)
                            if cth(layeri)==ref_nn(layerj)
                                layeri=layeri+1; 
                                cbh(layeri)=ref_nn(layerj+1);
                            end
                            if layeri >length(cth)
                                break
                            end
                        end
                        layer=length(cbh);
                        cth(layer)=ref_nn(length(ref_nn));
                    end

                    for layeri=1:layer
                        cbhi=cbh(layeri);
                        cthi=cth(layeri);
                        if cthi>=cbhi||length(ref_nn)<maxhlen
                            if isnan(min(refmask(cbhi:cthi,ti)))
                                refmask(cbhi:cthi,ti)=k;
                                if ti <maxtlen
                                    refmask(find(~isnan(ref(cbhi:cthi,ti+1)))+cbhi-1,ti:ti+1)=min(min(refmask(cbhi:cthi,ti:ti+1),[],'omitnan'),[],'omitnan');
                                end
                                k=k+1;
                            elseif min(refmask(cbhi:cthi,ti))<k
                                refmask(cbhi:cthi,ti)=min(refmask(cbhi:cthi,ti));
                                if ti <maxtlen
                                    refmask(find(~isnan(ref(cbhi:cthi,ti+1)))+cbhi-1,ti:ti+1)=min(min(refmask(cbhi:cthi,ti:ti+1),[],'omitnan'),[],'omitnan');
                                end
                            end
                            veltestrange=min(cthi,m300*3);
                            if layeri==1&&ref(m300,ti)>-15&&cbhi>veltestrange&&min(vel(cbhi:veltestrange,ti))<-1.5&&ti<maxtlen
                                velmask(cbhi:cthi,ti)=refmask(cbhi:cthi,ti);
                                velmask(find(~isnan(ref(cbhi:cthi,ti+1)))+cbhi-1,ti:ti+1)=min(min(refmask(cbhi:cthi,ti:ti+1),[],'omitnan'),[],'omitnan');
                            end
                        end
                    end
                end
                clear cbh cth cbhi cthi layer 
                clear ref_nn
            end

            %% From Right to Left Propagation
            for ti=maxtlen:-1:1
                ref_nn=find(~isnan(ref(:,ti)));
                if ~isempty(ref_nn)
                    cbhi=ref_nn(1);
                    cthi=0;
                    cth=ref_nn(diff(ref_nn)>1);
                    if isempty(cth)
                        cbh=ref_nn(1);
                        cth=ref_nn(length(ref_nn));
                        layer=1;
                    else 
                        cbh(1)=ref_nn(1);
                        layeri=1;
                        for layerj=1:length(ref_nn)
                            if cth(layeri)==ref_nn(layerj)
                                layeri=layeri+1; 
                                cbh(layeri)=ref_nn(layerj+1);
                            end
                            if layeri >length(cth)
                                break
                            end
                        end
                        layer=length(cbh);
                        cth(layer)=ref_nn(length(ref_nn));
                    end

                    for layeri=1:layer
                        cbhi=cbh(layeri);
                        cthi=cth(layeri);
                        if cthi>=cbhi
                            if ~isnan(min(refmask(cbhi:cthi,ti)))
                                refmask(cbhi:cthi,ti)=min(refmask(cbhi:cthi,ti));
                                if ti >1
                                    mask_min=min(min(refmask(cbhi:cthi,ti-1:ti),[],'omitnan'),[],'omitnan');
                                    mask_max=max(max(refmask(cbhi:cthi,ti-1:ti),[],'omitnan'),[],'omitnan');
                                    while mask_min~=mask_max
                                        refmask(refmask==mask_max)=mask_min;
                                        veltestrange=min(cthi,m300*3);
                                        if layeri==1&&ref(m300,ti)>-15&&cbhi>veltestrange&&min(vel(cbhi:veltestrange,ti))<-1.5&&ti>1n
                                            velmask(velmask==mask_max)=mask_min;
                                        end
                                        mask_max=max(max(refmask(cbhi:cthi,ti-1:ti),[],'omitnan'),[],'omitnan');
                                    end
                                end
                            end
                        end
                    end
                end
                clear cbh cth cbhi cthi
                clear ref_nn
            end
            fnl=length(fname);
            if isMMCR
                fday=str2num(fname(fnl-12:fnl-11));
                fmonth=str2num(fname(fnl-14:fnl-13));
                fyear=str2num(fname(fnl-18:fnl-15));
                fmd=str2num(fname(fnl-14:fnl-11));
                ymd=fname(fnl-18:fnl-11);
            else
                fday=str2num(fname(fnl-11:fnl-10));
                fmonth=str2num(fname(fnl-13:fnl-12));
                fyear=str2num(fname(fnl-17:fnl-14));
                fmd=str2num(fname(fnl-13:fnl-10));
                ymd=fname(fnl-17:fnl-10);
            end
            [pyr pmn pda]=paday(1,fyear,fmonth,fday);
            if ~isempty(find(~isnan(refmask(:,1))))
                findl=dir([matdir,'*',num2str(pyr),num2str(pmn,'%02d'),num2str(pda,'%02d'),'*mat']);
                if ~isempty(findl)
                    pfname=findl.name;
                    pfdir=findl.folder;
                    pfilen=[pfdir '/' pfname]
                    clear findl
                    pfilen=strcat(pfilen);
                    prev=load(pfilen);
                    minmask=min(refmask(:,1));
                    maxmask=max(refmask(:,1));

                    for hi=1:maxhlen
                        if ~isnan(prev.refmask(hi,maxtlen))&&~isnan(refmask(hi,1))
                            
                            if mod(refmask(hi,1),1)==0
                                
                                refmask(refmask==refmask(hi,1))=prev.refmask(hi,maxtlen);
                            else
                                min_mask=min(prev.refmask(hi,maxtlen),refmask(hi,1));
                                max_mask=max(prev.refmask(hi,maxtlen),refmask(hi,1));
                                if min_mask~=max_mask
                                    prev.refmask(prev.refmask==max_mask)=min_mask;
                                    refmask(refmask==max_mask)=min_mask;
                                    if isfield(prev,'velmask')&&~isnan(prev.velmask(hi,maxtlen))&&isnan(velmask(hi,1))
                                        velmask(refmask==refmask(hi,1))=prev.velmask(hi,maxtlen);
                                    elseif isfield(prev,'velmask')&&isnan(prev.velmask(hi,maxtlen))&&~isnan(velmask(hi,1))
                                        prev.velmask(refmask==refmask(hi,maxtlen))=velmask(hi,1);
                                    end
                                end
                            end
                        end
                    end
                    refmask_pres=refmask;
                    velmask_pres=velmask;
                    velmask=prev.refmask;
                    refmask=prev.refmask;
                    save(pfilen,'refmask','velmask','-append');
                    clear refmask velmask
                    velmask=velmask_pres;
                    refmask=refmask_pres;
                    clear refmask_pres;
                    clear velmask_pres
                    clear prev 
                    clear pfilen
                else
                    for hi=1:maxhlen
                        mask=refmask(hi,1);
                        refmask(refmask==mask)=NaN;
                        velmask(velmask==mask)=NaN;
                    end
                end
                
            end
            if ~isempty(find(~isnan(refmask(:,maxtlen))))
            [ayr amn ada]=paday(0,fyear,fmonth,fday);
                findl=dir([radar.datadir{radarindex},'*',num2str(ayr),num2str(amn,'%02d'),num2str(ada,'%02d'),'*nc']);
                afilen=[radar.datadir{radarindex} findl.name]
                clear findl
                if isempty(afilen)
                    for hi=1:maxhlen
                        mask=refmask(hi,maxtlen);
                        refmask(refmask==mask)=NaN;
                        velmask(velmask==mask)=NaN;
                    end
                end
                clear afilen
            end
            k=1+str2num(ymd)*10^(-8);

            for mi=floor(min(min(refmask))):floor(max(max(refmask)))
                if ~isempty(refmask==mi)
                    if length(find(refmask==mi))>4
                        refmask(refmask==mi&refmask-mi<0.1)=k;
                        k=k+1;
                    else
                        mi;
                        refmask(refmask==mi)=NaN;
                    end
                end
            end
            if min(nanlength)>nanthres
                for ni=1:maxtlen-1
                    if nanlength(ni)>nanthres&&nanlength(ni+1)<nanthres
                        for hi=1:maxhlen
                            refmask(refmask==refmask(hi,ni+1))=NaN;
                        end
                    elseif nanlength(ni)<nanthres&&nanlength(ni+1)>nanthres
                        for hi=1:maxhlen
                            refmask(refmask==refmask(hi,ni))=NaN;
                        end
                    end
                    if ni==1&&nanlength(ni)>nanthres
                    end
                end
            end
            disp(strcat(num2str(i),'/',num2str(fn),': ','outfile:',matdir,'/day_',ymd))
            if isMMCR
                save(strcat(matdir,'/MMCR_day_',ymd),'t','h','validmask','ref','refmask','ldr','vel','velmask','nanlength')
            else
                save(strcat(matdir,'/day_',ymd),'t','h','validmask','ref','refmask','ldr','vel','velmask','nanlength')
            end

        clear ref* *mask nanlength nanstart
        end
    catch
        disp(['error' fl(i,:)])
    end    
end
%ARM_maskcomp
%clear r* f* m*  t* h* n* 
ARM_rainmask
%echovelmask;
%h=[15:15:15000];
%t=[120:120:86400]/3600;
