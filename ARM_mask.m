%%Consecutive Cloud Detector CONCLUDE
site=siteo{1}
%site='SGP'
matdir=['./out/mat/' site '/'];
datadir=['~/CR_work/ARM/DATA/' site '/'];
flo=dir([datadir '*.nc']);
j=1;
for i=1:length(flo)
    if flo(i).bytes>1000000
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
    %% get vars
    fname=strcat([datadir ],fl(i,:))
    try
    t=ncread(fname,'time');
    h=ncread(fname,'height');
    maxtlen=length(t);
    maxhlen=length(h);
    if maxtlen==21600&&maxhlen==596
    refmask=NaN(maxhlen,maxtlen);
    ref_inst=ncread(fname,'reflectivity_best_estimate');
    ref_inst(ref_inst<-100)=NaN;
    ref=ref_inst;
    vel=ncread(fname,'mean_doppler_velocity');
    vel(vel<-100)=NaN;
    ldr=ncread(fname,'linear_depolarization_ratio');
    ldr(vel<-100)=NaN;
    k=0;
    iaf=ncread(fname,'instrument_availability_flag'); 
    iaf(iaf<-100)=NaN;
    nanmask=zeros(1,length(iaf));
    nanmask(mod(iaf,2)==1)=11;
    nansatart=0;
    nanlength=zeros(maxtlen,1);
    
    %% Check if the start of file is not valid
    if nanmask(1)<10
       nanstart=1;
       nanlength(1)=1;
    end
    %% Count non-valid columns
    for ti=2:maxtlen
        if nanmask(ti-1)<10&&nanmask(ti)<10
            nanlength(nanstart:ti)=nanlength(nanstart)+1;
        elseif nanmask(ti)<10
            nanstart=ti;
            nanlength(ti)=1;
        end
    end
    %% Fill nonvalid cells 
    for ti=1:maxtlen
        if ti==1&&nanmask(ti)<10&&nanlength(ti)<nanthres
            ref_inst(:,ti)=ref_inst(:,nanlength(ti)+1);
        elseif ti>1&&nanmask(ti)<10&&nanlength(ti)<nanthres
            ref_inst(:,ti)=ref_inst(:,ti-1);
        end
    end
    %% From left to right Propagation
    for ti=1:maxtlen
        ref_nn=find(~isnan(ref_inst(:,ti)));
        if ~isempty(ref_nn)
            cbhi=ref_nn(1);
            cthi=0;
            for nnhi=1:length(ref_nn)
                if nnhi==length(ref_nn)
                    cthi=ref_nn(nnhi);
                elseif ref_nn(nnhi)<ref_nn(nnhi+1)-1
                    cthi=ref_nn(nnhi);
                end
                
                if cthi>=cbhi
                    if isnan(min(refmask(cbhi:cthi,ti)))
                        refmask(cbhi:cthi,ti)=k;
                        if ti <maxtlen
                            refmask(find(~isnan(ref_inst(cbhi:cthi,ti+1)))+cbhi-1,ti:ti+1)=min(min(refmask(cbhi:cthi,ti:ti+1),[],'omitnan'),[],'omitnan');
                        end
                        k=k+1;
                    elseif min(refmask(cbhi:cthi,ti))<k
                        refmask(cbhi:cthi,ti)=min(refmask(cbhi:cthi,ti));
                        if ti <maxtlen
                            refmask(find(~isnan(ref_inst(cbhi:cthi,ti+1)))+cbhi-1,ti:ti+1)=min(min(refmask(cbhi:cthi,ti:ti+1),[],'omitnan'),[],'omitnan');
                        end
                    end
                    if nnhi<length(ref_nn)
                        cbhi=ref_nn(nnhi+1);

                    end
                end
                
            end
        end
        clear ref_nn
    end
    %% From Right to Left Propagation
    for ti=maxtlen:-1:1
        ref_nn=find(~isnan(ref_inst(:,ti)));
        if ~isempty(ref_nn)
            cbhi=ref_nn(1);
            cthi=0;
           
            for nnhi=1:length(ref_nn)
                if nnhi==length(ref_nn)
                    cthi=ref_nn(nnhi);
                elseif ref_nn(nnhi)<ref_nn(nnhi+1)-1
                    cthi=ref_nn(nnhi);
                end
                
                if cthi>=cbhi
                    if ~isnan(min(refmask(cbhi:cthi,ti)))
                        
                        refmask(cbhi:cthi,ti)=min(refmask(cbhi:cthi,ti));
                        if ti >1
                            mask_min=min(min(refmask(cbhi:cthi,ti-1:ti),[],'omitnan'),[],'omitnan');
                            mask_max=max(max(refmask(cbhi:cthi,ti-1:ti),[],'omitnan'),[],'omitnan');
                            while mask_min~=mask_max
                                refmask(refmask==mask_max)=mask_min;
                                mask_max=max(max(refmask(cbhi:cthi,ti-1:ti),[],'omitnan'),[],'omitnan');
                            end
                        end
                    end
                    if  (nnhi)<length(ref_nn)
                        cbhi=ref_nn(nnhi+1);
                    end
                end
            end
        end
        clear ref_nn
    end
    fnl=length(fname);
    fday=str2num(fname(fnl-11:fnl-10));
    fmonth=str2num(fname(fnl-13:fnl-12));
    fyear=str2num(fname(fnl-17:fnl-14));
    fmd=str2num(fname(fnl-13:fnl-10));
    [pyr pmn pda]=paday(1,fyear,fmonth,fday);
    if ~isempty(find(~isnan(refmask(:,1))))
        %if fmd==101
        %    findl=strcat(matdir,'/*',num2str(fyear-1),num2str(1231),'.mat');
        %else
        %    findl=strcat(matdir,'/day_',num2str(fday-1,'%03d'),'_',num2str(fyear),'*.mat');
        %end
        %pfilen=ls(findl);
        findl=dir(['*',matdir,num2str(pyr),num2str(pmn,'%02d'),num2str(pda,'%02d'),'*mat']);
        pfilen=[matdir findl.name]

        if ~isempty(findl)
            clear findl
            pfilen=strcat(matdir,'/',pfilen);
            prev=load(pfilen);
            minmask=min(refmask(:,1));
            maxmask=max(refmask(:,1));
            ie(:,i*2-1)=prev.refmask(:,maxtlen);
            ie(:,i*2)=refmask(:,1);

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
                        end
                    end
                end
            end
            ie2(:,i*2-1)=prev.refmask(:,maxtlen);
            ie2(:,i*2)=refmask(:,1);
            refmask_pres=refmask;
            refmask=prev.refmask;
            save(pfilen,'refmask','-append');
            clear refmask;
            refmask=refmask_pres;
            clear refmask_pres
            clear prev 
            clear pfilen
        else
            for hi=1:maxhlen
                mask=refmask(hi,1);
                refmask(refmask==mask)=NaN;
            end
        end
        
    end

    if ~isempty(find(~isnan(refmask(:,maxtlen))))
    [ayr amn ada]=paday(0,fyear,fmonth,fday);
        findl=dir([datadir,num2str(ayr),num2str(amn,'%02d'),num2str(ada,'%02d'),'*nc']);
        afilen=[datadir findl.name]
        clear findl
        %if fmd==1231
        %    findl=strcat('cilnc/*',num2str(fyear+1),num2str(101,'%04d'),'*.cfradial');
        %else
        %    findl=strcat('cilnc/*',num2str(fday+1,'%03d'),'_',num2str(fyear),'*.cfradial');
        %end
        %afilen=ls(findl);
        %clear findl
        if isempty(afilen)
            for hi=1:maxhlen
                mask=refmask(hi,maxtlen);
                refmask(refmask==mask)=NaN;
            end
        end
        clear afilen
    end
    k=1+str2num(fname(fnl-17:fnl-10))*10^(-8);

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
    disp(strcat('outfile:',matdir,'/day_',fname(fnl-17:fnl-10)))
    save(strcat(matdir,'/day_',fname(fnl-17:fnl-10)),'t','h','nanmask','ref','refmask','ldr','vel','nanlength')
clear ref* *mask nanlength nanstart
end
catch
    disp(['error' fl(i,:)])
end    
end
%echovelmask;
%h=[15:15:15000];
%t=[120:120:86400]/3600;
