clear
%%Consecutive Cloud Detector CONCLUDE
matdir=['.out/mat' site]
fl=ls('./cilnc/*.cfradial');
fn=length(fl);
j=0;
maxtlen=length(t);
maxhlen=length(h);
nanthres=10;
ref1mask=zeros(maxhlen,720);
ref2mask=zeros(maxhlen,720);
ldrthres=-10.5;
k=1;
for i=1:fn
    %% get vars
    refmask=NaN(maxhlen,720);
    fname1=strcat('./nqcnc/',fl(i,:));
    fname2=strcat('./cilnc/',fl(i,:));
    ref2inst=ncread(fname2,'reflectivity_h');
    ref=ref2inst;
    vel=ncread(fname2,'mean_doppler_velocity_h');
    ldr=ncread(fname2,'linear_depolarization_ratio');
    k=0;
    nanmask=ncread(fname1,'nyquist_velocity');
    fday=str2num(fl(i,26:28));
    fyear=str2num(fl(i,30:33));
    nansatart=0;
    nanlength=zeros(720,1);
    
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
            ref2inst(:,ti)=ref2inst(:,nanlength(ti)+1);
        elseif ti>1&&nanmask(ti)<10&&nanlength(ti)<nanthres
            ref2inst(:,ti)=ref2inst(:,ti-1);
        end
    end
    %% From left to right Propagation
    for ti=1:maxtlen
        ref_nn=find(~isnan(ref2inst(:,ti)));
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
                            refmask(find(~isnan(ref2inst(cbhi:cthi,ti+1)))+cbhi-1,ti:ti+1)=min(min(refmask(cbhi:cthi,ti:ti+1),[],'omitnan'),[],'omitnan');
                        end
                        k=k+1;
                    elseif min(refmask(cbhi:cthi,ti))<k
                        refmask(cbhi:cthi,ti)=min(refmask(cbhi:cthi,ti));
                        if ti <maxtlen
                            refmask(find(~isnan(ref2inst(cbhi:cthi,ti+1)))+cbhi-1,ti:ti+1)=min(min(refmask(cbhi:cthi,ti:ti+1),[],'omitnan'),[],'omitnan');
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
        ref_nn=find(~isnan(ref2inst(:,ti)));
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
    fday=str2num(fl(i,26:28));
    fyear=str2num(fl(i,30:33));
    fmd=str2num(fl(i,34:37));

    if ~isempty(find(~isnan(refmask(:,1))))
        if fmd==101
            findl=strcat(matdir,'/*',num2str(fyear-1),num2str(1231),'.mat');
        else
            findl=strcat(matdir,'/day_',num2str(fday-1,'%03d'),'_',num2str(fyear),'*.mat');
        end
        pfilen=ls(findl);
        clear findl

        if ~isempty(pfilen)
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
        if fmd==1231
            findl=strcat('cilnc/*',num2str(fyear+1),num2str(101,'%04d'),'*.cfradial');
        else
            findl=strcat('cilnc/*',num2str(fday+1,'%03d'),'_',num2str(fyear),'*.cfradial');
        end
        afilen=ls(findl);
        clear findl
        if isempty(afilen)
            for hi=1:maxhlen
                mask=refmask(hi,maxtlen);
                refmask(refmask==mask)=NaN;
            end
        end
        clear afilen
    end
    k=1+str2num(fl(i,30:37))*10^(-8);

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
    save(strcat(matdir,'/day_',fl(i,26:37)),'refmask','vel','ldr','nanlength')
clear ref* *mask nanlength nanstart
end
%echovelmask;
%h=[15:15:15000];
%t=[120:120:86400]/3600;
