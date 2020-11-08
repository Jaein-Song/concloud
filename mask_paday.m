if ~isempty(find(~isnan(refmask(:,1))))
    [pyr pmn pda]=paday(1,fyear,fmonth,fday);
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