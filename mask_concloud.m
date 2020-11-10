%% From left to right Propagation
for ti=1:maxtlen
    ref_nn=find(~isnan(ref(:,ti)));
    if ~isempty(ref_nn)&&length(ref_nn)<maxhlen
    %%Find Cloud Top and Bottom
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
            if cthi>cbhi||length(ref_nn)<maxhlen
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
                % veltestrange=min(cthi,m300);
                % if layeri==1&&ref(m300,ti)>-15&&cbhi<=veltestrange&&min(vel(cbhi:veltestrange,ti))<-1.5&&ti<maxtlen
                %     velmask(cbhi:cthi,ti)=refmask(cbhi:cthi,ti);
                %     velmask(find(~isnan(ref(cbhi:cthi,ti+1)))+cbhi-1,ti:ti+1)=min(min(refmask(cbhi:cthi,ti:ti+1),[],'omitnan'),[],'omitnan');
                % end
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
            if cthi>cbhi&&~isnan(min(refmask(cbhi:cthi,ti)))&&ti >1
                refmask(cbhi:cthi,ti)=min(refmask(cbhi:cthi,ti));
                mask_min=min(min(refmask(cbhi:cthi,ti-1:ti),[],'omitnan'),[],'omitnan');
                mask_max=max(max(refmask(cbhi:cthi,ti-1:ti),[],'omitnan'),[],'omitnan');
                while mask_min~=mask_max
                    refmask(refmask==mask_max)=mask_min;
                    % veltestrange=min(cthi,m300);
                    % if layeri==1&&ref(m300,ti)>-15&&cbhi<=veltestrange&&min(vel(cbhi:veltestrange,ti))<-1.5&&ti>1n
                    %     velmask(velmask==mask_max)=mask_min;
                    % end
                    mask_max=max(max(refmask(cbhi:cthi,ti-1:ti),[],'omitnan'),[],'omitnan');
                end
            end
        end
    end
    clear cbh cth cbhi cthi
    clear ref_nn
end
