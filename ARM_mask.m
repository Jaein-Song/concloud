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
    %try
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
            mask_initialize
            mask_concloud
            mask_timegap
            mask_paday
            mask_rain
            %% Masking renumber
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

            %%Nan values process
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
            
            %% Save
            disp(strcat(num2str(i),'/',num2str(fn),': ','outfile:',matdir,'/day_',ymd))
            if isMMCR
                save(strcat(matdir,'/MMCR_day_',ymd),'t','h','validmask','ref','refmask','ldr','vel','velmask','nanlength')
            else
                save(strcat(matdir,'/day_',ymd),'t','h','validmask','ref','refmask','ldr','vel','velmask','nanlength')
            end

            clear ref* *mask nanlength nanstart
        end
    %catch
    %    disp(['error' fl(i,:)])
    %end    
end
ARM_maskcomp
%clear r* f* m*  t* h* n* 
% ARM_rainmask
%echovelmask;
%h=[15:15:15000];
%t=[120:120:86400]/3600;
