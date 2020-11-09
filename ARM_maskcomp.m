

mdiv=4;
comp_initialize
nums=zeros(nh,nt);
for i=1:fn
    matfname=strcat([matdir fl(i,:)]);
    mfnl=length(matfname);
    ymd=matfname(mfnl-11:mfnl-4)
    y=str2num(ymd(1:4));
    m=str2num(ymd(5:6));
    d=str2num(ymd(7:8));
    for mfi = 1 : length(months)
        if max(m==months{mfi}.mlist)
            mi = mfi
        end
    end 
    r{ri}.ref2mask=zeros(nh,nt);
    load(matfname,'velmask');
    load(matfname,'refmask');
    load(matfname,'rainmask');
    load(matfname,'validmask');
    for ri=1:5
        if  ri==1
            loc = ~isnan(refmask);
            clear velmask refmask rainmask 
        elseif ri==2
            loc=~isnan(velmask);
            clear refmask velmask
        elseif ri==3
            loc = ~isnan(refmask)&isnan(velmask);
            clear refmask velmask rainmask
        elseif ri == 4
            loc = zeros(nh,nt);
            for ti = 1:length(rainloc)
                if rainmask(1,ti)==1
                    loc(:,ti) = ~isnan(refmask(:,ti));
                end
            end
            clear refmask rainmask velmask
        elseif ri == 5
            loc = zeros(nh,nt);
            for ti = 1:length(rainloc);
                if rainmask(1,ti)==0
                    loc(:,ti) = ~isnan(refmask(:,ti));
                end
            end
            clear refmask rainmask velmask
        end
        r{ri}.ref2mask(loc)=r{ri}.ref2mask(loc)+1;
        nums(1:nh,validmask>0)=nums(1:nh,validmask>0)+1;
        clear *inst *1m validmask 
    end
end

for ri = 1:5
    r{ri}.ref2mask=r{ri}.ref2mask./nums;

    if mdiv==1
        refmaskTOT=r{ri}.ref2mask;
        refmaskTOTLST(:,1:tz_0l)=refmaskTOT(:,tz_0h:nt);
        refmaskTOTLST(:,tz_0l+1)=refmaskTOT(:,1:tz_0h-1);
    elseif mdiv==4&&mi==1
        refmaskDJF=r{ri}.ref2mask;
        refmaskDJFLST(:,1:tz_0l)=refmaskDJF(:,tz_0h:nt);
        refmaskDJFLST(:,tz_0l+1:nt)=refmaskDJF(:,1:tz_0h-1);
    elseif mdiv==4&&mi==2
        refmaskMAM=r{ri}.ref2mask;
        refmaskMAMLST(:,1:tz_0l)=refmaskMAM(:,tz_0h:nt);
        refmaskMAMLST(:,tz_0l+1:nt)=refmaskMAM(:,1:tz_0h-1);
    elseif mdiv==4&&mi==3
        refmaskJJA=r{ri}.ref2mask;
        refmaskJJALST(:,1:tz_0l)=refmaskJJA(:,tz_0h:nt);
        refmaskJJALST(:,tz_0l+1:nt)=refmaskJJA(:,1:tz_0h-1);
    elseif mdiv==4&&mi==4
        refmaskSON=r{ri}.ref2mask;
        refmaskSONLST(:,1:tz_0l)=refmaskSON(:,tz_0h:nt);
        refmaskSONLST(:,tz_0l+1:nt)=refmaskSON(:,1:tz_0h-1);
    end
    clear nums
    t24 = [-0.5:1:24.5];
    h24 = [450:450:13500];
    h_num = 450/(h(2)-h(1));
    t_num = 3600/(t(2)-t(1));
    if mdiv==4
        for i=1:24
            for hi=1:length(h24)
                data=refmaskMAMLST;
                instmat=data((hi-1)*h_num+1:hi*h_num,(i-1)*t_num+1:i*t_num);
                refmaskMAMLST_24(hi,i+1)=mean(instmat(~isnan(instmat)));
                clear data instmat
                data=refmaskJJALST;
                instmat=data((hi-1)*h_num+1:hi*h_num,(i-1)*t_num+1:i*t_num);
                refmaskJJALST_24(hi,i+1)=mean(instmat(~isnan(instmat)));
                clear data instmat
                data=refmaskSONLST;
                instmat=data((hi-1)*h_num+1:hi*h_num,(i-1)*t_num+1:i*t_num);
                refmaskSONLST_24(hi,i+1)=mean(instmat(~isnan(instmat)));
                clear data instmat
                data=refmaskDJFLST;
                instmat=data((hi-1)*h_num+1:hi*h_num,(i-1)*t_num+1:i*t_num);
                refmaskDJFLST_24(hi,i+1)=mean(instmat(~isnan(instmat)));
                clear data instmat
            end
        end
        refmaskMAMLST_24(:,1)=refmaskMAMLST_24(:,25);
        refmaskMAMLST_24(:,26)=refmaskMAMLST_24(:,2);
        refmaskJJALST_24(:,1)=refmaskJJALST_24(:,25);
        refmaskJJALST_24(:,26)=refmaskJJALST_24(:,2);
        refmaskSONLST_24(:,1)=refmaskSONLST_24(:,25);
        refmaskSONLST_24(:,26)=refmaskSONLST_24(:,2);
        refmaskDJFLST_24(:,1)=refmaskDJFLST_24(:,25);
        refmaskDJFLST_24(:,26)=refmaskDJFLST_24(:,2);
    else
        for i=1:24
            for hi=1:length(h24)
                data=refmaskTOTLST;
                instmat=data((hi-1)*h_num+1:hi*h_num,(i-1)*t_num+1:i*t_num);
                refmaskTOTLST_24(hi,i+1)=mean(instmat(~isnan(instmat)));
                clear data instmat
            end
        end
        refmaskTOTLST_24(:,1)=refmaskTOTLST_24(:,25);
        refmaskTOTLST_24(:,26)=refmaskTOTLST_24(:,2);
    end
    comp_outfile
end