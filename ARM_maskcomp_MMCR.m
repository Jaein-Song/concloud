site=siteo{1}
matdir=['./out/mat/' site '/'];
datadir=['~/CR_work/ARM/DATA/' site '/'];
flo=dir([matdir 'MMCR*.mat']);
j=1;
for i=1:length(flo)
    if flo(i).bytes>10000
        fl(j,:)=flo(i).name;
        j=j+1;
    end
end
clear flo
timezone_list=[-6 -9 -1 -9 0 12 12];
if site=="SGP"
    tz=timezone_list(1);
elseif site=="NSA"
    tz=timezone_list(2);
elseif site=="ENA"
    tz=timezone_list(3);
elseif site=="OLI"
    tz=timezone_list(4);
elseif site=="ASI"
    tz=timezone_list(5);
elseif site=="AWR"
    tz=timezone_list(6);
elseif site=="TWP"
    tz=timezone_list(7);
end
load([matdir fl(1,:)],'h');
load([matdir fl(1,:)],'t');
nt=length(t);
nh=length(h);
tz_os=tz*3600;
tz_de=mod(24-tz,24)*3600;
lst=mod(t+tz*3600+24*3600,24*3600);
tz_0h=find(lst==0);
tz_0l=length(t(tz_0h:nt));

fn=length(fl);
j=0;
nanthres=10;
ldrthres=-10.5;

mdiv=4;
if mdiv==4
elseif mdiv==1
    mlist=[1:1:12];
end
for ri=1:7
ref1mask=zeros(nh,nt);
ref2mask=zeros(nh,nt);
    for mi=1:mdiv
        if mdiv==4&&mi==1
            mlist=[12 1 2];
        elseif mdiv==4&&mi==2
            mlist=[3 4 5];
        elseif mdiv==4&&mi==3
            mlist=[6 7 8];
        elseif mdiv==4&&mi==4
            mlist=[9 10 11];
        end
        nums=zeros(nh,nt);
        ri
        mi
        for i=1:fn
            matfname=strcat([matdir fl(i,:)]);
            mfnl=length(matfname);
            ymd=matfname(mfnl-11:mfnl-4)
            y=str2num(ymd(1:4));
            m=str2num(ymd(5:6));
            d=str2num(ymd(7:8));
            if(max(m==mlist))
                load(matfname,'ref');
                ref2inst=ref;
                clear ref;
                
                if ri==1
                    load(matfname,'rainmask')
                    ref2mask(rainmask==0)=ref2mask(rainmask==0)+1;
                    clear rainmask
                elseif ri==2
                    ref2mask(~isnan(ref2inst))=ref2mask(~isnan(ref2inst))+1;
                elseif ri==3
                    load(matfname,'rainflag')
                    for ti=1:nt
                        if rainflag(ti)==1
                            ref2mask(~isnan(ref2inst(:,ti)),ti)=ref2mask(~isnan(ref2inst(:,ti)),ti)+1;
                        end
                    end
                    clear rainflag
                elseif ri==4
                    load(matfname,'rainflag')
                    for ti=1:nt
                        if rainflag(ti)~=1
                            ref2mask(~isnan(ref2inst(:,ti)),ti)=ref2mask(~isnan(ref2inst(:,ti)),ti)+1;
                        end
                    end
                    clear rainflag
                elseif ri==5
                    load(matfname,'rainmask')
                    ref2mask(rainmask==1)=ref2mask(rainmask==1)+1;
                    clear rainmask
                elseif ri==6
                    load(matfname,'velmask')
                    ref2mask(velmask==0)=ref2mask(velmask==0)+1;
                    clear refmask velmask
                elseif ri==7
                    load(matfname,'velmask')
                    ref2mask(velmask==1)=ref2mask(velmask==1)+1;
                    clear refmask velmask
                end
                load(matfname,'nanmask');
                nums(1:nh,nanmask>10)=nums(1:nh,nanmask>10)+1;
             clear *inst *1m nanmask
            end
        end
        ref2mask=ref2mask./nums;

        if mdiv==1
            refmaskTOT=ref2mask;
            refmaskTOTLST(:,1:tz_0l)=refmaskTOT(:,tz_0h:nt);
            refmaskTOTLST(:,tz_0l+1)=refmaskTOT(:,1:tz_0h-1);
        elseif mdiv==4&&mi==1
            refmaskDJF=ref2mask;
            refmaskDJFLST(:,1:tz_0l)=refmaskDJF(:,tz_0h:nt);
            refmaskDJFLST(:,tz_0l+1:nt)=refmaskDJF(:,1:tz_0h-1);
        elseif mdiv==4&&mi==2
            refmaskMAM=ref2mask;
            refmaskMAMLST(:,1:tz_0l)=refmaskMAM(:,tz_0h:nt);
            refmaskMAMLST(:,tz_0l+1:nt)=refmaskMAM(:,1:tz_0h-1);
        elseif mdiv==4&&mi==3
            refmaskJJA=ref2mask;
            refmaskJJALST(:,1:tz_0l)=refmaskJJA(:,tz_0h:nt);
            refmaskJJALST(:,tz_0l+1:nt)=refmaskJJA(:,1:tz_0h-1);
        elseif mdiv==4&&mi==4
            refmaskSON=ref2mask;
            refmaskSONLST(:,1:tz_0l)=refmaskSON(:,tz_0h:nt);
            refmaskSONLST(:,tz_0l+1:nt)=refmaskSON(:,1:tz_0h-1);
        end
        clear nums
    end
    t24=[-0.5:1:24.5];
    h24=[450:450:13500];
    h_num=450/(h(2)-h(1));
    t_num=3600/(t(2)-t(1));
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
%% OUTFILE SETUP FOR PLOTTING
% define var, name, xax xlim xlb, yax, ylim, ylb, prefix_filename
%define cmin, cmax, cint, cbint, cticks
            out{1,1}.var=refmaskMAMLST_24*100;
            out{1,2}.var=refmaskJJALST_24*100;
            out{1,3}.var=refmaskSONLST_24*100;
            out{1,4}.var=refmaskDJFLST_24*100;
            out{2,1}.var=refmaskMAMLST*100;
            out{2,2}.var=refmaskJJALST*100;
            out{2,3}.var=refmaskSONLST*100;
            out{2,4}.var=refmaskDJFLST*100;
            contourflag=0;
            contourline=2;

            if ri==1
                prefix_filename{1}='hourly_norain';
                prefix_filename{2}='2min_norain';
            elseif ri==2
                prefix_filename{1}='hourly_rain';
                prefix_filename{2}='2min_rain';
            elseif ri==3
                prefix_filename{1}='hourly_rain_col';
                prefix_filename{2}='2min_rain_col';
            elseif ri==4
                prefix_filename{1}='hourly_norain_col';
                prefix_filename{2}='2min_norain_col';
            elseif ri==5
                prefix_filename{1}='hourly_raincell';
                prefix_filename{2}='2min_raincell';
            elseif ri==6
                prefix_filename{1}='hourly_norainsys';
                prefix_filename{2}='2min_norainsys';
            elseif ri==7
                prefix_filename{1}='hourly_rainsys';
                prefix_filename{2}='2min_rainsys';
            end
            for oi=1:2
                out{oi,1}.name='MAM';
                out{oi,2}.name='JJA';
                out{oi,3}.name='SON';
                out{oi,4}.name='DJF';
                for oj=1:4
                    out{oi,oj}.cmin=0;
                    out{oi,oj}.cmax=20;
                    out{oi,oj}.cint=200;
                    out{oi,oj}.cbint=5;
                    out{oi,oj}.cb_name='%';
                        k=1;
                    for cti=out{oi,oj}.cmin:out{oi,oj}.cbint:out{oi,oj}.cmax
                        out{oi,oj}.cticks{k}=num2str(cti);
                        k=k+1;
                    end
                    out{1,oj}.xax=t24;
                    out{2,oj}.xax=t;
                    out{oi,oj}.xlim=[0 24];
                    out{oi,oj}.xlb='Hour (LST)';
                    out{1,oj}.yax=h24/1000;
                    out{2,oj}.yax=h/1000;
                    out{oi,oj}.ylim=[0 15];
                    out{oi,oj}.ylb='Height (km)';
                    out{oi,oj}.ctl=[10 10; 1 1];
                end
            end
%% CALL PLOTTING FUCTION 
try
Figure_make
catch
end

        if ri==1
            save(['MMCR' 'norain'  site],'out','refmask*');
        elseif ri==2
            save(['MMCR' 'rain'  site],'out','refmask*');
        elseif ri==3
            save(['MMCR' 'raincol'  site],'out','refmask*');
        elseif ri==4
            save(['MMCR' 'noraincol'  site],'out','refmask*');
        elseif ri==5
            save(['MMCR' 'raincell'  site],'out','refmask*');
        elseif ri==6
            save(['MMCR' 'norainsys'  site],'out','refmask*');
        elseif ri==7
            save(['MMCR' 'rainsys'  site],'out','refmask*');
        end
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

        if ri==1
            save(['MMCR' 'norain' site],'out','refmask*');
        elseif ri==2
            save(['MMCR' 'rain' site],'out','refmask*');
        elseif ri==3
            save(['MMCR' 'raincol' site],'out','refmask*');
        elseif ri==4
            save(['MMCR' 'noraincol' site],'out','refmask*');
        elseif ri==5
            save(['MMCR' 'raincell' site],'out','refmask*');
        elseif ri==6
            save(['MMCR' 'norainsys' site],'out','refmask*');
        elseif ri==7
            save(['MMCR' 'rainsys' site],'out','refmask*');
        end
    end
end

%     subplot(2,2,1);contourf(t,h,refmaskMAMLST,'linestyle','none');set(gca,'CLim',[0 0.25]);colorbar
%     subplot(2,2,2);contourf(t,h,refmaskJJALST,'linestyle','none');set(gca,'CLim',[0 0.25]);colorbar
%     subplot(2,2,3);contourf(t,h,refmaskSONLST,'linestyle','none');set(gca,'CLim',[0 0.25]);colorbar
%     subplot(2,2,4);contourf(t,h,refmaskDJFLST,'linestyle','none');set(gca,'CLim',[0 0.25]);colorbar    
% subplot(2,2,1);pcolor(t24,h24,refmaskMAMLST_24);shading flat;set(gca,'CLim',[0 0.2]);colorbar
%     set(gca,'CLim',[0 0.2],'FontSize',14,'FontWeight','bold','XTick',...
%     [0 2 4 6 8 10 12 14 16 18 20 22 24],'YTick',[0 3000 6000 9000 12000 15000],'xlim',[0 24],'ylim',[0 15000]);
%     title('MAM','FontWeight','bold','FontSize',16);
% subplot(2,2,2);pcolor(t24,h24,refmaskJJALST_24);shading flat;set(gca,'CLim',[0 0.25]);colorbar
%     set(gca,'CLim',[0 0.2],'FontSize',14,'FontWeight','bold','XTick',...
%     [0 2 4 6 8 10 12 14 16 18 20 22 24],'YTick',[0 3000 6000 9000 12000 15000],'xlim',[0 24],'ylim',[0 15000]);
%     title('JJA','FontWeight','bold','FontSize',16);
% subplot(2,2,3);pcolor(t24,h24,refmaskSONLST_24);shading flat;set(gca,'CLim',[0 0.25]);colorbar
%     set(gca,'CLim',[0 0.2],'FontSize',14,'FontWeight','bold','XTick',...
%     [0 2 4 6 8 10 12 14 16 18 20 22 24],'YTick',[0 3000 6000 9000 12000 15000],'xlim',[0 24],'ylim',[0 15000]);
%     title('SON','FontWeight','bold','FontSize',16);
% subplot(2,2,4);pcolor(t24,h24,refmaskDJFLST_24);shading flat;set(gca,'CLim',[0 0.25]);colorbar    
%     set(gca,'CLim',[0 0.2],'FontSize',14,'FontWeight','bold','XTick',...
%     [0 2 4 6 8 10 12 14 16 18 20 22 24],'YTick',[0 3000 6000 9000 12000 15000],'xlim',[0 24],'ylim',[0 15000]);
%     title('DJF','FontWeight','bold','FontSize',16);
% 
% subplot(2,2,1);pcolor(t,h,refmaskMAMLST);shading flat;set(gca,'CLim',[0 0.25]);colorbar
%     set(gca,'CLim',[0 0.15],'FontSize',14,'FontWeight','bold','XTick',...
%     [0 2 4 6 8 10 12 14 16 18 20 22 24],'YTick',[0 3000 6000 9000 12000 15000],'xlim',[0 24],'ylim',[0 15000]);
%     title('MAM','FontWeight','bold','FontSize',16);
%     subplot(2,2,2);pcolor(t,h,refmaskJJALST);shading flat;set(gca,'CLim',[0 0.25]);colorbar
%     set(gca,'CLim',[0 0.15],'FontSize',14,'FontWeight','bold','XTick',...
%     [0 2 4 6 8 10 12 14 16 18 20 22 24],'YTick',[0 3000 6000 9000 12000 15000],'xlim',[0 24],'ylim',[0 15000]);
%     title('JJA','FontWeight','bold','FontSize',16);
% subplot(2,2,3);pcolor(t,h,refmaskSONLST);shading flat;set(gca,'CLim',[0 0.25]);colorbar
%     set(gca,'CLim',[0 0.15],'FontSize',14,'FontWeight','bold','XTick',...
%     [0 2 4 6 8 10 12 14 16 18 20 22 24],'YTick',[0 3000 6000 9000 12000 15000],'xlim',[0 24],'ylim',[0 15000]);
%     title('SON','FontWeight','bold','FontSize',16);
% subplot(2,2,4);pcolor(t,h,refmaskDJFLST);shading flat;set(gca,'CLim',[0 0.25]);colorbar   
%     set(gca,'CLim',[0 0.15],'FontSize',14,'FontWeight','bold','XTick',...
%     [0 2 4 6 8 10 12 14 16 18 20 22 24],'YTick',[0 3000 6000 9000 12000 15000],'xlim',[0 24],'ylim',[0 15000]);
%     title('DJF','FontWeight','bold','FontSize',16);
%  
% else
%     for i=1:24
%         for hi=1:50
%             data=refmaskTOTLST;
%             instmat=data((hi-1)*20+1:hi*20,(i-1)*30+1:i*30);
%             refmaskTOTLST_24(hi,i)=mean(instmat(~isnan(instmat)));
%         end
%     end
%     refmaskTOTLST_24(:,1)=refmaskTOTLST_24(:,25);
%     refmaskTOTLST_24(:,26)=refmaskTOTLST_24(:,2);
%     pcolor(t24,h24,refmaskTOTLST_24);shading flat;set(gca,'Clim',[0 0.25]);colorbar
% end
