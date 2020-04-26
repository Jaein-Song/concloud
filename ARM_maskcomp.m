clear
site=siteo{1}
matdir=['./out/mat/' site '/'];
datadir=['~/CR_work/ARM/DATA/' site '/'];
flo=dir([matdir '*.mat']);
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
ldrthres=-10.5;
h=load([matdir fl(1,:)],'h');
t=load([matdir fl(1,:)],'t');
nh=load([matdir fl(1,:)],'nh');
nt=load([matdir fl(1,:)],'nt');

mdiv=1;
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
            ymd=matfname(mfnl-11:mfnl-4);
            if(max(m==mlist))
                ref2inst2inst=load(matfname,'ref');
                %ref2inst=ncread(fname2,'reflectivity_h');
                %vel2inst=fname(matfname,'vel');
                %vel2inst=ncread(fname2,'mean_doppler_velocity_h');
                nanmask=ncread(fname2,'nyquist_velocity');
                
                if ri==1
                    load(matfname,'rainmask')
%                     ref1mask(~isnan(ref1inst)&(vel1inst>-1.5))=ref1mask(~isnan(ref1inst)&(vel1inst>-1.5))+1;
                    ref2mask(rainmask==0)=ref2mask(rainmask==0)+1;
                    clear rainmask
                elseif ri==2
%                     ref1mask(~isnan(ref1inst))=ref1mask(~isnan(ref1inst))+1;
                    ref2mask(~isnan(ref2inst))=ref2mask(~isnan(ref2inst))+1;
                elseif ri==3
                    matfname=strcat('mat/day',fname1(33:45));
                    load(matfname,'rainflag')

                    for ti=1:720
                        if rainflag(ti)==1
                            ref2mask(~isnan(ref2inst(:,ti)),ti)=ref2mask(~isnan(ref2inst(:,ti)),ti)+1;
                        end
                    end
                    clear rainflag
                elseif ri==4
                    matfname=strcat('mat/day',fname1(33:45));
                    load(matfname,'rainflag')
                    for ti=1:720
                        if rainflag(ti)~=1
                            ref2mask(~isnan(ref2inst(:,ti)),ti)=ref2mask(~isnan(ref2inst(:,ti)),ti)+1;
                        end
                    end
                    clear rainflag
                elseif ri==5
                    matfname=strcat('mat/day',fname1(33:45));
                    load(matfname,'rainmask')
                    ref2mask(rainmask==1)=ref2mask(rainmask==1)+1;
                    clear rainmask
                elseif ri==6
                    matfname=strcat('mat/day',fname1(33:45));
                    load(matfname,'velmask')
%                     for ti=1:720
%                         if vel2inst(20,ti)<-1
%                             ref2inst(refmask==refmask(20,ti))=NaN;
%                         end
%                     end
                    ref2mask(velmask==0)=ref2mask(velmask==0)+1;
                    clear refmask velmask
                elseif ri==7
                    matfname=strcat('mat/day',fname1(33:45));
                    load(matfname,'velmask')
%                     ref2inst2=NaN(1000,720);
%                     for ti=1:720
%                         if vel2inst(20,ti)<-1
%                             ref2inst2(refmask==refmask(20,ti))=ref2inst(refmask==refmask(20,ti));
%                         end
%                     end
                    ref2mask(velmask==1)=ref2mask(velmask==1)+1;
                    clear refmask velmask
                end
                nums(1:1000,nanmask>10)=nums(1:1000,nanmask>10)+1;
             clear *inst *1m nanmask
            end
        end
%         ref1mask=ref1mask./nums;
        ref2mask=ref2mask./nums;
        if mdiv==1
            refmaskTOT=ref2mask;
            refmaskTOTKST(:,1:270)=refmaskTOT(:,451:720);
            refmaskTOTKST(:,271:720)=refmaskTOT(:,1:450);
        elseif mdiv==4&&mi==1
            refmaskDJF=ref2mask;
            refmaskDJFKST(:,1:270)=refmaskDJF(:,451:720);
            refmaskDJFKST(:,271:720)=refmaskDJF(:,1:450);
        elseif mdiv==4&&mi==2
            refmaskMAM=ref2mask;
            refmaskMAMKST(:,1:270)=refmaskMAM(:,451:720);
            refmaskMAMKST(:,271:720)=refmaskMAM(:,1:450);
        elseif mdiv==4&&mi==3
            refmaskJJA=ref2mask;
            refmaskJJAKST(:,1:270)=refmaskJJA(:,451:720);
            refmaskJJAKST(:,271:720)=refmaskJJA(:,1:450);
        elseif mdiv==4&&mi==4
            refmaskSON=ref2mask;
            refmaskSONKST(:,1:270)=refmaskSON(:,451:720);
            refmaskSONKST(:,271:720)=refmaskSON(:,1:450);
        end
        clear nums
    end
    h=[15:15:15000];
    t24=[-0.5:1:24.5];
    h24=[300:300:15000];
    t=[120:120:86400]/3600;
    if mdiv==4
        for i=1:24
            for hi=1:50
                data=refmaskMAMKST;
                instmat=data((hi-1)*20+1:hi*20,(i-1)*30+1:i*30);
                refmaskMAMKST_24(hi,i+1)=mean(instmat(~isnan(instmat)));
                clear data instmat
                data=refmaskJJAKST;
                instmat=data((hi-1)*20+1:hi*20,(i-1)*30+1:i*30);
                refmaskJJAKST_24(hi,i+1)=mean(instmat(~isnan(instmat)));
                clear data instmat
                data=refmaskSONKST;
                instmat=data((hi-1)*20+1:hi*20,(i-1)*30+1:i*30);
                refmaskSONKST_24(hi,i+1)=mean(instmat(~isnan(instmat)));
                clear data instmat
                data=refmaskDJFKST;
                instmat=data((hi-1)*20+1:hi*20,(i-1)*30+1:i*30);
                refmaskDJFKST_24(hi,i+1)=mean(instmat(~isnan(instmat)));
                clear data instmat
            end
        end
        refmaskMAMKST_24(:,1)=refmaskMAMKST_24(:,25);
        refmaskMAMKST_24(:,26)=refmaskMAMKST_24(:,2);
        refmaskJJAKST_24(:,1)=refmaskJJAKST_24(:,25);
        refmaskJJAKST_24(:,26)=refmaskJJAKST_24(:,2);
        refmaskSONKST_24(:,1)=refmaskSONKST_24(:,25);
        refmaskSONKST_24(:,26)=refmaskSONKST_24(:,2);
        refmaskDJFKST_24(:,1)=refmaskDJFKST_24(:,25);
        refmaskDJFKST_24(:,26)=refmaskDJFKST_24(:,2);
%% OUTFILE SETUP FOR PLOTTING
% define var, name, xax xlim xlb, yax, ylim, ylb, prefix_filename
%define cmin, cmax, cint, cbint, cticks
            out{1,1}.var=refmaskMAMKST_24*100;
            out{1,2}.var=refmaskJJAKST_24*100;
            out{1,3}.var=refmaskSONKST_24*100;
            out{1,4}.var=refmaskDJFKST_24*100;
            out{2,1}.var=refmaskMAMKST*100;
            out{2,2}.var=refmaskJJAKST*100;
            out{2,3}.var=refmaskSONKST*100;
            out{2,4}.var=refmaskDJFKST*100;
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
                    out{oi,oj}.xlb='Hour (KST)';
                    out{1,oj}.yax=h24/1000;
                    out{2,oj}.yax=h/1000;
                    out{oi,oj}.ylim=[0 15];
                    out{oi,oj}.ylb='Height (km)';
                    out{oi,oj}.ctl=[10 10; 1 1];
                end
            end
%% CALL PLOTTING FUCTION 
%             Figure_make

        if ri==1
            save('norain','*ref*mask*','-append');
        elseif ri==2
            save('rain','*ref*mask*','-append');
        elseif ri==3
            save('raincol','*ref*mask*','-append');
        elseif ri==4
            save('noraincol','*ref*mask*','-append');
        elseif ri==5
            save('raincell','*ref*mask*','-append');
        elseif ri==6
            save('norainsys','*ref*mask*','-append');
        elseif ri==7
            save('rainsys','*ref*mask*','-append');
        end
    else
      
        
        for i=1:24
            for hi=1:50
                data=refmaskTOTKST;
                instmat=data((hi-1)*20+1:hi*20,(i-1)*30+1:i*30);
                refmaskTOTKST_24(hi,i+1)=mean(instmat(~isnan(instmat)));
                clear data instmat
            end
        end
        refmaskTOTKST_24(:,1)=refmaskTOTKST_24(:,25);
        refmaskTOTKST_24(:,26)=refmaskTOTKST_24(:,2);

        if ri==1
            save('norain','*ref*mask*','-append');
        elseif ri==2
            save('rain','*ref*mask*','-append');
        elseif ri==3
            save('raincol','*ref*mask*','-append');
        elseif ri==4
            save('noraincol','*ref*mask*','-append');
        elseif ri==5
            save('raincell','*ref*mask*','-append');
        elseif ri==6
            save('norainsys','*ref*mask*','-append');
        elseif ri==7
            save('rainsys','*ref*mask*','-append');
        end
    end
end

%     subplot(2,2,1);contourf(t,h,refmaskMAMKST,'linestyle','none');set(gca,'CLim',[0 0.25]);colorbar
%     subplot(2,2,2);contourf(t,h,refmaskJJAKST,'linestyle','none');set(gca,'CLim',[0 0.25]);colorbar
%     subplot(2,2,3);contourf(t,h,refmaskSONKST,'linestyle','none');set(gca,'CLim',[0 0.25]);colorbar
%     subplot(2,2,4);contourf(t,h,refmaskDJFKST,'linestyle','none');set(gca,'CLim',[0 0.25]);colorbar    
% subplot(2,2,1);pcolor(t24,h24,refmaskMAMKST_24);shading flat;set(gca,'CLim',[0 0.2]);colorbar
%     set(gca,'CLim',[0 0.2],'FontSize',14,'FontWeight','bold','XTick',...
%     [0 2 4 6 8 10 12 14 16 18 20 22 24],'YTick',[0 3000 6000 9000 12000 15000],'xlim',[0 24],'ylim',[0 15000]);
%     title('MAM','FontWeight','bold','FontSize',16);
% subplot(2,2,2);pcolor(t24,h24,refmaskJJAKST_24);shading flat;set(gca,'CLim',[0 0.25]);colorbar
%     set(gca,'CLim',[0 0.2],'FontSize',14,'FontWeight','bold','XTick',...
%     [0 2 4 6 8 10 12 14 16 18 20 22 24],'YTick',[0 3000 6000 9000 12000 15000],'xlim',[0 24],'ylim',[0 15000]);
%     title('JJA','FontWeight','bold','FontSize',16);
% subplot(2,2,3);pcolor(t24,h24,refmaskSONKST_24);shading flat;set(gca,'CLim',[0 0.25]);colorbar
%     set(gca,'CLim',[0 0.2],'FontSize',14,'FontWeight','bold','XTick',...
%     [0 2 4 6 8 10 12 14 16 18 20 22 24],'YTick',[0 3000 6000 9000 12000 15000],'xlim',[0 24],'ylim',[0 15000]);
%     title('SON','FontWeight','bold','FontSize',16);
% subplot(2,2,4);pcolor(t24,h24,refmaskDJFKST_24);shading flat;set(gca,'CLim',[0 0.25]);colorbar    
%     set(gca,'CLim',[0 0.2],'FontSize',14,'FontWeight','bold','XTick',...
%     [0 2 4 6 8 10 12 14 16 18 20 22 24],'YTick',[0 3000 6000 9000 12000 15000],'xlim',[0 24],'ylim',[0 15000]);
%     title('DJF','FontWeight','bold','FontSize',16);
% 
% subplot(2,2,1);pcolor(t,h,refmaskMAMKST);shading flat;set(gca,'CLim',[0 0.25]);colorbar
%     set(gca,'CLim',[0 0.15],'FontSize',14,'FontWeight','bold','XTick',...
%     [0 2 4 6 8 10 12 14 16 18 20 22 24],'YTick',[0 3000 6000 9000 12000 15000],'xlim',[0 24],'ylim',[0 15000]);
%     title('MAM','FontWeight','bold','FontSize',16);
%     subplot(2,2,2);pcolor(t,h,refmaskJJAKST);shading flat;set(gca,'CLim',[0 0.25]);colorbar
%     set(gca,'CLim',[0 0.15],'FontSize',14,'FontWeight','bold','XTick',...
%     [0 2 4 6 8 10 12 14 16 18 20 22 24],'YTick',[0 3000 6000 9000 12000 15000],'xlim',[0 24],'ylim',[0 15000]);
%     title('JJA','FontWeight','bold','FontSize',16);
% subplot(2,2,3);pcolor(t,h,refmaskSONKST);shading flat;set(gca,'CLim',[0 0.25]);colorbar
%     set(gca,'CLim',[0 0.15],'FontSize',14,'FontWeight','bold','XTick',...
%     [0 2 4 6 8 10 12 14 16 18 20 22 24],'YTick',[0 3000 6000 9000 12000 15000],'xlim',[0 24],'ylim',[0 15000]);
%     title('SON','FontWeight','bold','FontSize',16);
% subplot(2,2,4);pcolor(t,h,refmaskDJFKST);shading flat;set(gca,'CLim',[0 0.25]);colorbar   
%     set(gca,'CLim',[0 0.15],'FontSize',14,'FontWeight','bold','XTick',...
%     [0 2 4 6 8 10 12 14 16 18 20 22 24],'YTick',[0 3000 6000 9000 12000 15000],'xlim',[0 24],'ylim',[0 15000]);
%     title('DJF','FontWeight','bold','FontSize',16);
%  
% else
%     for i=1:24
%         for hi=1:50
%             data=refmaskTOTKST;
%             instmat=data((hi-1)*20+1:hi*20,(i-1)*30+1:i*30);
%             refmaskTOTKST_24(hi,i)=mean(instmat(~isnan(instmat)));
%         end
%     end
%     refmaskTOTKST_24(:,1)=refmaskTOTKST_24(:,25);
%     refmaskTOTKST_24(:,26)=refmaskTOTKST_24(:,2);
%     pcolor(t24,h24,refmaskTOTKST_24);shading flat;set(gca,'Clim',[0 0.25]);colorbar
% end
