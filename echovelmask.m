clear
% fl1=ls('./nqcnc/*20150712*.cfradial');
% fl2=ls('./cilnc/*20150712*.cfradial');
% fl1=ls('./nqcnc/*20170512*.cfradial');
% fl2=ls('./cilnc/*20170512*.cfradial');
% fl1=ls('./nqcnc/*20150619*.cfradial');
% fl2=ls('./cilnc/*20150619*.cfradial');
% fl1=ls('./nqcnc/*20170725*.cfradial');
% fl2=ls('./cilnc/*20170725*.cfradial');
% fl1=ls('./nqcnc/*20141228*.cfradial');
% fl2=ls('./cilnc/*20141228*.cfradial');
% fl1=ls('./nqcnc/*20150215*.cfradial');
% fl2=ls('./cilnc/*20150215*.cfradial');
fl1=ls('./nqcnc/*.cfradial');
fl2=ls('./cilnc/*.cfradial');
fn=length(fl1);
j=0;
ref1mask=zeros(1000,720);
ref2mask=zeros(1000,720);
ldrthres=-10.5;
k=1;
maxtlen=720;
nanthres=10;
rainl=0;
rainli=1;
totalnum=zeros(12,720);
cloud_num=zeros(12,720);
BB_num=zeros(12,720);
rain_intsty_bin=zeros(12,720,16);
BB_height_bin=zeros(12,720,50);
% for i=1:1%fn
for i=1:fn
    fname1=strcat('./nqcnc/',fl1(i,:));
    matfname=strcat('mat/day',fname1(33:45));
    fname2=strcat('./cilnc/',fl2(i,:));
    ref2inst=ncread(fname2,'reflectivity_h');    
    load(matfname)
    clear rainflag* rain_intsty BB_height BB 
    rainmask=NaN(1000,720);
    rainflag=NaN(1,720);
    rainflag0=NaN(1,720);
    rainflag10=NaN(1,720);
    rain_intsty=NaN(1,720);
    rainmask(~isnan(ref2inst))=0;
    BB=zeros(1,720);
    BB_height=NaN(1,720);
    vel=ncread(fname2,'mean_doppler_velocity_h');
    k=0;
    nanmask=ncread(fname1,'nyquist_velocity');
    ldr=ncread(fname2,'linear_depolarization_ratio');
    z=ncread(fname2,'reflectivity_h');
    rainflag(nanmask>10)=0;
    fday=str2num(fl1(i,26:28));
    fyear=str2num(fl1(i,30:33));
    fmonth=str2num(matfname(17:18));
%% Check Total Number and Cloud amount
    for ti=2:maxtlen
        if nanmask(ti)>10   
            totalnum(fmonth,ti)=totalnum(fmonth,ti)+1;
            z_nan=length(z(~isnan(z(:,ti)),ti));
            if z_nan>0
                cloud_num(fmonth,ti)=cloud_num(fmonth,ti)+1;
            end
            clear z_nan
        end
%% Check Rain type/intensity        
        if ref2inst(21,ti)>-15&&nanmask(ti)>10
            rainflag(ti)=1;
            rain_intsty(ti)=ref2inst(21,ti);
            binno=floor((ref2inst(21,ti)+15)/2.5)+1;
            if binno>16
                binno=16;
            end
            rain_intsty_bin(fmonth,ti,binno)=rain_intsty_bin(fmonth,ti,binno)+1;
            if ref2inst(21,ti)>0
                rainflag0(ti)=1;
                if ref2inst(21,ti)>10
                    rainflag10(ti)=1;
                end
            end
            if length(rainl(rainl==refmask(21,ti)))<1
                rainl(rainli)=refmask(21,ti);
                rainls(rainli,1)=str2num(fname1(38:45));
                rainls(rainli,2)=ti;
                rainli=rainli+1;             
            end
            rainmask(1:21,ti)=1;
            for hi=21:1000
                if hi>20&&rainmask(hi-1,ti)==1
                    if (ref2inst(hi,ti)<-15||vel(hi,ti)>vel(20,ti)*0.6)
                        rainmask(hi,ti)=0;
                    else
                        rainmask(hi,ti)=1;
                    end
                elseif rainmask(hi-1,ti)==0
                    rainmask(hi,ti)=0;
                end
            end
        end
%% Check Melting Layer Presence and Height
        ldr_nan=length(ldr(~isnan(ldr(:,ti)),ti));
        max_z=max(z(31:1000,ti));
        max_ldr=max(ldr(31:1000,ti));
        ldrthres=-20;
        if max_z>0&&max_ldr<-20
            ldrthres=-30;
        end
        if ldr_nan>0&&nanmask(ti)>10&&max_z>-10
            BB_flag=0;
            z_diff=diff(z(:,ti));
            vel_diff=diff(vel(:,ti));
            ldr_diff=diff(ldr(:,ti));
            hi=30;
            k=0;
            while hi<1000&&~isempty(z(~isnan(z(hi:1000,ti)),ti))
%                 if ldr_diff(hi)<0&&ldr_diff(hi-1)>0&&ldr(hi,ti)>-40
                if (z_diff(hi)<0&&z_diff(hi-1)>0||z_diff(hi)==0)&&max(ldr(hi-8:hi+8,ti))>ldrthres&&z(hi,ti)>-10
                    max_diff=z_diff(hi-8);
                    min_diff=z_diff(hi+8);
                    hi_max=hi-8;
                    hi_min=hi+8;
                    for hhi=hi-8:hi+8
                        if z_diff(hhi)>max_diff
                            max_diff=z_diff(hhi);
                            hi_max=hhi;
                        end
                        if z_diff(hhi)<min_diff
                            min_diff=z_diff(hhi);
                            hi_min=hhi+1;
                        end
                    end
                    clear max_diff min_diff
                    slope_left=(z(hi,ti)-z(hi_max,ti))/(hi-hi_max);
                    slope_right=(z(hi,ti)-z(hi_min,ti))/(hi_min-hi);
                    if hi==hi_max
                        slope_left=0;
                    elseif hi==hi_min
                        slope_right=0;
                    end
                    ldr_nan_partial=length(ldr(~isnan(ldr(hi-6:hi+6,ti)),ti));
                    if ldr_nan_partial>0&&slope_left>=0.05&&slope_right>=0.1
                        max_ldr=max(ldr(hi-8:hi+8,ti));
                        if max_ldr~=ldr(hi-8,ti)&&max_ldr~=ldr(hi+8,ti)
                            for hhi=hi-8:hi+8

                                if ldr(hhi,ti)==max_ldr
                                    hi_xldr=hhi;
                                end
                            end
                        max_diff=max(ldr_diff(hi_xldr-8:hi_xldr));
                        min_diff=min(ldr_diff(hi_xldr:hi_xldr+8));
                            for hhi=hi_xldr-8:hi_xldr+8
                                if ldr_diff(hhi)==max_diff
                                    hi_max=hhi;
                                end
                                if ldr_diff(hhi)==min_diff
                                    hi_min=hhi+1; 
                                end
                            end
                           
                            slope_left=(ldr(hi_xldr,ti)-ldr(hi_max,ti))/(hi_xldr-hi_max);
                            slope_right=(ldr(hi_xldr,ti)-ldr(hi_min,ti))/(hi_min-hi_xldr);
%                             if max_ldr<-20
%                                 ldr_slope_thresh=0.3;
%                             elseif max_ldr<-10
                                ldr_slope_thresh=0.2;
%                             else
%                                 ldr_slope_thresh=0.1;
%                             end
%                             if slope_left>=ldr_slope_thresh&&slope_right>ldr_slope_thresh&&max(vel_diff(hi_xldr-2:hi_xldr+2))>0.06
                            if slope_left>=ldr_slope_thresh&&slope_right>=ldr_slope_thresh&&(nanmax(vel(hi_xldr:hi_xldr+2,ti))-nanmin(vel(hi_xldr-2:hi_xldr,ti)))>0.3
%                             if slope_left>=ldr_slope_thresh&&slope_right>=ldr_slope_thresh&&max(vel_diff(hi_xldr-2:hi_xldr+2))>0.1
                                k=k+1;
                                ldr_val(k)=max_ldr;
                                BB_loc(k)=hi_xldr;
                                ldr_slope(k)=slope_left+slope_right;
                                BB_flag=1;
                                BB_num(fmonth,ti)=BB_num(fmonth,ti)+1;
                                BB_h=hi_xldr*15;
                                BB_height(ti)=BB_h;
                                BB_bin=floor(BB_h/300)+1;
                                BB_height_bin(fmonth,ti,BB_bin)=BB_height_bin(fmonth,ti,BB_bin)+1;
                                BB(ti)=1;
                            end
                        end
                    end
%                     end
                end
                hi=hi+1;
            end 
            if k>1
                ldr_val_max=ldr_val(1);
                ldr_slope_max=ldr_slope(1);
                BB_h=BB_loc(1)*15;
                BB_height(ti)=BB_h;
                for ki=2:k
                    if ldr_val(ki)>ldr_val_max||(ldr_val(ki)==ldr_val_max&&ldr_slope(ki)>ldr_slope_max)
                        BB_h=BB_loc(k)*15; 
                        BB_height(ti)=BB_h;
                    end
                end
            end

            clear k ldr_val_max ldr_val ldr_slope BB_loc ldr_slope_max BB_h 
        end
    end
    rainmask(isnan(ref2inst))=NaN;
    
    save(matfname,'rainmask','rainflag','-append')
    save(matfname,'rainmask','rainflag0','-append')
    save(matfname,'rainmask','rainflag10','-append')
    save(matfname,'rainmask','rain_intsty','-append')
    save(matfname,'rainmask','BB_height','-append')
    save(matfname,'rainmask','BB','-append')
    clear rainflag* rain_intsty BB_height BB 

    clear BB BB_height rain_intsty rainflag rainflag0 rainflag10 z ldr vel
end
save('rainlist','rainl','rainls');
save('Cloud_BB_Rain','totalnum','cloud_num','BB_num','rain_intsty_bin','BB_height_bin');
% era5comp
BB_contour_figs