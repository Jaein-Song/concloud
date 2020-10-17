clear
% fl1=ls('./nqcnc/*.cfradial');
% fl2=ls('./cilnc/*.cfradial');
% fn=length(fl1);
j=0;
ldrthres=-10.5;

mdiv=4;
if mdiv==4
elseif mdiv==1
    mlist=[1:1:12];
end
    h=[15:15:15000];
    t24=[-0.5:1:24.5];
    h24=[300:300:15000];
    t=[120:120:86400]/3600;
for ri=1:7
    if ri==1
        load('norain');
    elseif ri==2
        load('rain');
    elseif ri==3
        load('raincol');
    elseif ri==4
        load('noraincol');
    elseif ri==5
        load('raincell');
    elseif ri==6
        load('norainsys');
    elseif ri==7
        load('rainsys');
    end
%% OUTFILE SETUP FOR PLOTTING
% define var, name, xax xlim xlb, yax, ylim, ylb, prefix_filename
%define cmin, cmax, cint, cbint, cticks

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
       
        for oi=1:2
            out{oi,1}.name='MAM';
            out{oi,2}.name='JJA';
            out{oi,3}.name='SON';
            out{oi,4}.name='DJF';
            for oj=1:4
                out{oi,oj}.cmin=0;
                out{oi,oj}.cmax=20;
                out{oi,oj}.cint=200;
                out{oi,oj}.cbint=2;
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
                ymax=15;
                yintv=2.5;
%                     out{oi,oj}.ctl=[15 15; 10 10;5 5; 1 1];
                out{oi,oj}.ctl=[ 10 10; 1 1];
            end
        end
%% CALL PLOTTING FUCTION 
        Figure_make

end


