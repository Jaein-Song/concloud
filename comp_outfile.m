for mi = 1:4
r{1,mi}.outFileName = [radartype 'total' site];
r{2,mi}.outFileName = [radartype 'rainsys' site];
r{3,mi}.outFileName = [radartype 'clear' site];
end
%r{4}.outFileName = [radartype 'raincol' site];
%r{5}.outFileName = [radartype 'clearcol' site];
if mdiv==4
    %% OUTFILE SETUP FOR PLOTTING
    % define var, name, xax xlim xlb, yax, ylim, ylb, prefix_filename
    %define cmin, cmax, cint, cbint, cticks
    out{1,1}.var=refmaskMAMLST_24*100;
    out{1,2}.var=refmaskJJALST_24*100;
    out{1,3}.var=refmaskSONLST_24*100;
    out{1,4}.var=refmaskDJFLST_24*100;
    out{2,1}.var=refmaskMAMLST{ri}*100;
    out{2,2}.var=refmaskJJALST{ri}*100;
    out{2,3}.var=refmaskSONLST{ri}*100;
    out{2,4}.var=refmaskDJFLST{ri}*100;
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
    save(r{ri}.outFileName,'out','refmask*','contourflag','contourline');
else
    save(r{ri}.outFileName,'out','refmask*','contourflag','contourline');
end
