for ti = 1 : maxtlen
    if isnan(velmask(m300,ti))&&ref(m300,ti)>-15&&vel(m300,ti)<-1.5
        refmaskRain = refmask(m300,ti);
        velmask(refmask == refmaskRain) = refmaskRain;
    end
end