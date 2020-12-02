site=siteo{1}
fdir=['./out/mat/' site '/'];
BBdir = ['/home/jaein/CR_work/ARM/mat/' site '/'];
% fdir = './matfiles';
flist = dir([fdir '/day*.mat']);
listIndex = 1;
list2Index = 1;
hIntv = 15;
for fi =  1 : length(flist)
%     fname = flist(fi,:)
    try
        fname = [flist(fi).name]
        load([flist(fi).folder '/'  fname],'refmask');
        refSize = size(refmask);
        tLen = refSize(2);
        hLen = refSize(1);
        ymd = fname(5:12)
        try
            load([flist(fi).folder '/'  ymd '.mat'],'BB_height');
            BBflag = 1;
        catch
            BBflag = 0;
        end
        minRefMask = min(min(refmask));
    %     maxRefMask = max(max(refmask));     
        if ~isnan(minRefMask) 
            refMasks = refmask;
            refMasks(refMasks ~= minRefMask) = NaN;
            refMaskLength = length(find(refMasks==minRefMask));
            %refMaskLength = length(refMasks);
            while refMaskLength>0
                cloudTopIndex = 0;
                cloudBaseIndex = hLen;
                cloudTopPrev = 0;
                startIndex = 0;
                endIndex = 0;
                cloudTop = 0;
                cloudBase = 0;
                precip = 0;
                if BBflag
                    BBPresence = 0;
                else
                    BBPresence = NaN;
                end
                for ti = 1 : tLen
                    hi = hLen;
                    while isnan(refMasks(hi,ti)) && hi>1
                        hi = hi - 1;
                    end
                    if cloudTopIndex < hi && hi > 0
                        cloudTopIndex = hi;
                    end
                    if cloudTopPrev < 2  && cloudTopIndex > 1 && startIndex == 0
                        startIndex = ti;
                    elseif cloudTopPrev > 1 && hi == 1
                        endIndex = ti;
                    end
                    cloudTopPrev = hi;
                    hi = 1;
                    while isnan(refMasks(hi,ti)) && hi < hLen 
                        hi = hi + 1;
                    end
                    if hi <= 20 
                        precip = 1;
                    end
                    if BBflag && BBPresence == 0 && BB_height(ti)/hIntv >=hi && BB_height(ti)/hIntv <= cloudTopPrev
                        BBPresence =1;
                    end
                    if cloudBaseIndex > hi && hi > 20 && hi < hLen
                        cloudBaseIndex = hi;
                    end
                end
                if cloudBaseIndex == hLen
                    cloudBaseIndex = 20;
                end
                if endIndex == 0 
                    endIndex = 720;
                end
                maskYMD=(minRefMask - floor(minRefMask))*1e8;
                maskYY = floor(maskYMD/10000);
                maskMD = maskYMD - maskYY * 10000;
                maskMM = floor(maskMD/100);
                maskDD = floor(maskMD - maskMM*100);
                if endIndex/tLen - floor(endIndex/tLen) == 0 && endIndex > 0 
                    maxhits(list2Index,:) = [listIndex, minRefMask];
                    list2Index = list2Index + 1;
                end
                if list2Index > 1 && listIndex ~= maxhits(list2Index - 1,1) && ~isempty(find(maxhits(:,2)==minRefMask))
                    index = maxhits(find(maxhits(:,2)==minRefMask,1),1);
                  %  statistics(index,12) = max(statistics(index,12),BBPresence);
                    statistics(index,11) = max(statistics(index,11),precip);
                    statistics(index,10)=statistics(index,10)+endIndex;
                    statistics(index,8) = min(statistics(index,9),cloudBaseIndex);
                    statistics(index,7) = max(statistics(index,8),cloudTopIndex);
                else
                    statistics(listIndex,:) = [str2num(ymd),minRefMask,floor(minRefMask),maskYY, maskMM, maskDD , cloudTopIndex, cloudBaseIndex, startIndex, endIndex, precip];%, BBPresence];
                    listIndex = listIndex +1 ;
                end

                %%Update
                minRefMask=min(min(refmask(refmask>minRefMask)));
                if isempty(minRefMask)
                    break
                end
                refMasks = refmask;
                refMasks(refMasks ~= minRefMask) = NaN;
                %refMaskLength = length(refMasks);
                refMaskLength = length(find(refMasks==minRefMask));
            end
        end
    catch
    end
end
save([site 'statistics.mat'],'statistics')
