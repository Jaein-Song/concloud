gapInMin = 4;
gap = gapInMin * 60 / 86400 * maxtlen;

% minRefMask = min(min(refmask))
% maxRefMask = max(max(refmask))
% for rmi = minRefMask : maxRefMask
%     refmaskLength = length(find(refmask(refmask == rmi)));
%     if refmaskLength > 0

%     end
% end

for hi = 1 : maxhlen
    for ti = 1 : maxtlen-1
        tRightEnd = min(maxtlen, ti + gap);
        if ~isnan(refmask(hi,ti))&&isnan(refmask(hi,ti+1)) && ~isnan(min(refmask(hi,ti+1:tRightEnd)))
            for tj = ti : tRightEnd
                minRefMask = min(refmask(hi,ti+1:tRightEnd));
                if refmask(hi,tj) ~= minRefMask
                    refmask(refmask==refmask(hi,tj))=minRefMask;
                end
            end
        end
    end
end