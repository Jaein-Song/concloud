function [ ry, rm, rd ] = paday (pa_flag, yr, mn, da)
%% pa_flag==1: past else after
    if pa_flag
        if da~=1
            ry=yr;rm=mn;rd=da-1;
        elseif mn==2||mn==4||mn==6||mn==8||mn==9||mn==11
            ry=yr;rm=mn-1;rd=31;
        elseif mn==3
            ry=yr;rm=2;
            if mod(yr,4)==0
                rd=29;
            else
                rd=28;
            end
        elseif mn==1
            ry=yr-1;rm=12;rd=31;
        else
            ry=yr;rm=mn-1;rd=30;
        end
    else
        eoy=0;
       if mn==2&&mod(yr,4)==0
           eday=29;
       elseif mn==2&&mod(yr,4)~=0
           eday=28;
       elseif mn==4||mn==6||mn==9||mn==11
           eday=30;
       else
           eday=31;
       end
       if da~=eday
           ry=yr;rm=mn;rd=da+1;
       elseif da==eday&&mn<12
           ry=yr;rm=mn+1;rd=1;
       else
           ry=yr+1;rm=1;rd=1;
       end
    end
end

