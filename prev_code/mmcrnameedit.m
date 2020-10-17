%%Consecutive Cloud Detector CONCLUDE
site=siteo{1}
%site='SGP'
matdir=['./out/mat/' site '/'];
datadir=['~/ARM_CRML/MMCR/' site '/'];
flo=dir([datadir '*.cdf']);
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
k=1;
for i=1:fn
    %% get vars
    fname=strcat([datadir ],fl(i,:))
    fnl=length(fname);
    ymd=fname(fnl-18:fnl-11);
%    try
    movefile(strcat(matdir,'day_',ymd,'.mat'),[matdir,'MMCR_day_',ymd,'.mat'])
clear ref* *mask nanlength nanstart
%catch
%    disp(['error: ' fl(i,:)])
%end    
end
