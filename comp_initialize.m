site=siteo{1}
matdir=['./out/mat/' site '/'];
datadir=['~/CR_work/ARM/DATA/' site '/'];
if isMMCR
    flo=dir([matdir 'MMCR*.mat']);
    radartype='MMCR';
else
    flo=dir([matdir 'day*.mat']);
    radartype = 'KAZR';
end
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

r{1}.prefix_filename{1}='hourly_total';
r{1}.prefix_filename{2}='2min_total';
r{2}.prefix_filename{1}='hourly_rainsys';
r{2}.prefix_filename{2}='2min_rainsys';
r{3}.prefix_filename{1}='hourly_clear';
r{3}.prefix_filename{2}='2min_clear';
r{4}.prefix_filename{1}='hourly_raincol';
r{4}.prefix_filename{2}='2min_raincol';
r{5}.prefix_filename{1}='hourly_clearcol';
r{5}.prefix_filename{2}='2min_clearcol';

r{1}.outFileName = [radartype 'total' site];
r{2}.outFileName = [radartype 'rainsys' site];
r{3}.outFileName = [radartype 'clear' site];
r{4}.outFileName = [radartype 'raincol' site];
r{5}.outFileName = [radartype 'clearcol' site];


if mdiv==1
    months{1}.mlist=[1:1:12];
else
    months{1}.mlist=[12 1 2];
    months{2}.mlist=[3 4 5];
    months{3}.mlist=[6 7 8];
    months{4}.mlist=[9 10 11];
end
