#!/bin/sh
sites=(`cat sitelist_KAZR`)
i=0
sitenum=${#sites[*]}
while [ $i -lt $sitenum ]; do
    echo ${sites[$i]}
    site=\"${sites[$i]}\"
    echo $site
    logfile=KAZR_${sites[$i]}.log
    echo $logfile 
    errmsg="\"Error occurred: site = ${sites[$i]}\""
    cmdline="siteo=$site;isMMCR=0;ARM_mask;exit"
    #cmdline="siteo=$site;isMMCR=0;try;ARM_rainmask;catch;end;exit"
    #cmdline="siteo=$site;try;rainmask_KAZR;catch;end;exit"
   echo $cmdline 
    nohup matlab -r $cmdline -nodesktop -nosplash -nojvm -logfile $logfile  & 
    let i++
done
sites=(`cat sitelist_MMCR`)
i=0
sitenum=${#sites[*]}
while [ $i -lt $sitenum ]; do
    echo ${sites[$i]}
    site=\"${sites[$i]}\"
    echo $site
    logfile=MMCR_${sites[$i]}.log
    echo $logfile 
    errmsg="\"Error occurred: site = ${sites[$i]}\""
    cmdline="siteo=$site;isMMCR=1;try;ARM_mask;catch;end;exit"
    #cmdline="siteo=$site;isMMCR=1;try;ARM_rainmask;catch;end;exit"
    #cmdline="siteo=$site;try;rainmask_KAZR;catch;end;exit"
   echo $cmdline 
    nohup matlab -r $cmdline -nodesktop -nosplash -nojvm -logfile $logfile & 
    let i++
done
