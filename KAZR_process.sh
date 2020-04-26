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
   # cmdline="siteo=$site;try;ARM_mask;catch;end;exit"
    cmdline="siteo=$site;try;rainmask_KAZR;catch;end;exit"
   echo $cmdline 
    nohup matlab -r $cmdline -nodesktop -nosplash -nojvm -logfile $logfile > ${sites[$i]}_KAZR.log & 
    let i++
done
