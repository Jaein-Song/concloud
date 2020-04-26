#!/bin/sh
sites=(`cat sitelist_MMCR`)
i=0
sitenum=${#sites}
while [ $i -lt $sitenum ]; do
    echo ${sites[$i]}
    site=\"${sites[$i]}\"
    echo $site
    logfile=MMCR_${sites[$i]}.log
    echo $logfile 
    errmsg="\"Error occurred: site = ${sites[$i]}\""
    echo $errmsg
    cmdline="siteo=$site;try;mmcrnameedit;catch;end;exit"
   echo $cmdline 
    nohup matlab -r $cmdline -nodesktop -nosplash -nojvm -logfile $logfile > ${sites[$i]}_MMCR.log & 
    let i++
done
