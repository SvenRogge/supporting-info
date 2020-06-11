#!/bin/sh

prefix='m'
suffix='MPa'
sign=''

for i in m1000MPa m500MPa m400MPa m300MPa m200MPa m100MPa 0MPa 100MPa 200MPa 300MPa 400MPa 500MPa 600MPa 700MPa 800MPa 900MPa 1000MPa 1100MPa 1200MPa 1300MPa 1400MPa 1500MPa 1600MPa 1700MPa 1800MPa 1900MPa 2000MPa

do
    if [[ $i == ${prefix}* ]]; then 
        sign='-'
        press=${i#$prefix}
    else
        sign=''
        press=${i}
    fi

    press=${press%$suffix}

    mkdir ../${i}
    cp ./ymd.py ../${i}/
    sed -i "s/press=1e6\*pascal/press=${sign}${press}e6*pascal/" ../${i}/ymd.py
    cp ./job.sh ../${i}/
done

for i in m1000MPa m500MPa m400MPa m300MPa m200MPa m100MPa 0MPa 100MPa 200MPa 300MPa 400MPa 500MPa 600MPa 700MPa 800MPa 900MPa 1000MPa 1100MPa 1200MPa 1300MPa 1400MPa 1500MPa 1600MPa 1700MPa 1800MPa 1900MPa 2000MPa

do
    cd ../${i}
    qsub job.sh
done

