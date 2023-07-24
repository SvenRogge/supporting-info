#!/bin/sh

vol_start=7000
vol_stop=9800
vol_step=25
path=./structures/


all_found=true
for (( i=${vol_start}; i<=${vol_stop}; i=i+${vol_step} ))
do
    if [ ! -f ${path}/init_${i}.chk ]
    then
        echo "File with volume ${i} not found"
        all_found=false
    fi
done

if $all_found
then
    echo "All files in the range ${vol_start}..${vol_stop}..${vol_step} were found"
fi
