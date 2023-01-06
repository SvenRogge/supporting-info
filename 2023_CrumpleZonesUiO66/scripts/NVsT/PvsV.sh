#!/bin/sh

module load yaff

#for j in CN10-5_typexx  CN10-5_typexy  CN10_typeyy  CN11-5_type0  CN11_type5  CN12
#for j in CN10-5_typexx_av  CN10-5_typexy_av  CN10_typeyy_av  CN11_type5_av
#for j in CN12_av
for j in CN12_av
#for j in CN10-5_typexx_av CN10-5_typexy_av
do
    cd ../${j}
    for (( i=5000 ; i < 9800; i += 25 ))
    do
        if [ -d ${i} ]; then
            cd ${i}
            cp ../../scripts/PvsV.py .
            python PvsV.py
            rm PvsV.py
            cd ..
        fi
    done
done
