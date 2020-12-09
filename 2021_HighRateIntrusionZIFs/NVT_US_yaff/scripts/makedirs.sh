LC_NUMERIC=C

cv_start=-11.5
cv_stop=11.5
cv_step=0.5

N_pore1=5
N_pore2=10
dir_prefix="${N_pore1}H2O_${N_pore2}H2O"

run=$( ls .. | grep ${dir_prefix} | sort -V | tail -1 ) 
run=$(( ${run:$(( ${#dir_prefix} + 1 ))} + 1 ))

K=25
temp=300

for i in $(seq ${cv_start} ${cv_step} ${cv_stop})
do

    if (( $( echo "${i} < 0" | bc -l ) == 1 )) 
    then
       j=$( echo "${i}*(-1)" | bc )
       j="m$( printf "%.1f\n" ${j} )"
    else
       j="$( printf "%.1f\n" ${i} )"
    fi

    cd ..

    mkdir -p ${dir_prefix}_${run}/${j}
    cp inputfiles/ymd_ghost.py ${dir_prefix}_${run}/${j}/ymd_ghost.py
    cp structures/init_${N_pore1}H2O_${N_pore2}H2O_US${i}.chk ${dir_prefix}_${run}/${j}/init.chk
    cp inputfiles/job.sh ${dir_prefix}_${run}/${j}

    cd ${dir_prefix}_${run}/${j}

    job_id="ZIF8_US_${N_pore1}H2O_${N_pore2}H2O_${run}_300K_${j}"

    sed -i "s/job_id/${job_id}/g" job.sh
    #sed -i "s/job_id/${job_id}/g" job_restart.sh

    sed -i "s/temp = temp\*kelvin/temp = ${temp}\*kelvin/g" ymd_ghost.py
    sed -i "s/DUMMY_NPORE1/${N_pore1}/g" ymd_ghost.py
    sed -i "s/DUMMY_NPORE2/${N_pore2}/g" ymd_ghost.py
    sed -i "s/DUMMY_CV0\*angstrom/${i}\*angstrom/g" ymd_ghost.py
    sed -i "s/DUMMY_K\*kjmol/${K}\*kjmol/g" ymd_ghost.py

    qsub job.sh

    cd ../../inputfiles

done
