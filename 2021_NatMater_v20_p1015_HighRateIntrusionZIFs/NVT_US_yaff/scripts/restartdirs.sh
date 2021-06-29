LC_NUMERIC=C

cv_start=-11.5
cv_stop=11.5
cv_step=0.5

N_pore1=5
N_pore2=10
dir_prefix="${N_pore1}H2O_${N_pore2}H2O"

fn_prefix="traj"
run=$( ls ../${dir_prefix}/0.0 | grep "${fn_prefix}" | sort -V | tail -1 )
run=$(( ${run:$(( ${#fn_prefix} + 1 )):-3} + 1 ))
run_old=$(( ${run} - 1 ))
run_older=$(( ${run} - 2 ))

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

    cp inputfiles/ymd_ghost_restart.py ${dir_prefix}/${j}/ymd_ghost_restart.py
    cp inputfiles/job_restart.sh ${dir_prefix}/${j}

    cd ${dir_prefix}/${j}

    job_id="ZIF8_US_${N_pore1}H2O_${N_pore2}H2O_300K_${j}_restart"

    sed -i "s/job_id/${job_id}/g" job_restart.sh
    sed -i "s/restart_0/restart_${run_old}/g" job_restart.sh

    sed -i "s/DUMMY_NPORE1/${N_pore1}/g" ymd_ghost_restart.py
    sed -i "s/DUMMY_NPORE2/${N_pore2}/g" ymd_ghost_restart.py
    sed -i "s/DUMMY_CV0\*angstrom/${i}\*angstrom/g" ymd_ghost_restart.py
    sed -i "s/DUMMY_K\*kjmol/${K}\*kjmol/g" ymd_ghost_restart.py
    sed -i "s/traj_1.h5/traj_${run}.h5/g" ymd_ghost_restart.py
    sed -i "s/restart_1.h5/restart_${run}.h5/g" ymd_ghost_restart.py
    sed -i "s/restart_0.h5/restart_${run_old}.h5/g" ymd_ghost_restart.py

    qsub job_restart.sh

    cd ../../inputfiles

done
