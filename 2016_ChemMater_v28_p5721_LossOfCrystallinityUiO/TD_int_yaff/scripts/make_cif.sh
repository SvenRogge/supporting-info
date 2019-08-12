module load con3f/version_CIF2-intel-2015a-Python-2.7.10

for i in {8000..9700..10}
do
    cd ${i}
    c3f.py convert traj_conc12.h5:1000:-1:1 traj_conc12.cif --average
    cd ..
done
