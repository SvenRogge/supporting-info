%chk=gaussian.chk
%nproc=1
%mem=1GB
#p B3LYP/6-311G(d,p) freq=noraman int(grid=ultrafine) maxdisk=50GB

title

0 1
 S     0.000000     0.000000     0.000000
 F     0.000000     0.000000     1.603354
 F     1.603354     0.000000     0.000000
 F     0.000000     0.000000    -1.603354
 F    -1.603354     0.000000     0.000000
 F     0.000000    -1.603354     0.000000
 F     0.000000     1.603354     0.000000

