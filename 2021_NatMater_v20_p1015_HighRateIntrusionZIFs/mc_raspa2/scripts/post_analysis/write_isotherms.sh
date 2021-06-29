rm gcmc_results.txt
rm gcmcdesorption_results.txt
fws=(ZIF8_opengate ZIF8_closedgate)
ffs=(qff-dreiding-mbis-tip4p)
guests=(water argon N2)
shopt -s extglob

for ff in ${ffs[@]}
do
    for fw in ${fws[@]}
    do 
        for guest in ${guests[@]}
        do 
            cat ${fw}/${ff}/${guest}/t*-p+([0-9])/Output/System_0/*.summary >> gcmc_results.txt 2> /dev/null
            cat ${fw}/${ff}/${guest}/t*-p+([0-9])-desorption/Output/System_0/*.summary >> gcmcdesorption_results.txt 2> /dev/null
            grep "Conversion factor" ${fw}/${ff}/${guest}/t0-p0/Output/System_0/*.data > conversion_${fw}_${guest}.txt
        done
    done
done
