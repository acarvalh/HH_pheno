#!/bin/bash

# arrays on parameters 
declare -a MR=(270 300 350 400 450 500 550 600 650 700 750 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 2000 2100 2200 2300 2400 2500 2600 2700 2800 2900 3000) 
declare -a ImFmeio=(0. 0. 0. 0.595708 1.10968 1.40622 1.55 1.59889 1.59138 1.5515 1.5515 1.42767 1.28839 1.1555 1.03591 0.930789 0.839232 0.759665 0.690437 0.630027 0.577111 0.530565 0.489445 0.452963 0.42046 0.391386 0.365279 0.341751 0.320474 0.30117 0.283601 0.267565 0.252888 0.239419)
declare -a ReFmeio=(-1.58565 -1.67995 -2. -2.31658 -2.1282 -1.81233 -1.49395 -1.21145 -0.972651 -0.774884 -0.774884 -0.479396 -0.280785 -0.146203 -0.0540262 0.00968438 0.053976 0.0848092 0.10618 0.120813 0.130599 0.136867 0.140566 0.142379 0.142806 0.142216 0.140881 0.139005 0.136741 0.134206 0.131488 0.128652 0.125749 0.122819) 


###############################################
# turn off the automatic htmal open
sed -i -e 's/automatic_html_opening = True/automatic_html_opening = False/g' Cards/me5_configuration.txt

#change the number of events to 3
sed -i -E -e 's/[0-9]+ = nevents/20000 = nevents/g' Cards/run_card.dat

#change the minimum distance between jets
sed -i -E -e 's/ [0-9]+?\.?[0-9]+?.* \= drjj    \! min distance between jets/0.0001 \= drjj    \! min distance between jets/g' Cards/run_card.dat
sed -i -E -e 's/ [0-9]+?\.?[0-9]+?.* \= draa    \! min distance between gammas/ 0.0001 \= draa    \! min distance between gammas/g' Cards/run_card.dat
sed -i -E -e 's/ [0-9]+?\.?[0-9]+?.* \= draj    \! min distance between gamma and jet/ 0.0001 \= draj    \! min distance between gamma and jet/g' Cards/run_card.dat

#change the CM energy to 2Tev
sed -i -E -e 's/[0-9]+.*=.*ebeam1/4000     = ebeam1/g' Cards/run_card.dat
sed -i -E -e 's/[0-9]+.*=.*ebeam2/4000     = ebeam2/g' Cards/run_card.dat

# change the Higgs mass to 125 Gev
#sed -i -E -e 's/25 [0-9]+?\.?[0-9]+?.* .*MH/25 1.25000e+02 \# MH/g' Cards/param_card.dat

# change LR
#sed -i -E -e 's/25 [0-9]+?\.?[0-9]+?.* .*MH/25 1.25000e+02 \# MH/g' Cards/param_card.dat
################################################

mkdir Radion_bbtata_LR3tevLHC8

#vary on xi
#for xi in 0.166667 1 1.1;
#do

  #change the parameters       9000004 1.000000e-06 # xi
 # sed -i -E -e "s/ 9000004 [0-9]+?\.?[0-9]+?.* \# xi/ 9000004 $xi \# xi/g" Cards/param_card.dat

  #vary on the masses
  for (( i = 0 ; i <${#MR[@]} ; i++ )); do
# ${#MR[@]} ; i++ )); do
    #echo "${aa[$i]}"  "${bb[$i]}"  "${b3[$i]}" "${d3m5[$i]}"  "${d4m5[$i]}"

    #change the parameters
    sed -i -E -e "s/ 9000001 [0-9]+?\.?[0-9]+?.* \# ImFmeio/ 9000001 ${ImFmeio[$i]} \# ImFmeio/g" Cards/param_card.dat
    sed -i -E -e "s/ 9000003 [0-9]+?\.?[0-9]+?.* \# ReFmeio/ 9000003 ${ReFmeio[$i]} \# ReFmeio/g" Cards/param_card.dat
    sed -i -E -e "s/ 35 [0-9]+?\.?[0-9]+?.* \# MH02/ 35 ${MR[$i]} \# MH02/g" Cards/param_card.dat

    #generate the events
    ./bin/generate_events 0 MR_${MR[$i]}_on run1

    gunzip Events/MR_${MR[$i]}_on/unweighted_events.lhe.gz 
    mv Events/MR_${MR[$i]}_on/unweighted_events.lhe Radion_bbtata_LR3tevLHC8/MR_${MR[$i]}_on.lhe


    #take the CX
    #echo "a" "${aa[$i]}" "$(grep -E '\#  Integrated weight \(pb\)  \:  [0-9]+?.*?' Events/a_${aa[$i]}_m4/a_${aa[$i]}_m4.lhe)"
    #echo "a" "${aa[$i]}" "$(grep -E '\#  Integrated weight \(pb\)  \:  [0-9]+?.*?' Events/a01/a01.lhe)" >> teste.txt

    sed -i -e 's/ 8   0/ 8   100/g' Radion_bbtata_LR3tevLHC8/MR_${MR[$i]}_on.lhe
    sed -i -e 's/ 9   0/ 9   100/g' Radion_bbtata_LR3tevLHC8/MR_${MR[$i]}_on.lhe

    #take the CX
    echo "MR" "${MR[$i]}" "$(grep -E '\#  Integrated weight \(pb\)  \:  ' Radion_bbtata_LR3tevLHC8/MR_${MR[$i]}_on.lhe)" >> Radion_bbtata_LR3tevLHC8/CX_LR3tev_LHC100.txt

  done
#done

