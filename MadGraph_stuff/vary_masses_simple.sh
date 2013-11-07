#!/bin/bash

# turn off the automatic htmal open
sed -i -e 's/automatic_html_opening = True/automatic_html_opening = False/g' Cards/me5_configuration.txt

#change the number of events to 3
sed -i -E -e 's/[0-9]+ = nevents/500000 = nevents/g' Cards/run_card.dat

#change the CM energy to 8Tev
sed -i -E -e 's/[0-9]+.*=.*ebeam1/4000     = ebeam1/g' Cards/run_card.dat
sed -i -E -e 's/[0-9]+.*=.*ebeam2/4000     = ebeam2/g' Cards/run_card.dat

# change the Higgs mass to 125 Gev
sed -i -E -e 's/25 [0-9]+?\.?[0-9]+?.* .*MH/25 1.250000e+02 \# MH/g' Cards/param_card.dat

#make results directory
mkdir Graviton_Parton

# change the Graviton/mass  550 
for i in 250 550 750 1000 1500 2000 2500 3000; do

sed -i -E -e "s/111 [0-9]+?\.?[0-9]+?.* \# MH02/111 $i \# MH02/g" Cards/param_card.dat
sed -i -E -e "s/39 [0-9]+?\.?[0-9]+?.* \# MGr/39 $i \# MGr/g" Cards/param_card.dat

#Graviton or graviton?
grep -E '111 [0-9]+?\.?[0-9]+?.* \# MH02' Cards/param_card.dat
var1=$?
grep -E '39 [0-9]+?\.?[0-9]+?.* \# MGr' Cards/param_card.dat
var2=$?
echo $var1 $var2

#generate the events
./bin/generate_events 0 M_$i run1

gunzip Events/M_$i/unweighted_events.lhe.gz 
mv Events/M_$i/unweighted_events.lhe Graviton_Parton/MGraviton_$i.lhe

#change the events ID
#on the header
#pcregrep -M ' 0\n\<\/init\>' Graviton_Parton/MGraviton_250.lhe
# nao funciona mudar a mao ...
#sed -i -e 's/\(0\)\n\<\/init\>/100\n\<\/init\>/g' Graviton_Parton/MGraviton_$i.lhe
#change the events
sed -i -e 's/ 8   0/ 8   100/g' Graviton_Parton/MGraviton_$i.lhe
sed -i -e 's/ 9   0/ 9   100/g' Graviton_Parton/MGraviton_$i.lhe

done
