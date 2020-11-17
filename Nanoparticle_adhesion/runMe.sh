#!/bin/bash

module load anaconda3/personal
module load matlab

l=9
m=4
r=2
s=1
p=1

mkdir L${l}M${m}R${r}
cp * L${l}M${m}R${r}
cd L${l}M${m}R${r}

repstr1="REPLACEL"
repstr2="REPLACEM"
repstr3="REPLACER"
repstr4="REPLACES"
repstr5="REPLACEP"

sed -i -e "s/$repstr1/$l/g" -e "s/$repstr2/$m/g" -e "s/$repstr3/$r/g" -e "s/$repstr4/$s/g" -e "s/$repstr5/$p/g"  doMinimize.py 

sed -i -e "s/$repstr1/$l/g" -e "s/$repstr2/$m/g" -e "s/$repstr3/$r/g" -e "s/$repstr4/$s/g" -e "s/$repstr5/$p/g"  Characterize.m 

matlab -nodisplay -nosplash -nodesktop -r "run('Characterize.m');exit;" | tail -n +11
python doMinimize.py
qsub run.pbs
rm runMe.sh


