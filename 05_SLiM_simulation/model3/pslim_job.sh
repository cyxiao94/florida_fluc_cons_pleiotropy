sleep 4h
echo "job start"
../pslim.sh -r -X ../slim.tar.gz -i 03_selection.slim -n 1000 -p 0parameters -f 0files/02pleiotropic_effects,0files/01vcffile -o out -M 10240 -T ./tempf 2>&1 | tee > output.log


echo "job finished"

