#!/bin/bash
# BioFET-SIM CMD range execution
input=$1
for i in {100..300}
do
    a=`python -c "print $i/100.0"`
    # Prevent printing print statements of bio_run.py script
    b=`python bio_run.py --set L_d $a $input`
    # Assemble result string, delete ('-d') whitespace (' ')
    res=`python bio_run.py --calc $input|grep -i "sensitivity"|cut -d ":" -f 2|tr -d " "`
    echo $a $res
done
