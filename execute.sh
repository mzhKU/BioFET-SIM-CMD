#!/bin/bash
input="1ssk-07.40-pos.bfs"
for i in {100..300}
do
    python -c "print $i/100.0"
    #a=`python -c "print $i/100.0"`
    #python bio_run.py --set L_d $a $input
    #python bio_run.py --calc $input
done
