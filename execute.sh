#!/bin/bash
<<<<<<< HEAD
input="1ssk-11.00-neg.bfs"
for i in {1..300}
=======
input="1ssk-07.40-pos.bfs"
for i in {100..300}
>>>>>>> 34eb046b3f79295d7a86b24aa23e11411c7d4bf5
do
    python -c "print $i/100.0"
    #a=`python -c "print $i/100.0"`
    #python bio_run.py --set L_d $a $input
    #python bio_run.py --calc $input
done
