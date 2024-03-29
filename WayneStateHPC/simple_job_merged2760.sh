#!/bin/bash
#      0  1  2  3  4  5  6  7   8   9   10  11  12  13  14  15  16  17  18  19  20  21  22  23   24   25   26   27   28   29   30   31   32   
#33  34  35    36  37    38  39    40   41  42   43    44  45   46   47   48   49   50   51   52   53    54    55    56    57    58    59   60
#    61    62    63    64    65    
pTHat=(1  2  3  4  5  7  9  11  13  15  17  20  25  30  35  40  45  50  55  60  70  80  90  100  110  120  130  140  150  160  170  180  190  200  210  220  230  240  250  260  270  280  290  300  350  400  450  500  550  600  700  800  900  1000  1380 )
#pTHat=(1  2  3  4  5  7  9  11  13  15  17  20  25  30  35  40  45  50  55  60  70  80  90  100  110  120  130  140  150  160  170  180  190  200  210  220  230  240  250  260  270  280  290  300  350  400  450  500  550  600  700  800  900  1000  1100  1200  1300  1400  1500  1600  1700  1800  1900  2000  2200  2400  2510)

DirInput=$1
DirOutput=$2
HadronORParton=$3
Index1=$4
Index2=$((1+${Index1}))

cat  ${DirInput}/Jetscape${HadronORParton}ListBin${pTHat[${Index1}]}_${pTHat[${Index2}]}_Run*.out >  ${DirOutput}/Jetscape${HadronORParton}ListBin${pTHat[${Index1}]}_${pTHat[${Index2}]}.out
