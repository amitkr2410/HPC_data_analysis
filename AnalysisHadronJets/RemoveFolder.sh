#!/bin/bash                                                                                                                                                   

for (( Trial=1; Trial<13; Trial++ ))
do
echo "removing code"$Trial
rm -rf "Code$Trial" &
done
wait
for (( Trial=13; Trial<22; Trial++ ))
do
echo "removing code"$Trial
rm -rf "Code$Trial" &
done
wait

for (( Trial=22; Trial<30; Trial++ ))
do
echo "removing code"$Trial
rm -rf "Code$Trial" &
done
