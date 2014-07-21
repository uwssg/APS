make mcmcBrute

for ((i=0 ; i<100 ; i++)) do
./mcmcBrute -1 5 3
sleep 2
done

for ((i=0 ; i<100 ; i++)) do
./mcmcBrute -1 5 2
sleep 2
done

for ((i=0 ; i<100 ; i++)) do
./mcmcBrute -1 5 4
sleep 2
done
