#~/bin/bash
echo "Running Toggle Switch Script:"
for i in {1..1}; do
	echo $i
	./bin/prog -s repressilator -t 0:.005:5 -i 0:.005:5 -u  216,0.216,2,0.2 -o .5,.5,.5,.5 -y 5,2,5,2,5,2 -k 1 -p 12 -r > results/ts_results/ts$i.txt
done
echo "Finished Running Toggle Switch Script"	