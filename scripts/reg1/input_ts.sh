#~/bin/bash
echo "Running Toggle Switch Script:"
for i in {1..1}; do
	echo $i
	./bin/prog -s toggle_switch -t 0:.005:6 -i 0:.005:6 -u 1,1,10,2,1,1,1,1,10,2,1,1 -o .5,.5,5,.5,.5,.5,.5,.5,5,.5,.5,.5 -y 5,2,5,2 -k 1 -p 12 -r > results/ts_results/ts$i.txt
done
echo "Finished Running Toggle Switch Script"	