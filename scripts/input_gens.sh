#~/bin/bash
echo "Running Toggle Switch Script:"
for i in {1..1}; do
	echo $i
	./bin/prog -s toggle_switch -t 0:.05:12 -i 0:.05:12 -u 1,1,10,2,1,1,1,1,10,2,1,1 -o .5,.5,5,.5,.5,.5,.5,.5,5,.5,.5,.5 -y 5,2,5,2 -k 5 -p 4 -n $i -r -b > results/ts_results/ts$i.txt
done
echo "Finished Running Toggle Switch Script"
