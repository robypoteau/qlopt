#~/bin/bash
echo "Running Bistable Script:"
for i in {1..5}; do
	echo $i
	./bin/prog -s bistable_switch -t 0:.005:3 -i 0:.005:3 -u 150,3.2,2 -o 125,2,1 -y 25,25 -k 10 -p 5 -n $i -b -r > results/bs_results/bs$i.txt
done
echo "Finished Running BS Script"
