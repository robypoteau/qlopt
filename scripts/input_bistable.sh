#~/bin/bash
echo "Running Bistable Script:"
for i in {0..0}; do
	echo $i
	./bin/main -s bistable_switch -i 0:1:3 -u 150,3.2,2.0 -o 125,2.2,1.0 -y 14,13 -k 1 -p 5 -d 1 -r 0 > results/bs_results/bs$i.txt
done
echo "Finished Running BS Script"
#2-4 1 div
#5-8 2 div
# -n $i
