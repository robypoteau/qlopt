#~/bin/bash
echo "Running Bistable Script:"
for i in {6..9}; do
	echo $i
	./bin/prog -s bistable_switch -t 0:.1:3 -i 0:.1:3 -u 150,3.2,2 -o 135,2.2,1.5 -y 25,25 -k 10 -p 15 -n $i -b -d 1 > results/bs_results/bs$i.txt
done
echo "Finished Running BS Script"
#2-4 1 div
#5-8 2 div
