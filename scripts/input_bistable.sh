#~/bin/bash
echo "Running Bistable Script:"
for i in {11..12}; do
	echo $i
	./bin/prog -s bistable_switch -t 0:.001:10 -i 0:.001:7 -u 150,3.2,2 -o 135,2.2,1.5 -y 25,25 -k 10 -p 15 -n $i -b -d 2 > results/bs_results/bs$i.txt
done
echo "Finished Running BS Script"
#2-4 1 div
#5-8 2 div
