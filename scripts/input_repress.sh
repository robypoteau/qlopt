#~/bin/bash
echo "Running Toggle Switch Script:"
for i in {1..1}; do
	echo $i
	./bin/prog -s repressilator -t 0:.005:10 -i 0:.005:10 -u 216,0.216,2,0.2 -o 150,.5,.5,.5 -y 5,2,5,2,5,2 -k 1 -p 12 > results/repress_results/repress$i.txt
done
echo "Finished Running Toggle Switch Script"	