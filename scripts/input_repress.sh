#~/bin/bash
echo "Running Repressilator Script:"
for i in {1..1}; do
	echo $i
	./bin/prog -s repressilator -t 0:.01:10 -i 0:.01:10 -u 216,0.216,2,0.2 -o 150,.5,1,.5 -y 5,1,10,1,15,1 -k 1 -p 12 > results/repress_results/repress$i.txt
done
echo "Finished Running Repressilator Script"
