#~/bin/bash
echo "Running General Repressilator Script:"
for i in {1..1}; do
	echo $i
	./bin/prog -s general_repressilator -t 0:.01:10 -i 0:.01:10 -u  216,0.216,2,0.2,0.2,216,0.216,2,0.2,0.2,216,0.216,2,0.2,0.2 -o 150,.5,.5,.5,.5,150,.5,.5,.5,.5,,150,.5,.5,.5,.5 -y 5,1,10,1,15,1 -k 1 -p 12 > results/gr_results/gr$i.txt
done
echo "Finished Running General Repressilator Script"	