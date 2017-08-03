#~/bin/bash
echo "Running Repressilator Script:"
for i in {1..1}; do
	echo $i
	./bin/main \
	-s repressilator \
	-t 0:.01:10 \
	-i 0:.01:10 \
	-u 216,0.216,2,0.2 \
	-o 150,.5,1,.5 \
	-y 5,1,10,1,15,1 \
	-k 100 \
	-p 16 \
	-n $i \
	> results/repress_results/repress$i.txt
done
echo "Finished Running Repressilator Script"
