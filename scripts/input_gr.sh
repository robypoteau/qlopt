#~/bin/bash
echo "Running General Repressilator Script:"
for i in {1..2}; do
	echo $i
	./bin/main \
	-s general_repressilator \
	-t 0:.01:10 \
	-i 0:.01:10 \
	-u  216,0.216,2,0.2,0.2,216,0.216,2,0.2,0.2,216,0.216,2,0.2,0.2 \
	-o 150,.5,.5,.5,.5,150,.5,.5,.5,.5,,150,.5,.5,.5,.5 \
	-y 5,1,10,1,15,1 \
	-k 4 \
	-p 16 \
	-r \
	> results/gr_results/gr$i.txt
done
echo "Finished Running General Repressilator Script"
#	-n $i \
