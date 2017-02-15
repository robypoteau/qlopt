#~/bin/bash
echo "Running General Lotka Volterra Script:"
for i in {1..1}; do
	echo $i
	./bin/prog -s general_lv -t 0:.001:2 -i 0:.001:2 -u .48,0,-.026,-.53,.026,0 -o .48,0,-.026,-.93,.026,0 -y 35,4 -k 1  -p 7 -d 1 > results/glv_results/glv$i.txt
done
echo "Finished Running General Lotka Volterra Script"
