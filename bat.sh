for n in 40 
do
	for p in 0.15
	do
		for b in 0.2 0.3
		do
		   echo "numnodes" $n > parfile.txt
		   echo "gtype er"  >> parfile.txt
		   echo "p " $p >> parfile.txt
		   echo "edge_density 0.2"   >> parfile.txt
		   echo "ws_param 0.1"  >> parfile.txt
		   echo "itemcosts 1"  >> parfile.txt
		   echo "budgetprop" $b >> parfile.txt
		   echo "dbug 0"  >> parfile.txt
		   echo "test 0"  >> parfile.txt
		   echo "ccut 1" >> parfile.txt
		   echo "os m" >> parfile.txt
		   echo "file 2"  >> parfile.txt
		   for i in {1..10}
		   do
		     python BiCNDP.py $i >> result_new.log
		   done
		done
	done
done
	   	  
