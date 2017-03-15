g++ -W -Wall -O3 two_choice_simulate.cc

number_of_simulation() {
    N=$1
    C=`expr 500 + $N \* $N / 100`
    echo $C
}

launch_for_rho (){
    rho=$1
    finished=()
    for i in `seq 1 100`; do
	for N in 3 5 10 20 30 40 50 100 200 300 400 500 1000; do
	    if ! [ ${finished[$N]} ]; then 
		fileName=results/exp_${rho}_${N}
		if [ -e $fileName ]; then A=`wc -l $fileName`; else A=0; fi
		simulations_performed=`echo $A | sed 's/ .*//'`
		simulations_todo=`number_of_simulation $N`
		if [ $simulations_performed -le $simulations_todo ] ;
		then
		    # We launch 400 additional simulations 
		    echo "Simulations for $N 0.${rho} not done    ($simulations_performed < $simulations_todo)"
		    #{ time  ./a.out r0.${rho} N${N} e500 >> ${fileName};} 2>&1 | grep real
		else
		    finished[$N]=1
		    echo "N=$N rho=0.$rho done     ($simulations_performed < $simulations_todo)";
		fi
	    fi
	done
    done
}

launch_for_rho_q05 (){
    echo "PID=$$"
    rho=$1
    for i in `seq 1 100`; do
	for N in 10 20 50 100 200 500 1000; do
	    #for rho in 50 60 70 75 80 85 90 95 99; do
	    fileName=results/exp_q05_${rho}_${N}
	    if [ -e $fileName ]; then A=`wc -l $fileName`; else A=0; fi
	    B=`echo $A | sed 's/ .*//'`
	    C=`expr 100 + $N \* $N / 10`
	    if [ $B -le $C ] ;
	    then
		echo "$N 0.$rho not done \t ($B < $C)"
		{ time  ./a.out r0.${rho} q0.05 N${N} e400 >> ${fileName};} 2>&1 | grep real
	    else
		echo "$N 0.$rho done \t ($B >= $C)";
	    fi
	done
    done
}

for rho in 50 75 85 ; do
    launch_for_rho $rho
    #launch_for_rho_q05 $rho &
done
