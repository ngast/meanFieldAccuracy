foo() {
    rho=$1
    echo $rho
}

for rho in .80 .95 .99; do
    foo $rho 
done
