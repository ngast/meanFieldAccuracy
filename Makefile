all: simulator simulations paper

paper:
	pdflatex mf_rate_convergence.tex
	bibtex mf_rate_convergence.aux
	pdflatex mf_rate_convergence.tex
	pdflatex mf_rate_convergence.tex
	pdflatex mf_rate_convergence.tex

simulator:
	cd jsq2_simulate/ && g++ -W -Wall -O2 two_choice_simulate.cc -o simulator

simulations:
	echo "this might take a while"
	cd jsq2_simulate/ && echo "TBD"
	cd simu && echo "TBD2"
