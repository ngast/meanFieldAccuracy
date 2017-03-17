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
	echo "\n### this might take a while ###\n"
	cd jsq2_simulate/ && jupyter-3.5 nbconvert --execute  "Rate of convergence two choice.ipynb"
	cd simu && jupyter-3.5 nbconvert --execute --ExecutePreprocessor.timeout=-1 mean_field_error_simulation_1D.ipynb

