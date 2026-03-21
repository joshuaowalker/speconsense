.PHONY: figures paper clean

figures:
	python compute_critical_rates.py

paper: figures
	pdflatex variant_significance.tex

clean:
	rm -f *.aux *.log *.out figures/*.pdf
