# Make file of plots.
Assignment4.pdf: Assignment4.tex plot1.pdf plot2.pdf plot3.pdf plot4.pdf plot5.pdf plot6.pdf plot7.pdf plot8.pdf plot9.pdf plot10.pdf
	pdflatex Assignment4.tex
	pdflatex Assignment4.tex

# Make plots.
#.PHONY : plots
#plots : plot1.pdf plot2.pdf plot3.pdf plot4.pdf plot5.pdf plot6.pdf plot7.pdf plot8.pdf plot9.pdf plot10.pdf
plot1.pdf : as3.py
	python $< 1 $@
plot2.pdf : as3.py
	python $< 2 $@
plot3.pdf : as3.py
	python $< 3 $@
plot4.pdf : as3.py
	python $< 4 $@
plot5.pdf : as3.py
	python $< 5 $@
plot6.pdf : as3.py
	python $< 6 $@
plot7.pdf : as3.py
	python $< 7 $@
plot8.pdf : as3.py
	python $< 8 $@
plot9.pdf : as3.py
	python $< 9 $@
plot10.pdf : as3.py
	python $< 10 $@

#.PHONY : clean
#clean :
#	rm -f *.dat
#	rm -f results.txt
