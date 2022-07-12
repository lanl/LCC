# Makefile for the LCC documentation

DIR=$(shell pwd)

all:
	(make doxy)	

clean:
	rm -f *.backup *.log *.bbl *.aux *.blg *.out
	rm -f MyBib.bib~ lcc.tex~ Makefile~
	rm -f tmp* 
	rm -rf ./doxy-out

doxy:
	(doxygen ./doxy-src/Doxyfile)	
	(cd $(DIR)/doxy-out/latex ; make)	
	(cp $(DIR)/doxy-out/latex/refman.pdf lcc_refman_devel.pdf)
	(doxygen ./doxy-src/Doxyfile_user)	
	(cd $(DIR)/doxy-out/latex ; make)	
	(cp $(DIR)/doxy-out/latex/refman.pdf lcc_refman_user.pdf)

         
