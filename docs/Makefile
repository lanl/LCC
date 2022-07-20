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

sphinx:
	(cp ../README.md ./sphinx-src/source/README-main.md)
	(cp ./README.md ./sphinx-src/source/)
	(cd sphinx-src/source; mkdir _static; sphinx-build -b html . ../../sphinx-out)
	(cp -r figures sphinx-out/docs/)
	(rm lcc_user.html; ln -s ./sphinx-out/index.html ./lcc_user.html)

