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
	(cp -r doxy-out/html ./sphinx-out/_static/doxy)

sphinx:
	(cp ../README.rst ./sphinx-src/source/README-main.rst)
	(cp ./README-doc.rst ./sphinx-src/source/README-doc.rst)
	(cd sphinx-src/source; sphinx-build -b html . ../../sphinx-out)
	(cd sphinx-src/source; sphinx-build -b latex . ../../sphinx-out)
	(cd ./sphinx-out ; make ; cp lcc.pdf ../lcc_user.pdf)
	(rm lcc_user.html; ln -s ./sphinx-out/index.html ./lcc_user.html)

