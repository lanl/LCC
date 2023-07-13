Building LCC documentation
===========================

The folder (:code:`src/docs`:) contains all the documentation relevant to both users and
developpers.

Prerequisites
---------------

  - `[pdflatex] <http://pdftex.org>`_ Latex GNU compiler. pdfTeX is an extension of TeX which can produce PDF directly from TeX source, as well as original DVI files. pdfTeX incorporates the e-TeX extensions.  

  - `[doxygen] <https://www.doxygen.nl/index.html>`_  Doxygen  is a documentation system for C++, C, Java, Objective-C, IDL (Corba and Microsoft flavors) and to some extent PHP, C#, and D.
  
  - `[sphinx] <https://www.sphinx-doc.org/en/master/usage/quickstart.html>`_ Sphinx is a documentation generator or a tool that translates a set of plain text source files into various output formats, automatically producing cross-references, indices, etc. That is, if you have a directory containing a bunch of reStructuredText or Markdown documents, Sphinx can generate a series of HTML files, a PDF file (via LaTeX), man pages and much more.

  - Any pdf viewer.

  - Any web browser. 


These programs can be installed as follows::

  sudo apt-get install pdflatex 
  sudo apt-get install doxygen 
  sudo apt-get install dot2tex
  sudo apt-get install python3-sphinx
  pip3 install PSphinxTheme
  pip3 install recommonmark

Build the full documentation
------------------------------
 
This will build all three types of docs (Sphinx, Doxygen, and latex)::

  make  

The documentation that is build with Sphinx can be tested as follows::

  firefox lcc.html

The file can be explored using any web browser.  

One can also build any of the documentations separatly. For example, to build 
the Sphinx documentation, we can do::

  make sphinx 

Documenting 
------------

In order to add a documentation using Sphinx follow these steps: 
  1) make a file with a proper name under :code:`./sphinx-src/source/`. For example: :code:`MYPAGE.md`. 
  2) Add the documentation inside the file using "markdown" syntax. 
  3) Modify the file in :code:`./sphinx-src/source/index.txt` to include the documentation.

After modyfing this file, recompile Sphinx by typing :code:`make sphinx`.
