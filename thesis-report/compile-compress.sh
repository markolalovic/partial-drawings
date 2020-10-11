#!/bin/bash
pdflatex -output-directory=./build 0-main.tex
bibtex ./build/0-main.aux 
pdflatex -output-directory=./build 0-main.tex
pdflatex -output-directory=./build 0-main.tex

# compress to screen size (72dpi)
gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=/screen -dNOPAUSE -dQUIET -dBATCH -sOutputFile=../thesis-report.pdf ./build/0-main.pdf
