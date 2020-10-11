#!/bin/bash
pdflatex -output-directory=./build seminar-slides.tex
mv ./build/seminar-slides.pdf ..
