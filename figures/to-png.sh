#!/bin/bash

pdftoppm -r 450 figure-crop.pdf | pnmtopng > figure.png
