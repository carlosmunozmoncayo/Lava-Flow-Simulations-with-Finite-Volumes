#!/bin/bash

latexmk -pdf -bibtex first_report
open first_report.pdf
