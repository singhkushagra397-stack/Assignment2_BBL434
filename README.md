# Assignment2_BBL434
# Global DNA Sequence Alignment Tool

## Overview

This project is a simple Python program that performs **global alignment**
of two DNA sequences using **affine gap penalties**.

The program reads sequences from FASTA files, applies user-defined
scoring values, and outputs the best alignment along with the score.

---

## What this program does

- Reads two nucleotide sequences from FASTA files
- Uses match, mismatch, gap opening, and gap extension scores
- Computes the best global alignment
- Displays the alignment and final score

---

## Requirements

- Python 3.x
- NumPy

If NumPy is not installed:

```bash
pip install numpy
