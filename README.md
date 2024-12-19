# LIQA2
LIQA2 is a computational tool for detecting and quantifying isoform expression (including novel isoforms) from long read sequencing. The tool uses a meta-regression model to quantify differential alternative splicing and incorporates a targeted approach to correct for 3' and 5' biases.

## Installation

## Pipeline
The main function for running the LIQA2 pipeline is `main_preprocessing.py`. There are 4 steps, marked by the --task [step] command.
### 1. Prepare annotation file
--task annotation
### 2. Generate compatible matrix
--task compatible matrix
### 3. Summarize novel gene annotations
--task summary
### 4. Generate count matrix
--task count matrix
