# LIQA2
LIQA2 is a computational tool for detecting and quantifying isoform expression (including novel isoforms) from long read sequencing. The tool uses a meta-regression model to quantify differential alternative splicing and incorporates a targeted approach to correct for 3' and 5' biases.

## Installation
```
git clone https://github.com/WGLab/LIQA2.git
cd LIQA2
```

## Pipeline
The main function for running the LIQA2 pipeline is `main_preprocessing.py`. There are 4 steps, marked by the --task [step] command.

  
### 1. Prepare annotation file
The first step of the LIQA2 pipeline involves generaing annotation files from the reference genome and indexed/tagged BAM files. This is done with `--task annotation`. 

Arguments include:
- `--target`: path to output location(s)
- `--bam`: path to bam location(s)
- annotation mode (see below)
- `--reference`: path to GTF annotation file
- `reference_pkl`: path to picke gene annotation file generated by LIQA2 (geneStructureInformation.pkl)
- `--workers`: threads for parallel processing
- `--coverage_threshod_exon`: coverage threshold for exon discovery as a percent of maximum coverage. Default is 0.02, larger values are more conservative

Below is a sample command that can be run.

```
python3 src/main_preprocessing.py \
--task annotation \
--target path/to/output/folder/of/sample1 path/to/output/folder/of/sample2 \
--bam path/to/bam/file/or/bamfolder/sample1 path/to/bam/file/or/bamfolder/sample2 \
--update_gtf \
--reference path/to/reference/genes.gtf \
--reference_pkl data/geneStructureInformation.pkl \
--workers 30
```
There are 3 different modes of annotation.
1. Annotation-Only: 
2. Enhanced-Annotation:
3. Annotation-Free: 

### 2. Generate compatible matrix
This step creates a compatible matrix that aligns reads to existing gene isoforms and identifies novel isoforms based on previously generated annotation files.
Arguments include:
- `--target`: path to output location(s)
- `--bam`: path to bam location(s)
- `--reference`: path to GTF annotation file

Below is a sample command that can be run.
```
python3 src/main_preprocessing.py \
--task 'compatible matrix' \
--target path/to/output/folder/of/sample1 path/to/output/folder/of/sample2 \
--bam path/to/bam/file/or/bamfolder/sample1 path/to/bam/file/or/bamfolder/sample2 \
--reference path/to/reference/genes.gtf
```
### 3. Summarize novel gene annotations
After annotation and compatible matrix generation, the LIQA2 pipeline generates a new summary annotation file that also includes novel isoforms:
```
python3 src/main_preprocessing.py \
--task 'summary' \
--target path/to/output/folder/
```
### 4. Generate count matrix
The final step generates a count matrix.
Arguments include:
- `--target`: path to output location(s)
- `--bam`: path to bam location(s)
- `--reference`: path to GTF annotation file
- `--workers`: threads for parallel processing
-- `save_csv/--save_mtx`: set up output format
--`group_novel`: whether group some novel isoforms that are potentially generated by read truncations together as one novel isoform.
--`platform`: 10x-ont or 10x-pacbio

Below is an example:
```
python3 src/main_preprocessing.py \
--task 'count matrix' \
--target path/to/output/folder/ \
--platform 10x-ont
--workers 8 --group_novel
```
With these steps, you should be able to generate counts that you can use for expression analysis or further study.

## Citation
Gouru, A., Xu, Z., Wang, K. “LIQA2: Efficient Isoform Detection and Quantification from Long Read Sequencing,” October 2024. International Conference on Intelligent Biology and Medicine, Houston, TX.
