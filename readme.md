# scl-ppin-reweighting

### Todo

- add sample usage
- add format project with ruff

## Abstract

Protein complexes are important features of a protein-protein interaction network (PPIN). Protein complexes aid in understanding organisms at the cellular level [7] and in mapping diseases [5]. Previous researchers have done protein complex prediction through variations of Markov Clustering (MCL) [2, 6], P5COMP [4] and MCODE [1]. The algorithms tackled propose improvements to the post-processing step of the prediction pipeline such as MCL-CAw [6], a refined version of MCL. Functional Similarity weighting (FS-weighting) [3] is an algorithm used in the pre-processing of the complex prediction process [2]. A paper on predicting protein functions using protein-protein interactions proposed utilizing indirect neighbors to improve the performance of FS-weighting [3]. This paper proposes using FS-weighting including higher level indirect neighbors to improve the precision and recall of existing protein complex prediction pipelines.

## Requirements

**Important: All commands that will be provided in this documentation must be run from the repository root. All commands assume that the `bash` shell is used.**

This software was written in Python 3.10, and requires some external dependencies, provided in the `requirements.txt` file.

To set up the environment for running the software:
```
# Set up virtual env
python3 -m venv .venv

# Activating virtual env (must be done every new instance of the shell)
source .venv

# installing the necessary packages
pip install -r requirements.txt
```

## Usage

### Reweighting

The main entry point for the reweighting algorithm is in `main.py`. Running `python3 main.py -h` will display the following help message:

```
$ python3 main.py -h

usage: SCL 2024-2025 Thesis [-h] --input_graph INPUT_GRAPH --database_path DATABASE_PATH [--depth_from DEPTH_FROM]
                            [--depth_to DEPTH_TO] [--verbose] [--path PATH] [--cache]

PPIN Reweighter using Topology and Experimenal Source

options:
  -h, --help            show this help message and exit
  --input_graph, -I INPUT_GRAPH
                        path of input graph
  --database_path, -D DATABASE_PATH
                        path to BioGRID database
  --depth_from, -df DEPTH_FROM
                        include starting from to depth n (must be geq than 0)
  --depth_to, -dt DEPTH_TO
                        include up to depth n (must be geq than 0)
  --verbose, -v
  --path, -p PATH       path to save reweighted graph
  --cache, -c           cache computations?

```


### Visualizing

The clusters outputted from the `P5COMP` pipeline can be visualized with the provided `visualize_clusters.py` file. Running `python3 visualize_clusters.py -h` will display the following help message:

```
$ python3 visualize_clusters.py

Visualize clusters produced by p5comp algo

positional arguments:
  file        path for the cluster file

options:
  -h, --help  show this help message and exit
```

It is to be noted that the program will produce an `index.html` file and a `lib` directory within the root directory. These files are the ones that will be used to display the visualizations of the clusters. You can use these files to host the visualizations webpage yourself (i.e. running `python3 -m http.server` within the root directory to be able to view the html page in the browser).


# References

[1] Gary D Bader and Christopher WV Hogue. An automated method for finding molecular complexes
in large protein interaction networks. BMC Bioinformatics, 4:2, 2003.

[2] Jerome Beltran, Catalina Montes, John Justine Villar, and Adrian Roy Valdez. A hybrid method
for protein complex prediction in weighted protein-protein interaction networks. Philippine Science
Letters, 10:50–57, 2017.

[3] Hon Nian Chua, Wing-Kin Sung, and Limsoon Wong. Exploiting indirect neighbours and topological
weight to predict protein function from protein–protein interactions. Bioinformatics, 22(13):1623–
1630, 04 2006.

[4] Yenzy Urson S. Hebron, Francis Donald P. Alfonso, and John Justine S. Villar. P5COMP: Parameter-
free Pipeline for Predicting Problematic Protein Complexes. PhD thesis, University of the Philippines
Diliman, 2024.

[5] Md. Shahidul Islam, Md. Rafiqul Islam, and A.B.M. Shawkat Ali. Protein complex prediction in large
protein–protein interaction network. Informatics in Medicine Unlocked, 30:100947, 2022.

[6] Sriganesh Srihari, Kang Ning, and Hon Wai Leong. Mcl-caw: a refinement of mcl for detecting yeast
complexes from weighted ppi networks by incorporating core-attachment structure. BMC Bioinfor-
matics, 11, 10 2010.

[7] Javad Zahiri, Abbasali Emamjomeh, Samaneh Bagheri, Asma Ivazeh, Ghasem Mahdevar, Hessam
Sepasi Tehrani, Mehdi Mirzaie, Barat Ali Fakheri, and Morteza Mohammad-Noori. Protein complex
prediction: A survey. Genomics, 112(1):174–183, 2020.