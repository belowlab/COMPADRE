# COMPADRE

COMPADRE integrates genome-wide IBD sharing estimates from [PRIMUS](https://primus.gs.washington.edu/primusweb/index.html) 
and shared segments length and distribution data from [ERSA](https://hufflab.org/software/ersa) to improve 
relationship estimation accuracy in family networks ahead of pedigree generation. COMPADRE aims to extend the number and variety of constructed pedigrees derived from populations with increased data heterogeneity.



## Updates

1. COMPADRE integrates shared segments-based relationship estimation ahead of pedigree reconstruction. We used [GERMLINE2](https://github.com/gusevlab/germline2) to identify shared segments for our benchmarking, but there are many other tools that can generate these data. Specifically, COMPADRE expects a file with the following columns: `id1 id2 start end length chrom`. A file formatting example for this input can be found in the `example_data` folder.

    ```bash
    --segment_data example_data/simulations/EUR/eur_size20_segments.txt
    ```

    Note: COMPADRE does not <i>require</i> segment-specific IBD status ('half'/1 or 'full'/2) as part of this input; however, inclusion of this information can improve the composite algorithm's performance. COMPADRE will check for the presence of an `ibd` column containing values 1 or 2 at the last index of the `--segment_data` input file. We have provided a generic script to identify IBD2 segments from standard IBD detection output here: `tools/determine_ibd.py`. It is also possible to run COMPADRE without shared segments input if desired.

2. COMPADRE supports optional PADRE computation after completion of standard network reconstruction. Use the `--run_padre` flag at runtime. 

3. COMPADRE utilizes 1000 Genomes Project genetic reference data to generate pairwise IBD estimates. This update also leverages a support vector machine (SVM) algorithm trained on pre-processed PLINK PCA results to predict ancestry ahead of IBD estimation and reconstruction.



## Installation

Git: Click the green `Code` button at the top of this page and select a download option

Direct download: https://github.com/belowlab/compadre/archive/refs/heads/main.zip



## Execution

We have provided a Dockerfile to help install dependencies and reference data ahead of running COMPADRE. First, however, you must install and launch the Docker client on your machine. Instructions to install Docker Engine on your system can be found [here](https://docs.docker.com/engine/install/).

Navigate into the `compadre` directory downloaded from GitHub:

```bash
# Assuming you cloned the repository in your current directory path

cd compadre
```

Copy your input data into the `input` folder (to be automatically copied into the Docker image):

```bash
# cd into the COMPADRE directory first

cp /example/local/folder/test.bed ./input/
cp /example/local/folder/test.bim ./input/
cp /example/local/folder/test.fam ./input/
cp /example/local/folder/segments.txt ./input/
```

Build the Docker image:

```bash
docker build -t compadre .
```

Run (interactive mode):

```bash
# Set entrypoint to bring you into the Docker image location

docker run -v \
    /local/path/to/compadre_repo/output:/usr/src/output \
    -p 4000:4000 -it --entrypoint /bin/bash compadre:latest 

# Run COMPADRE (replace inputs with your own from the /usr/src/input/ folder)

perl run_COMPADRE.pl \
    --file ../example_data/simulations/EUR/size20_0missing/eur_20_0 \
    --segment_data ../example_data/simulations/EUR/eur_size20_segments.txt \
    --output ../output/eur_test \
    --genome --verbose 1 --run_padre --port_number 4000

## NOTE: This example uses a non-default port, hence the --port_number flag.
## ALL ports, however, require the docker '-p PORT:PORT' flag in the first command.
```

Run (non-interactive mode):

```bash
docker run -v /local/path/to/compadre_repo/output:/usr/src/output -p 4000:4000 compadre \
    --file ../example_data/simulations/AMR/size20_0missing/amr_20_0 \
    --segment_data ../example_data/simulations/AMR/amr_size20_segments.txt \
    --output ../output/amr_test \
    --genome --verbose 1 --run_padre --port_number 4000
```

<u><strong>NOTE</strong></u>: The "Run" examples above perform all steps of COMPADRE: (1) input data quality control, (2) identification of an unrelated set, and (3) pedigree reconstruction. While performing the first two of these steps is encouraged in most instances, if you already have PLINK *.genome formatted data (and performed quality control), you can skip to pedigree reconstruction by using both `--no_IMUS` and `--plink_ibd <yourfile.genome>` flags. Conversely, if you only want to generate IBD estimates per network and the overall unrelated set from your standard input data, you can use the `--no_PR` flag.



### Execution notes
---
- In order to easily access COMPADRE results on your local machine, use the `-v` flag in the Docker entrypoint step to link your local COMPADRE repository folder path (specifically, the `output` folder). For example, on macOS, this might be `/Users/yourname/Downloads/compadre/output` if you cloned this repository into your Downloads folder. 
- Additional computation now takes place over an open socket. COMPADRE defaults to port 6000; if you need to use a different port, please indicate as such with the following COMPADRE flag `--port_number <INT>` AND Docker flag `--publish <INT>:<INT>`. See the "Run" examples for more details.
- All other runtime options are detailed [here](https://compadre.dev/docs). 



## Additional Resources

The source code for generating family genetic data simulations can be found [here](https://github.com/belowlab/unified-simulations). 

More details on PRIMUS, ERSA, and PADRE can be found in their respective documentation:
- [PRIMUS](https://primus.gs.washington.edu/primusweb/res/documentation.html)
- [ERSA](https://hufflab.org/software/ersa/)
- [PADRE](https://hufflab.org/software/padre/)

Please visit the [official COMPADRE website](https://compadre.dev/about) for publication updates and other information. 



## Questions?

Please email <strong><i>contact AT compadre DOT dev</strong></i> with the subject line "COMPADRE Help" or [submit an issue report/pull request on GitHub](https://github.com/belowlab/compadre/issues). 

If you use COMPADRE in your research, please cite the following:
```
Evans GF, Baker JT, Petty LE, Petty AS, Polikowsky HG, Bohlender RJ, Chen HH, Chou CY, 
Viljoen KZ, Beilby JM, Kraft SJ, Zhu W, Landman JM, Morrow AR, Bian D, Scartozzi AC, 
Huff CD, Below JE. Combined Pedigree-aware Distant Relatedness Estimation (COMPADRE): 
improved pedigree reconstruction using integrated relationship estimation approaches. 
[Publication details forthcoming]
```



## License

COMPADRE was developed by the [Below Lab](https://thebelowlab.com) in the Division of Genetic Medicine at Vanderbilt University Medical Center, Nashville, TN, USA. 

COMPADRE is distributed under the following APACHE 2.0 license: https://compadre.dev/licenses/compadre_license.txt