# COMPADRE

COMPADRE integrates genome-wide IBD sharing estimates from [PRIMUS](https://primus.gs.washington.edu/primusweb/index.html)
and shared segments length and distribution data from [ERSA](https://hufflab.org/software/ersa) to improve
relationship estimation accuracy in family networks ahead of pedigree generation. COMPADRE aims to extend the number and variety of constructed pedigrees derived from populations with increased data heterogeneity.

## Key Improvements

1. COMPADRE integrates shared segments-based relationship estimation ahead of pedigree reconstruction. We used [GERMLINE2](https://github.com/gusevlab/germline2) to identify shared segments for our benchmarking, but there are many other tools that can generate these data. Specifically, COMPADRE expects a file with the following columns: `id1 id2 start end length chrom`. COMPADRE can read a .gz or zstd compressed segment file. All standard ERSA runtime options are customizable during this step; simply use the applicable flag(s) at runtime. Descriptions of each ERSA flag can be found using the `--help` message as well as on the [COMPADRE website](https://compadre.dev/docs).

    Note: COMPADRE does not <i>require</i> segment-specific IBD status ('half'/1 or 'full'/2) as part of this input; however, inclusion of this information can improve the composite algorithm's performance. COMPADRE will check for the presence of an `ibd` column containing values 1 or 2 at the last index of the `--segment_data` input file. We have provided a generic script to identify IBD2 segments from standard IBD detection output here: `tools/determine_ibd.py`. It is also possible to run COMPADRE without shared segments input if desired.

    Another note: ERSA uses a 2.5 cM default threshold for analyzing IBD segments. Please use the `--min_cm` flag at runtime to specify a different threshold.

2. COMPADRE supports optional PADRE computation after completion of standard network reconstruction. Use the `--run_padre` flag at runtime.

3. COMPADRE utilizes 1000 Genomes Project genetic reference data to generate pairwise IBD estimates. This update also leverages a support vector machine (SVM) algorithm dynamically trained on 1000 Genomes Project subpopulation PCA results to predict ancestry ahead of IBD estimation and reconstruction.

## Installation

Git: Click the green `Code` button at the top of this page and select a cloning option.

## Execution

### Option 1: Docker

We have provided a Dockerfile to help install dependencies and reference data ahead of running COMPADRE. First, however, you must install and launch the Docker client on your machine. Instructions to install Docker Engine on your system can be found [here](https://docs.docker.com/engine/install/).

Navigate into the `compadre` directory downloaded from GitHub:

```bash
# Assuming you cloned the repository in your current directory path

cd compadre
```

First, copy your input data into the `input` folder (to be automatically copied into the Docker image):

```bash
# (cd into the COMPADRE directory first)

cp /example/local/folder/fileset.bed input/
cp /example/local/folder/fileset.bim input/
cp /example/local/folder/fileset.fam input/
cp /example/local/folder/segments.txt input/
```

Next, build the Docker image:

```bash
docker build -f Dockerfile.github -t compadre .
```

Finally, run (in interactive mode):

```bash
# Step 1: Set entrypoint to bring you into the Docker image location
docker run -v \
    /local/path/to/compadre_repo/output:/usr/src/output \
    -p 4000:4000 -it --entrypoint /bin/bash compadre:latest 

# Step 2: Run COMPADRE (replace inputs with your own from the /usr/src/input/ folder)
perl run_COMPADRE.pl \
    --file ../example_data/input \
    --segment_data ../example_data/segments.txt \
    --output ../output/compadre_test \
    --genome --verbose 1 --port_number 4000

## NOTE: This example uses a non-default port, hence the --port_number flag.
## ALL ports, however, require the docker '-p PORT:PORT' flag in the first command.

# Step 3: Once complete, press CTRL+A+D to exit the image.
```

... or run in non-interactive mode:

```bash
# This example also uses a non-default port
docker run -v /local/path/to/compadre_repo/output:/usr/src/output -p 4000:4000 compadre \
    --file ../example_data/input \
    --segment_data ../example_data/segments.txt \
    --output ../output/compadre_test \
    --genome --verbose 1 --port_number 4000
```

<u><strong>NOTE</strong></u>: The "Run" examples above perform all steps of COMPADRE: (1) input data quality control, (2) identification of an unrelated set, and (3) pedigree reconstruction. While performing the first two of these steps is encouraged in most instances, if you already have PLINK *.genome formatted data (and performed quality control), you can skip to pedigree reconstruction by using both `--no_IMUS` and `--plink_ibd <yourfile.genome>` flags. Conversely, if you only want to generate IBD estimates per network (step 1) and the overall unrelated set (step 2) from your standard input data, you can use the `--no_PR` flag to stop execution before pedigree reconstruction.

### Option 2: Singularity

COMPADRE also offers Singularity support. This option is recommended for use in HPC environments without Docker permissions.

```bash
# Build image file using Docker Hub link
singularity pull compadre.sif docker://grahamebelowlab/compadre:latest

# Run with bind mounts to connect the image to your local data and output folders
singularity run \
    --bind /local/path/for/your/input/data:/compadre_data \
    --bind /local/path/for/your/output/folder:/compadre_output \
    compadre.sif \
    --file /data/your_plink_fileset \
    --segment_data /compadre_data/your_segments.txt \
    --genome \
    --output /compadre_output/results
```

### Execution Notes

- In order to easily access COMPADRE results on your local machine, use the `-v` flag in the Docker entrypoint step to link your local COMPADRE repository folder path (specifically, the `output` folder). For example, on macOS, this might be `/Users/yourname/Downloads/compadre/output` if you cloned this repository into your Downloads folder.
- Additional computation now takes place over an open socket. COMPADRE defaults to port 6000; if you need to use a different port, please indicate as such with the following COMPADRE flag `--port_number <INT>` AND Docker flag `--publish <INT>:<INT>`. See the "Run" examples for more details.
- All other runtime options are detailed [here](https://compadre.dev/docs).

## Data Simulation Resources

The source code for generating family genetic data simulations can be found [here](https://github.com/belowlab/unified-simulations).

A simple script to simulate sample ages for a given pedigree can be found in the `tools/` folder in this repository. Briefly, the first step is assignment of generations to each sample by identifying the pedigree founders as generation 0, then recursively assigning their children to generation 1, grandchildren to generation 2, and so on. Next, the script generates ages by giving each generation a base age that's 25 years older than the generation below it (founders at ~100 years, their children at ~75 years, etc.), then iteratively adjusts these ages to maintain realistic parent-child age gaps of at least 14 years and similar ages between spouses.

More details on PRIMUS, ERSA, and PADRE can be found in their respective documentation:

- [PRIMUS](https://primus.gs.washington.edu/primusweb/res/documentation.html)
- [ERSA](https://hufflab.org/software/ersa) (under maintenance)
- [PADRE](https://hufflab.org/software/padre) (under maintenance)

## Running Integration Test (This section is for individuals developing COMPADRE or adding features)

Following standard PERL convention, we have added integration test to COMPADRE in the subdirectory, t/. Currently, these test only check to make sure that the code output the same files and results as the commands in the "Execution" section. Future progress will add more rigorous testing. These integration test are run by on a Github action automatically, but we ask that developers run these test locally to ensure their feature/changes are working as expected. Users can run the following command using "prove" (bundled with the PERL runtime) to run the test suite.

```bash
# We assume that you have already activated the appropriate Conda environment and that you are running this from the home directory of the COMPADRE repository after cloning it from github
prove -r t/ -v
```

*note*: These tests are not currently implemented in the docker images.

## Questions?

Please email <strong><i>contact AT compadre DOT dev</strong></i> with the subject line "COMPADRE Help" or [submit an issue report/pull request on GitHub](https://github.com/belowlab/compadre/issues).

If you use COMPADRE in your research, please cite the following:

```
Evans, G. F., Baker, J. T., Petty, L. E., ... & Below, J. E. (2025). COMPADRE: 
Combined pedigree-aware distant relatedness estimation for improved pedigree reconstruction. 
The American Journal of Human Genetics. DOI: 10.1016/j.ajhg.2025.09.011
```

## License

COMPADRE was developed by the [Below Lab](https://thebelowlab.com) in the Division of Genetic Medicine at Vanderbilt University Medical Center, Nashville, TN, USA.

COMPADRE is distributed under the following APACHE 2.0 license: <https://compadre.dev/licenses/compadre_license.txt>
