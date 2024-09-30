# COMPADRE

COMPADRE integrates genome-wide IBD sharing estimates from [PRIMUS](https://primus.gs.washington.edu/primusweb/index.html) 
and shared segment length- and quantity-based data from [ERSA](https://hufflab.org/software/ersa) to improve 
relationship estimation accuracy in family networks ahead of pedigree generation. 
We aim to extend the number and variety of constructed pedigrees derived from populations with increased data heterogeneity.

You can read more about these individual methods in the original [PRIMUS](https://compadre.dev/publications/primus.pdf) 
and [ERSA manuscripts](https://compadre.dev/publications/ersa.pdf).



## Updates

1. Added support for using 1000 Genomes Project genetic reference data to generate pairwise IBD estimates.
2. Added support for shared segments-based relationship estimation. This requires a file with pairwise shared segment data provided as input via the `segment_data` flag. We used [GERMLINE2](https://github.com/gusevlab/germline2) in our benchmarking, but there are a large number of tools that can generate these data. File formatting examples for this input can be found in the `example_data` folder(s).

    ```bash
    --segment_data example_data/simulations/AMR/amr_size20_segments.txt
    ```

    Note: COMPADRE does not require segment-specific IBD[1/2] status as part of the `--segment_data` input; however, inclusion of this information can improve the composite algorithm's performance. We have provided a generic script to identify IBD2 segments in standard IBD detection output: `tools/determine_ibd.py`. COMPADRE will check for the presence of an `ibd` column with values 1 or 2 at the last index of the `--segment_data` input file. 


## Installation

***Subject to change after adding Github Releases support ~ might not require regular cloning steps

Git: Click the green `Code` button at the top of this page and select a download option

Direct download: https://github.com/belowlab/compadre/archive/refs/heads/main.zip



## Execution

We have provided a Dockerfile to assist with installing dependencies and reference data. You must install and start the Docker client on your machine prior to building and running COMPADRE. Instructions to install Docker Engine on your system can be found [here](https://docs.docker.com/engine/install/).

Navigate into the compadre directory:

```bash
cd compadre
```

Copy your data inputs into the `input` folder (which will then be present in the Docker image):

```
cp -r /example/local/data/folder input/
```


Build the Docker image:

```bash
docker build -t compadre .
```


Run (interactive mode):

```bash
# Set entrypoint with port mapping from 
docker run -v /local/path/to/compadre_repo/output:/usr/src/output -p 6000:6000 -it --entrypoint /bin/bash compadre:latest 

# Run COMPADRE
perl run_COMPADRE.pl --file ../example_data/simulations/AMR/size20_0missing/amr_20_0 --segment_data ../example_data/simulations/AMR/amr_size20_segments.txt --genome --output ../output/amr_test --verbose 3
```


Run (non-interactive mode):

```bash
docker run -v /local/path/to/compadre_repo/output:/usr/src/output -p 6000:6000 compadre --file ../example_data/simulations/AMR/size20_0missing/amr_20_0 --segment_data ../example_data/simulations/AMR/amr_size20_segments.txt --genome --output ../output/amr_test --verbose 3
```


### Important execution notes
---
- In order to easily access COMPADRE results on your local machine, make sure to update the `-v` flag in the entrypoint step to reflect your own machine's path to the COMPADRE repository folder (specifically, the `output` folder). For example, on MacOS, this might be `/Users/yourname/compadre/output`. 
- Additional computation now takes place over an open socket, set to use port 6000 as default. If you need to use a different port, please indicate as such with the `--port_number=<INT>` flag (COMPADRE option) as well as the `-p` Docker option at runtime. 
- Please use standard PRIMUS runtime flags as detailed in the original [PRIMUS documentation](https://primus.gs.washington.edu/primusweb/res/documentation.html). 



## Additional Resources

The source code for generating family genetic data simulations can be found [here](https://github.com/belowlab/unified-simulations). 

More documentation:
- [Original PRIMUS docs](https://primus.gs.washington.edu/primusweb/res/documentation.html)
- [Original ERSA docs](https://hufflab.org/software/ersa/)

Please visit the [official COMPADRE website](https://compadre.dev/about) for publication updates and other details. 



## License

COMPADRE was developed by the [Below Lab](https://thebelowlab.com) at Vanderbilt University School of Medicine, and distributed under the following APACHE 2.0 license: https://compadre.dev/licenses/compadre_license.txt



## Questions?

Please email <strong><i>contact AT compadre DOT dev</strong></i> with the subject line "COMPADRE Help" or [submit an issue report](https://github.com/belowlab/compadre/issues). 