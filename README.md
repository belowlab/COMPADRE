# COMPADRE (Alpha)

COMPADRE uses both genome-wide IBD estimates from the program [PRIMUS](https://primus.gs.washington.edu/primusweb/index.html) 
and shared segment length- and quantity-based data from the program [ERSA](https://hufflab.org/software/ersa) to build more accurate 
relationship estimates for each sample pair in a given non-directional graph of possible family members. 
We aim to extend the number and variety of constructed pedigrees derived from populations with increased admixture and sample missingness.

You can read more about these individual methods in the original [PRIMUS](https://compadre.dev/publications/primus.pdf) 
and [ERSA manuscripts](https://compadre.dev/publications/ersa.pdf).



## Changelog

1. Added support for 1000 Genomes Project reference data.
2. Added support for ERSA-derived shared segments relationship estimation ahead of pedigree reconstruction. This mode requires a [GERMLINE2](https://github.com/gusevlab/germline2)-generated .match file as input. Use the flag `--segment_data=/your/germline/path/here.match` 


## Installation

Via git:

```bash
# HTTPS
git clone https://github.com/belowlab/primus-ersa-v2.git

# SSH
git clone git@github.com:belowlab/primus-ersa-v2.git
```
Via direct download:

https://github.com/belowlab/primus-ersa-v2/archive/refs/heads/main.zip



## Execution

We have provided a Dockerfile to assist with installing dependencies and reference data. Instructions to install Docker Engine on your system can be found [here](https://docs.docker.com/engine/install/).

```bash
cd primus-ersa-v2

# Build
docker build -t compadre .

# Run (interactive)
docker run -it --entrypoint /bin/bash compadre:latest 

# Run (non-interactive):
docker run compadre --file ../example_data/MEX_pop --genome -o ./testoutput -v 3
```


## Additional Resources

The source code for generating family genetic data simulations can be found [here](https://github.com/belowlab/unified-simulations). 

More documentation:
- [Original PRIMUS docs](https://primus.gs.washington.edu/primusweb/res/documentation.html)
- [Original ERSA docs](https://hufflab.org/software/ersa/)

Please visit the [official COMPADRE website](https://compadre.dev/about) for publication information and other details. 



## License

COMPADRE was developed by the [Below Lab](https://thebelowlab.com) at Vanderbilt University School of Medicine, and distributed under the following APACHE 2.0 license: https://compadre.dev/licenses/compadre_license.txt



## Questions?

Please email <strong><i>contact AT compadre DOT dev</strong></i> with the subject line "COMPADRE Help" or submit an issue in this repository. 
