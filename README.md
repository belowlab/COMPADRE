# COMPADRE (Alpha)

COMPADRE uses both genome-wide IBD estimates from the program [PRIMUS](https://primus.gs.washington.edu/primusweb/index.html) 
and shared segment length- and quantity-based data from the program [ERSA](https://hufflab.org/software/ersa) to build more accurate 
relationship estimates for each sample pair in a given non-directional graph of possible family members. 
We aim to extend the number and variety of constructed pedigrees derived from populations with increased admixture and sample missingness.

You can read more about these individual methods in the original [PRIMUS](https://compadre.dev/publications/primus.pdf) 
and [ERSA manuscripts](https://compadre.dev/publications/ersa.pdf).


## Installation

Via git:

```bash
# HTTPS
git clone https://github.com/belowlab/primus-ersa-v2.git

# SSH
git clone git@github.com:belowlab/primus-ersa-v2.git

cd primus-ersa-v2
```
Via zip file:

https://github.com/belowlab/primus-ersa-v2/archive/refs/heads/main.zip



## Execution

We have provided a Dockerfile to assist with installing dependencies and reference data. You must first install Docker Engine on your system [here](https://docs.docker.com/engine/install/).

```bash
# Build
docker build -t compadre .

# Run (interactive)
docker run -it --entrypoint /bin/bash compadre:latest 

# Run (non-interactive):
docker run compadre --file ../example_data/MEX_pop --genome -o ./testoutput -v 3
```


## Questions?

Please email grahame.f.evans AT vanderbilt DOT edu, or submit an issue in this repository. 



## License

COMPADRE was developed by the [Below Lab](https://thebelowlab.com) at Vanderbilt University School of Medicine, and licensed under the following APACHE 2.0 license: https://compadre.dev/software/license.txt