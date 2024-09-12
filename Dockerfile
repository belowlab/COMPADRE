FROM ubuntu:20.04

# Download and install Perl, PLINK1.9, R, Python3
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
    perl \
    plink1.9 \
    r-cran-devtools \
    wget \
    unzip \
    python3 \
    python3-pip \
    r-base

WORKDIR /usr/src

COPY . .

# this might not be necessary anymore -- symbolic link to 'old' plink that primus might be expecting
RUN ln -s /bin/plink1.9 /bin/plink

# Install python packages 
RUN pip3 install --no-cache-dir -r requirements.txt

#######
## IN PROGRESS: install reference data from github releases URL 
RUN wget https://github.com/belowlab/compadre/releases/download/v1.0.0/reference-data.dat

## after downloading, move it into the lib folder accordingly 
#######

# Install the KernSmooth R package
RUN Rscript -e "install.packages('KernSmooth', repos='http://cran.rstudio.com/')"

# Add perl path to where primus is expecting it 
CMD ["mkdir", '-p', "/usr/src/perl"]
ENV PERL_PATH=/usr/src/perl
ENV PERL5LIB=$PERL_PATH:$PERL_PATH/lib/perl5:/usr/src/lib/perl_modules/:/usr/src/lib/perl_modules/PRIMUS/:$PERL5LIB

# Set primus executable entrypoint
WORKDIR /usr/src/bin
ENTRYPOINT ["perl", "./run_PRIMUS.pl"]