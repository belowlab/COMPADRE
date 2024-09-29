FROM ubuntu:20.04

# Download and install dependencies
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
    perl \
    plink1.9 \
    plink2 \
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

## IN PROGRESS: install reference data from github releases URL 
RUN wget https://github.com/belowlab/compadre/releases/download/pre-release/compadre_data.tgz -O example_data.zip && unzip example_data.zip && rm example_data.zip

## after downloading, move things according accordingly 
RUN mv -r example_data/1KG /usr/src/lib/ && mv -r example_data/hapmap3 /usr/src/lib/ mv -r example_data/KDE_data /usr/src/lib/ 
RUN mkdir output

# Install the KernSmooth R package
RUN Rscript -e "install.packages('KernSmooth', repos='http://cran.rstudio.com/')"

# Add perl path to where primus is expecting it 
CMD ["mkdir", '-p', "/usr/src/perl"]
ENV PERL_PATH=/usr/src/perl
ENV PERL5LIB=$PERL_PATH:$PERL_PATH/lib/perl5:/usr/src/lib/perl_modules/:/usr/src/lib/perl_modules/PRIMUS/:$PERL5LIB

# Set primus executable entrypoint
WORKDIR /usr/src/bin
ENTRYPOINT ["perl", "./run_COMPADRE.pl"]