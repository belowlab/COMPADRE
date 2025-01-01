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

# Install Ghostscript -- TBD

WORKDIR /usr/src

COPY . .

# symbolic link to 'old' plink that primus is expecting
RUN ln -s /bin/plink1.9 /bin/plink

# Install python packages 
RUN pip3 install -r requirements.txt

## IN PROGRESS: install reference data from github releases URL 
RUN wget https://github.com/belowlab/compadre/releases/download/pre-release/compadre_data.zip && unzip compadre_data.zip && rm compadre_data.zip

## after downloading, move things according accordingly 
RUN mv compadre_data/1KG /usr/src/lib && mv compadre_data/hapmap3 /usr/src/lib && mv compadre_data/KDE_data /usr/src/lib && mv compadre_data/example_data .

# clean up unzipped empty folders
RUN rm -r __MACOSX && rm -r compadre_data

# Install the KernSmooth R package
RUN Rscript -e "install.packages('KernSmooth', repos='http://cran.rstudio.com/')"

# Add perl path to where primus is expecting it 
CMD ["mkdir", '-p', "/usr/src/perl"]
ENV PERL_PATH=/usr/src/perl
ENV PERL5LIB=$PERL_PATH:$PERL_PATH/lib/perl5:/usr/src/lib/perl_modules/:/usr/src/lib/perl_modules/PRIMUS/:$PERL5LIB

# Change this if you need to use a different port
EXPOSE 6000
ENV COMPADRE_HOST=0.0.0.0

# Set primus executable entrypoint
WORKDIR /usr/src/bin
ENTRYPOINT ["perl", "./run_COMPADRE.pl"]