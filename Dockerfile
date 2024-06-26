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

# Download and install PLINK2
RUN wget https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64.zip && \
    unzip plink2_linux_x86_64.zip && \
    mv plink2 /usr/local/bin/ && \
    rm plink2_linux_x86_64.zip

# this might not be necessary anymore -- symbolic link to 'old' plink that primus might be expecting
RUN ln -s /bin/plink1.9 /bin/plink

# Download reference data
RUN wget https://compadre.dev/api/data/primus_reference_data.tgz -O lib/reference_data/primus_reference_data.tgz && \
    tar -xzvf lib/reference_data/primus_reference_data.tgz -C lib/reference_data && \
    rm lib/reference_data/primus_reference_data.tgz
    ## ^ This assumes that the unzipped folders in the .tgz download are correctly named, need to double check this 

# Install python packages 
RUN pip3 install --no-cache-dir -r requirements.txt

# Install the KernSmooth R package
RUN Rscript -e "install.packages('KernSmooth', repos='http://cran.rstudio.com/')"

# Add perl path to where primus is expecting it 
CMD ["mkdir", '-p', "/usr/src/primus-ersa-v2/perl"]
ENV PERL_PATH=/usr/src/perl
ENV PERL5LIB=$PERL_PATH:$PERL_PATH/lib/perl5:/usr/src/lib/perl_modules/:/usr/src/lib/perl_modules/PRIMUS/:$PERL5LIB

# Set primus executable entrypoint
WORKDIR /usr/src/bin
ENTRYPOINT ["perl", "./run_PRIMUS.pl"]
