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

# Install python packages 
RUN pip3 install --no-cache-dir -r requirements.txt

# Install the KernSmooth R package
RUN Rscript -e "install.packages('KernSmooth', repos='http://cran.rstudio.com/')"

# Add perl path to where primus is expecting it 
WORKDIR /usr/src/primus-ersa-v2
RUN ln -s /bin/plink1.9 /bin/plink
CMD ["mkdir", '-p', "/usr/src/primus-ersa-v2/perl"]
ENV PERL_PATH=/usr/src/primus-ersa-v2/perl
ENV PERL5LIB=$PERL_PATH:$PERL_PATH/lib/perl5:/usr/src/primus-ersa-v2/lib/perl_modules/:/usr/src/primus-ersa-v2/lib/perl_modules/PRIMUS/:$PERL5LIB

# Set primus executable entrypoint
WORKDIR /usr/src/primus-ersa-v2/bin
ENTRYPOINT ["perl", "./run_PRIMUS.pl"]

# example run: 
# -----------------------------------------
# docker build -t compadre .
# docker run -it --entrypoint /bin/bash compadre:latest 

# run not in interactive mode:
# docker run compadre --file ../example_data/MEX_pop --genome -o ./testoutput --keep_inter_files -v 3

# --------------------------------

# perl run_PRIMUS.pl --file ../../final --genome --keep_inter_files --int_likelihood_cutoff 1 -v 3 -o ../output/test825
