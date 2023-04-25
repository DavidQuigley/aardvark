# syntax=docker/dockerfile:1
FROM centos:7
LABEL version="1.0" maintainer="David Quigley <david.quigley@ucsf.edu>"
WORKDIR /opt/software
RUN yum install -y epel-release

# Enable the Extra Packages for Enterprise Linux (EPEL) repository
RUN yum -y install https://dl.fedoraproject.org/pub/epel/epel-release-latest-7.noarch.rpm
# On RHEL 7, enable the Optional repository
#RUN subscription-manager repos --enable "rhel-*-optional-rpms"
# Download and install R
RUN curl -O https://cdn.rstudio.com/r/centos-7/pkgs/R-4.1.3-1-1.x86_64.rpm
RUN yum -y install R-4.1.3-1-1.x86_64.rpm
RUN yum -y install libpng*
RUN yum -y install libxml2-devel
RUN yum -y install openssl-devel
RUN yum -y install wget
# Create symlink to R
RUN ln -s /opt/R/4.1.3/bin/R /usr/local/bin/R
RUN ln -s /opt/R/4.1.3/bin/Rscript /usr/local/bin/Rscript

RUN wget https://github.com/DavidQuigley/aardvark/releases/download/v-0.2.5-beta/aardvark_0.2.5.tar.gz

RUN Rscript -e 'install.packages("getopt", repos="https://cran.case.edu")'
RUN Rscript -e 'install.packages("optparse", repos="https://cran.case.edu")'
RUN Rscript -e 'install.packages("stringr", repos="https://cran.case.edu")'
RUN Rscript -e 'install.packages("BiocManager", repos="https://cran.case.edu")'
RUN Rscript -e 'install.packages("testthat", repos="https://cran.case.edu")'
RUN Rscript -e 'BiocManager::install("GenomicRanges", site_repository="https://cran.rstudio.com")'
RUN Rscript -e 'BiocManager::install("VariantAnnotation", site_repository="https://cran.rstudio.com")'
RUN Rscript -e 'BiocManager::install("GenomicAlignments", site_repository="https://cran.rstudio.com")'
RUN Rscript -e 'BiocManager::install("Rsamtools", site_repository="https://cran.rstudio.com")'
RUN Rscript -e 'BiocManager::install("Biostrings", site_repository="https://cran.rstudio.com")'
RUN Rscript -e 'BiocManager::install("BSgenome.Hsapiens.UCSC.hg38", site_repository="https://cran.rstudio.com")'
RUN Rscript -e 'BiocManager::install("BSgenome.Hsapiens.UCSC.hg19", site_repository="https://cran.rstudio.com")'
RUN Rscript -e 'BiocManager::install("biomaRt", site_repository="https://cran.rstudio.com")'
RUN Rscript -e 'BiocManager::install("ensembldb", site_repository="https://cran.rstudio.com")'
RUN Rscript -e 'BiocManager::install("EnsDb.Hsapiens.v86", site_repository="https://cran.rstudio.com")'

RUN yum -y install libjpeg*
RUN yum -y install centos-release-scl-rh
RUN yum -y install devtoolset-11
RUN scl enable devtoolset-11 bash
RUN Rscript -e 'BiocManager::install("Gviz", site_repository="https://cran.rstudio.com")'
RUN Rscript -e 'install.packages("/opt/software/aardvark_0.2.5.tar.gz", source=TRUE)'
