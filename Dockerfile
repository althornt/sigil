FROM nfcore/base

WORKDIR /mnt/bin
COPY ./environment.yml .
RUN conda env create -f environment.yml && conda clean -a
ENV PATH /opt/conda/envs/sigl_prep/bin:$PATH
ENV PATH="$PATH:/mnt/bin"

# Install git and pip
RUN apt-get update && apt-get install -y \
  git \
  python3-pip \
  python3-dev \
  zlib1g-dev \
  rsync

# Install directory of stand-alone UCSC scripts ; can reduce the few scripts I use
RUN cd /mnt/bin \
        && rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/ ./

# Clone and install current MESA git repo
RUN cd /mnt/bin \
           && git clone https://github.com/BrooksLabUCSC/mesa.git \
           && cd mesa/ \
           && pip install . 