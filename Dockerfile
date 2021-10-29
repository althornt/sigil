FROM nfcore/base

WORKDIR /mnt/bin
COPY ./environment.yml .
RUN conda env create -f environment.yml && conda clean -a
ENV PATH /opt/conda/envs/sigl_prep/bin:$PATH

#Install git
RUN apt-get update && apt-get install -y \
  git \
  python3-pip \
  python3-dev \
  zlib1g-dev

#clone and install mesa git repo
RUN cd /mnt/bin \
           && git clone https://github.com/BrooksLabUCSC/mesa.git \
           && cd mesa/ \
           && pip install . \
