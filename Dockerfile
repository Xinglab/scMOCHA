FROM chunjiesamliu/jrocker:latest

RUN apt-get update && apt-get install -y \
  wget \
  curl \
  samtools \
  sqlite3 \
  libdbd-sqlite3-perl \
  && rm -rf /var/lib/apt/lists/* \
  && pip install sinto \
  && pip install mgatk \
  && pip install pysam \
  && pip install matplotlib \
  && pip install numpy==1.23.5 \
  && mkdir -p /scMOCHA \
  && mkdir -p /opt/ \
  && cd /opt/ \
  && wget -O cellranger-7.1.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.1.0.tar.gz?Expires=1691655506&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci03LjEuMC50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2OTE2NTU1MDZ9fX1dfQ__&Signature=PVRqashGFzHsyqomOJu6KYl3cwj195ChzQHFUu7f5HkscDjlJOcHX8dSuPiYb-n5ZG81ASEkBpuJ8GY6hETioKmpQkze~8d1lQKWGAhUPDENNm8bq6ClqNxgEI5CHYSEeKQv9ZJN4lw2O~tXKYqRdjc4wEfzn4l8~dUZTCRBpbf~KJaRse42Cml7Su6QkyV2POazOSpRZQvpsbK5Hc~U3igrikqXTudpMdX8p99YRezGuFiDvkgnwHsRb0NOgQBKPAnHo1P0Pl~zMtg4g-F8MQ~HmD92fq3AFntp9ehdfXaciGLs2tHycmvRN-bKCW3KXkuIovNvSfmIfUTcE-IDuQ__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA" \
  && tar -xzvf cellranger-7.1.0.tar.gz \
  && ln -s /opt/cellranger-7.1.0/bin/cellranger /usr/bin/cellranger \
  && rm cellranger-7.1.0.tar.gz \
  && cd /opt/ \
  && wget -O bamtools-2.5.2.tar.gz https://github.com/pezmaster31/bamtools/archive/refs/tags/v2.5.2.tar.gz \
  && tar -xzvf bamtools-2.5.2.tar.gz \
  && cd bamtools-2.5.2/ \
  && mkdir build \
  && cd build \
  && cmake .. \
  && make \
  && make install \
  && cd /opt/ \
  && rm bamtools-2.5.2.tar.gz \
  && mkdir haplogrep3 \
  && cd haplogrep3 \
  && wget -O haplogrep3-3.2.1-linux.zip https://github.com/genepi/haplogrep3/releases/download/v3.2.1/haplogrep3-3.2.1-linux.zip \
  && unzip haplogrep3-3.2.1-linux.zip \
  && rm haplogrep3-3.2.1-linux.zip \
  && ln -s /usr/bin/python3 /usr/bin/python

COPY . /opt/scMOCHA
RUN R -f /opt/scMOCHA/packages.R
ENV PATH /opt/scMOCHA/bin:$PATH

WORKDIR /scMOCHA


# wget
# cellranger

# apt-get
# samtools

# pip
# sinto
# mgatk