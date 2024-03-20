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
  && wget -O cellranger-8.0.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-8.0.0.tar.gz?Expires=1711015460&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=NNDHxw14LOdP14k5NhlkMLvE9tWutdtw0C3NgOeBQ4O-s8I3a0FO0L45eXYcnHm8N0IKDS4YIH49ryxUVNhveaJeFZIsjMirMZ92SaA72i7y6Xfr~i9EKE42jfTUVk8~wWby9BaWiIj75ZB9ERfhKYkYXirYT0mkef114-8x1hGBZWtFwks0wDk6TLjGiOHDlLOJIlb3lRGqysVzxGAtE9HSpFEpN4Na8CbKT6qWEt4SCEDoAv6SL8xqP0Z4ZaS4kUaMoS2Qv5SnZrVFqojRYs61P3wcYSiHm8sKAXABQp8EKoFj8-QHHe-zUe4qQwY3IYl5gUn2djQBB9ew7KQEng__" \
  && tar -xzvf cellranger-8.0.0.tar.gz \
  && ln -s /opt/cellranger-8.0.0/bin/cellranger /usr/bin/cellranger \
  && rm cellranger-8.0.0.tar.gz \
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