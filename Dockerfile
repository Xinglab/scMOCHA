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
  && pip install matplotlib \
  && mkdir -p /scMOCHA \
  && mkdir -p /opt/ \
  && cd /opt/ \
  && wget -O cellranger-7.1.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.1.0.tar.gz?Expires=1691585485&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci03LjEuMC50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2OTE1ODU0ODV9fX1dfQ__&Signature=eLFGzxQwPtbq3F4pYhzVT-91FX2Xvsfc3-xXRb0lojMyFOEUaDaUXPZN-nKIwqtA4I2arFkd71aLn03EnldxGVLjAFL3WS47kaD6iIb0NmyoD5i8X5oJN5zYi9NR20qhsGAdqQitM2g8iw5R5uI4ZsCim9KoXmjEsO11wzvRqaj22tAgC26iMv5U-0wYtSwLIFG6qf8eeSOADsWEUJCpESxJd-DFNl7s1TqR6OHU10yYpH1X4XizChtTgV9xz8Diar1QjyLmakBKkmjmyrEu5zxuvv7teD2Sa0T4KSQw9~b2tPLPiO9m33ZoUG5gLYISmoZCsrQT9v-aOp7lS60~BA__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA" \
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