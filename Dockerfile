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
  && wget -O cellranger-7.1.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.1.0.tar.gz?Expires=1698133810&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=K03-tm~GbfexOYcuPPbpIN8m0jGtCE51Jajc0TVjp1S3Yn9w0raTWtQkLyWtEvOG8nNxK1km5o-I5SQPK9VzVMHXCyE5BthwGSDiLEd7Iwju2kx56u203Wy8BX2tho7IZxl7UWy8-GF1CvtpLnDsYJaIgozSMxH1cmEAAIoV0C5f0kF2BHfmNQGos2Ytvk2DOJMv5KvbrbshPjlgRNu81qa9kpO9aQJvZN-UFB58TD0Hod3OGvLT1FxxSowOLjFhsJg0bShCWvJ7GkaArIdxpHSVH8ET8CwljrbMLMo9miDgFGQSP-e40awxELgNz5RfCru1dRxD8xWR5q0NFQYFuA__" \
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