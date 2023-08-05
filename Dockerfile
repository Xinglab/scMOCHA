FROM chunjiesamliu/jrocker:latest

RUN apt-get update && apt-get install -y \
  wget \
  samtools \
  && rm -rf /var/lib/apt/lists/* \
  && pip install sinto \
  && pip install mgatk \
  && makdir -p /scMOCHA \
  && mkdir -p /opt/ \
  && cd /opt/ \
  && wget -O cellranger-7.1.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.1.0.tar.gz?Expires=1691310147&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci03LjEuMC50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2OTEzMTAxNDd9fX1dfQ__&Signature=kN3pjE12j02vmerw3oZPz8tZjbeE1YvYIZUbxNsaaCp3BtvB06ABzMVn0v7Kqi6QRnSlYvUeP8aN6oE3qHDp8jKMD85GjmxCgyojGI7Tj94w9DanYDCJz0943~tZGXFBdDnOagB2EIhrrv-AzNGSnpa1w8a37Gi0c8~QznVCB5JLucAPAFL3Y-PX7nSR52sXZUaZ5zqIwaszpYnPo07hJ8ab51eJpUPhnLvEKWIrgTt~dTOoCit2JNlB1rXILeykTc4SlIr6fPyDrhLzSPGkJUn~cqdISYMgaBUu9UOAl8kusC0DaBj7Ddk4FLqEZnl9dvfTqi7T4YOoWLnchCCEtA__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA" \
  && tar -xzvf cellranger-7.1.0.tar.gz \
  && ln -s /opt/cellranger-7.1.0/bin/cellranger /usr/bin/cellranger \
  && rm cellranger-7.1.0.tar.gz \
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