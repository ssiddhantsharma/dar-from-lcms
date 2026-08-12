# Ubuntu 24.04 satisfies both prebuilt UniDec binaries: native GCC-14 libstdc++ for isogen.so,
# and libhdf5-103-1t64 (HDF5 1.10, libhdf5_serial.so.103) for unideclinux. Do NOT use libhdf5-dev
# here -- it pulls HDF5 1.14 (.so.310), which unideclinux can't load.
FROM ubuntu:24.04
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        python3 python3-pip git ca-certificates \
        libhdf5-103-1t64 libhdf5-hl-100t64 libfftw3-dev libgomp1 \
    && rm -rf /var/lib/apt/lists/*
# requirements before COPY so editing dar_auto.py never re-triggers the (slow) install
COPY requirements.txt .
RUN pip3 install --break-system-packages --no-cache-dir -r requirements.txt
# Arial-metric font for publication-style figures
RUN apt-get update && apt-get install -y --no-install-recommends fonts-liberation \
    && rm -rf /var/lib/apt/lists/*
COPY dar_auto.py /opt/dar_auto.py
WORKDIR /work
