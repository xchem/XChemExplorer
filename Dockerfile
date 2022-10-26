FROM gcr.io/diamond-privreg/xchem/ccp4:7.1 as ccp4

FROM registry.hub.docker.com/library/debian:bullseye as xce 

COPY --from=ccp4 /ccp4-7.1 /ccp4-7.1

ARG XCE_DIR=/xce
WORKDIR ${XCE_DIR}

COPY . ${XCE_DIR}

RUN apt-get update \
    && apt-get upgrade -y \
    && apt-get install -y \
        libxrender1 libfontconfig libxext6 \
        libglib2.0-0 libsm6 libxi6 libxrandr2 libxfixes3 libxcursor1 libxinerama1 \
        libgomp1 libxdamage1 libxcb-shm0 libxcb-render0 \
    && apt-get clean -y && rm -rf /var/lib/apt/lists/*

ENV QT_X11_NO_MITSHM=1
    XChemExplorer_DIR=${XCE_DIR}

ENTRYPOINT /ccp4-7.1/bin/ccp4-python /xce/XChemExplorer.py
