FROM gcr.io/diamond-privreg/xchem/ccp4:7.1 as ccp4
FROM gcr.io/diamond-privreg/xchem/phenix:1.20 as phenix

FROM buildpack-deps:bullseye-curl as devel

COPY --from=ccp4 /ccp4-7.1 /ccp4-7.1
COPY --from=phenix /phenix-1.20 /phenix-1.20

RUN apt-get update \
    && apt-get upgrade -y \
    && apt-get install -y \
        libxrender1 libfontconfig libxext6 \
        libglib2.0-0 libsm6 libxi6 libxrandr2 libxfixes3 libxcursor1 libxinerama1 \
        libgomp1 libxdamage1 libxcb-shm0 libxcb-render0 \
        python3-pip \
    && apt-get clean -y && rm -rf /var/lib/apt/lists/* \
    && echo . /ccp4-7.1/bin/ccp4.setup-sh | tee ~/.bashrc ~/.zshrc \
    && curl https://bootstrap.pypa.io/pip/2.7/get-pip.py -o /tmp/get-pip.py \
    && /ccp4-7.1/bin/ccp4-python /tmp/get-pip.py \
    && /ccp4-7.1/bin/ccp4-python -m pip install flake8 \
    && pip3 install black[python2]

ENV QT_X11_NO_MITSHM=1 \
    PATH=$PATH:/ccp4-7.1/bin:/phenix-1.20/build/bin