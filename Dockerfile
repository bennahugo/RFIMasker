FROM kernsuite/base:2
MAINTAINER bhugo@ska.ac.za

RUN docker-apt-install python-casacore python-pip

RUN pip install pip -U

ADD RFIMasker /src/RFIMasker/RFIMasker
ADD MANIFEST.in /src/RFIMasker/MANIFEST.in
ADD requirements.txt /src/RFIMasker/requirements.txt
ADD setup.py /src/RFIMasker/setup.py
ADD README.md /src/RFIMasker/README.md

RUN pip install /src/RFIMasker

ENTRYPOINT ["mask_ms.py"]
CMD ["--help"]
