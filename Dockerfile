FROM radioastro/base
MAINTAINER bhugo@ska.ac.za

RUN apt-get update
RUN apt-get install -y software-properties-common
RUN add-apt-repository -y -s ppa:radio-astro/main
RUN apt-get update

RUN apt-get install -y python-casacore
RUN apt-get install -y python-pip
RUN pip install pip -U

ADD RFIMasker /src/RFIMasker/RFIMasker
ADD MANIFEST.in /src/RFIMasker/MANIFEST.in
ADD requirements.txt /src/RFIMasker/requirements.txt
ADD setup.py /src/RFIMasker/setup.py
ADD README.md /src/RFIMasker/README.md

RUN pip install /src/RFIMasker

ENTRYPOINT ["mask_ms.py"]
CMD ["--help"]
