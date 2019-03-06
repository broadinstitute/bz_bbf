FROM python:3.7

RUN pip3 install -U setuptools
RUN pip3 install numpy scipy 
RUN pip3 install h5py pandas ArgumentParser
RUN pip3 install lxml

#copy contents of tasks/<taskname>/build/src on the build host into /opt/src on the docker
COPY src/ /opt/src/

WORKDIR /opt/src

