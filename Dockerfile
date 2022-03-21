FROM python:3.6

ENV HOME /tmp
WORKDIR /tmp

RUN python -m pip install --upgrade pip
RUN python -m pip install --no-cache notebook 
RUN python -m pip install matplotlib numpy astropy panel param scipy
RUN python -m pip install --index-url https://test.pypi.org/simple/ --no-deps stktools
RUN python -m pip install --index-url https://casa-pip.nrao.edu/repository/pypi-casa-release/ casatools casatasks casadata casatestutils
RUN apt-get update && apt-get -y install git

CMD ['bash']
