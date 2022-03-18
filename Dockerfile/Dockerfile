FROM python:3.6

ENV HOME /root
WORKDIR /root

RUN python -m pip install --upgrade pip
RUN python -m pip install notebook matplotlib numpy astropy panel param scipy
RUN python -m pip install --index-url https://test.pypi.org/simple/ --no-deps stktools
RUN python -m --index-url http://casa-pip.nrao.edu/repository/pypi-casa-release/simple casatools casatasks casadata casatestutils
RUN apt-get update && apt-get -y install git

CMD ['bash']