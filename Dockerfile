FROM python:3.6

WORKDIR /

RUN python -m pip install --upgrade pip
RUN python -m pip install --no-cache jupyterlab
RUN python -m pip install matplotlib numpy astropy panel param scipy
RUN python -m pip install --index-url https://test.pypi.org/simple/ --no-deps stktools
RUN python -m pip install --index-url https://casa-pip.nrao.edu/repository/pypi-casa-release/simple casatools casatasks casadata casatestutils

