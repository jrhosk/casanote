FROM python:3.6

RUN adduser casauser
USER casauser
WORKDIR /home/casauser

ENV PATH="/home/casauser/.local/bin:${PATH}"

RUN python -m pip install --upgrade pip
RUN python -m pip install --no-cache jupyterlab
RUN python -m pip install matplotlib numpy astropy panel param scipy
RUN python -m pip install --index-url https://test.pypi.org/simple/ --no-deps stktools
RUN python -m pip install casatools casatasks casadata casatestutils

COPY --chown=casauser:casauser . .

