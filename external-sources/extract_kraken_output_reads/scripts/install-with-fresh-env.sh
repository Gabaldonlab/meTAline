#!/bin/sh

deactivate \
    && rm -rf ./extract-kraken-output-reads-env \
    && python3 -m venv extract-kraken-output-reads-env \
    && source ./extract-kraken-output-reads-env/bin/activate \
    && pip install -e .