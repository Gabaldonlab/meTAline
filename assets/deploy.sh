#!/bin/bash

CLUSTER_USER=""
CLUSTER_HOST=""
CLUSTER_PATH=""
CURRENT_VERSION_TAG="v1.2"

rm -rf metaline.sif
docker image rm -f cgenomics/metaline:latest cgenomics/metaline:${CURRENT_VERSION_TAG}

set -e

docker login
docker image build -t cgenomics/metaline:latest . && docker push cgenomics/metaline:latest &
docker image build -t cgenomics/metaline:${CURRENT_VERSION_TAG} . && docker push cgenomics/metaline:${CURRENT_VERSION_TAG} &

# Sign up & login to Sylabs
# 1. Go to https://cloud.sylabs.io/
# 2. Create an account or sign in.
# 3. Generate an authentication token:
#    - Go to https://cloud.sylabs.io/auth/tokens.
#    - Click "Generate New Token" and copy it.

# 4. Step 2: Login via Singularity CLI
singularity remote login && \
make singularity && \
rsync -vaP ./metaline.sif ${CLUSTER_USER}@${CLUSTER_HOST}:${CLUSTER_PATH}

singularity push -U metaline.sif library://cgenomics/metaline/metaline:${CURRENT_VERSION_TAG}

wait
