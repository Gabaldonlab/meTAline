#!/bin/bash

# Check if spython is installed
if python3 -c 'import pkgutil; exit(0) if pkgutil.find_loader("spython") else exit(1)'; then
    echo "==> spython is already installed."
else
    echo "==> spython is not installed. Installing..."
    pip3 install spython
fi

echo "==> translating ./Dockerfile to metaline-singularity.def"
spython recipe --entrypoint /bin/sh ./Dockerfile > metaline-singularity.def

echo "==> executing post-fixs"

# Fix the translation error due the difference between how Singularity and Docker copies the files.
sed -i 's|./external-sources /bin/|./external-sources/* /bin/|g' metaline-singularity.def

# Remove the default runscript and startscript sections.
sed -i '/^%runscript$/,$d' metaline-singularity.def
