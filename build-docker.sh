#!/bin/bash

set -e

VERSION=$(git describe --tags `git rev-list --tags --max-count=1` --always)

# dynamically pull more interesting stuff from latest git commit
HASH=$(git show-ref --head --hash=8 head)  # first 8 letters of hash should be enough; that's what GitHub uses

# Build standalone docker image
docker build -t quartet-rseqc-report:${VERSION}-${HASH} . && \
docker tag quartet-rseqc-report:${VERSION}-${HASH} ghcr.io/chinese-quartet/quartet-rseqc-report:${VERSION}-${HASH} && \
docker push ghcr.io/chinese-quartet/quartet-rseqc-report:${VERSION}-${HASH}
