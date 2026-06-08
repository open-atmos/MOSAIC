# Dockerfile to build the MOSAIC box model (MOSAIC_v2).
#
# Build:  docker build -t mosaic .
# Run:    docker run --rm -it mosaic        (mosaic.x is on the PATH)

# Pull Fedora from its canonical registry rather than Docker Hub, which
# intermittently returns 500 errors (e.g. in GitHub Actions).
FROM registry.fedoraproject.org/fedora:41

RUN dnf -y update \
    && dnf -y install \
        less \
        tmux \
        make \
        gcc \
        gcc-gfortran \
        tcsh \
    && dnf clean all

# NOTE: see .dockerignore for which files/directories are copied.
COPY . /mosaic/

WORKDIR /mosaic

# Build mosaic.x and libmosaic.a via the project's build script (gfortran).
# compile.sh stages sources into compile/, which .dockerignore strips on COPY,
# so create it first.
RUN mkdir -p compile \
    && ./compile.sh gfortran

ENV PATH="/mosaic/compile:${PATH}"

CMD ["/bin/bash"]
