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
    && dnf clean all

# NOTE: see .dockerignore for which files/directories are copied.
COPY . /mosaic/

WORKDIR /mosaic

# Stage all sources into compile/ and build mosaic.x and libmosaic.a
# (mirrors compile.sh).
RUN mkdir -p compile \
    && cd compile \
    && rm -f *.f90 *.f Makefile* \
    && cp -p ../datamodules/*.f90 . \
    && cp -p ../main/*.f90 . \
    && cp -p ../gas/*.f90 . \
    && cp -p ../solver/*.f . \
    && cp -p ../aerosol/*.f90 . \
    && cp -p ../cloud/*.f90 . \
    && cp -p ../Makefile . \
    # module_mosaic_support.f90 must pass through the C preprocessor first
    && rm -f module_mosaic_support.f90 \
    && cpp -P -traditional-cpp ../main/module_mosaic_support.f90 module_mosaic_support.f90 \
    && make mosaic.x   COMP=gfortran FLAGS="-g -fallow-argument-mismatch" --always-make \
    && make libmosaic.a COMP=gfortran FLAGS="-g -fallow-argument-mismatch" --always-make

ENV PATH="/mosaic/compile:${PATH}"

CMD ["/bin/bash"]
