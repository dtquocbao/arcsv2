# Linux runtime for Windows, macOS, and Linux. Rust + maturin build the
# optional arcsv2-native wheel during image build. See docs/CROSS_PLATFORM.md.
FROM mambaorg/micromamba:2.0.5

LABEL org.opencontainers.image.title="arcsv2"

USER root
COPY environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean --all --yes

COPY . /opt/arcsv2
# OPENSSL_NO_VENDOR=1: link conda libssl (vendor openssl-src needs FindBin).
# PKG_CONFIG_PATH / LIBCLANG_PATH / patchelf: openssl.pc + bindgen + wheel RPATH.
ENV OPENSSL_NO_VENDOR=1 \
    PKG_CONFIG_PATH=/opt/conda/lib/pkgconfig \
    LIBCLANG_PATH=/opt/conda/lib
RUN micromamba run -n base pip install --no-cache-dir "/opt/arcsv2[dev]" \
    && micromamba run -n base env OPENSSL_NO_VENDOR=1 \
        PKG_CONFIG_PATH=/opt/conda/lib/pkgconfig \
        LIBCLANG_PATH=/opt/conda/lib \
        PATH=/opt/conda/bin:$PATH \
        maturin build --release -m /opt/arcsv2/native/arcsv2-native/Cargo.toml \
        -o /tmp/wheels \
    && micromamba run -n base pip install --no-cache-dir /tmp/wheels/arcsv2_native-*.whl \
    && rm -rf /tmp/wheels \
    && mkdir -p /data \
    && chown -R "$MAMBA_USER:$MAMBA_USER" /opt/conda /opt/arcsv2 /data
# Login shells (bash -lc, as in CI) reset PATH. Restore the env bin directory.
RUN printf 'export PATH="/opt/conda/bin:$PATH"\n' > /etc/profile.d/conda-path.sh

USER $MAMBA_USER
WORKDIR /work
# The base image activates conda in its own entrypoint. CI replaces that
# entrypoint with pytest, snakemake, or bash, so the env bin directory has to
# be on PATH itself.
ENV PATH="/opt/conda/bin:$PATH"
ENTRYPOINT ["arcsv2"]
