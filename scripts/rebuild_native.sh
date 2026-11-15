#!/bin/bash
# Rebuild and install arcsv2_native into the active micromamba env.
# Uses system OpenSSL (OPENSSL_NO_VENDOR=1) — vendor builds need FindBin.
# See docs/R1_STATUS.md.
set -euo pipefail
cd "$(dirname "$0")/.."
export PYTHONUNBUFFERED=1
mkdir -p target/wheels

if ! find /opt/conda -name "libclang.so*" 2>/dev/null | head -1 | grep -q .; then
  echo "installing clang/libclang via micromamba..."
  micromamba install -y -n base -c conda-forge clangdev libclang
fi

export OPENSSL_NO_VENDOR=1
export PKG_CONFIG_PATH="/opt/conda/lib/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"
if [[ -z "${LIBCLANG_PATH:-}" ]]; then
  for cand in /opt/conda/lib /usr/lib/llvm-*/lib /usr/lib/x86_64-linux-gnu; do
    for so in ${cand}/libclang.so*; do
      if [[ -e "$so" ]]; then
        export LIBCLANG_PATH="$(dirname "$so")"
        break 2
      fi
    done
  done
fi
if [[ -z "${LIBCLANG_PATH:-}" ]]; then
  found=$(find /opt/conda -name "libclang.so*" 2>/dev/null | head -1 || true)
  if [[ -n "$found" ]]; then
    export LIBCLANG_PATH="$(dirname "$found")"
  fi
fi

echo "OPENSSL_NO_VENDOR=$OPENSSL_NO_VENDOR LIBCLANG_PATH=${LIBCLANG_PATH:-unset} patchelf=$(command -v patchelf || echo missing)"
echo "openssl.pc: $(pkg-config --modversion openssl 2>/dev/null || echo missing)"

OPENSSL_NO_VENDOR=1 micromamba run -n base \
  env OPENSSL_NO_VENDOR=1 \
      LIBCLANG_PATH="${LIBCLANG_PATH:-}" \
      PKG_CONFIG_PATH="${PKG_CONFIG_PATH:-}" \
      PATH="/opt/conda/bin:$PATH" \
  maturin build --release -m native/arcsv2-native/Cargo.toml \
    -o target/wheels

WHEEL=$(ls -t target/wheels/arcsv2_native-*.whl | head -1)
echo "WHEEL=$WHEEL"
command -v patchelf >/dev/null || echo "WARN: patchelf not on PATH"
micromamba run -n base pip install --force-reinstall --no-deps "$WHEEL"
micromamba run -n base python -c "from arcsv2 import native; print(native.native_version(), native.rust_available())"
