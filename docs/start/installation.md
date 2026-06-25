# Installation

XGrammar Python Package can be installed directly from a prebuilt package or built from source.

## Method 1: Prebuilt Package

XGrammar supports various platforms:
* Operating Systems: Linux, macOS, and Windows
* Hardware: CPU, NVIDIA GPUs, AMD GPUs, Apple Silicon, TPU, etc.
* Python: 3.9 and later.

We provide Python wheels for XGrammar via pip.

```bash
python -m pip install xgrammar
```

For use with MPS on Apple Silicon, install with:

```bash
python -m pip install "xgrammar[metal]"
```

We also provide conda packages for XGrammar:

```bash
conda install -c conda-forge xgrammar
```

Use the following command to verify installation:

```bash
python -c "import xgrammar; print(xgrammar)"
# Prints: <module 'xgrammar' from '/path-to-env/lib/python3.11/site-packages/xgrammar/__init__.py'>
```

### Optional Dependency: PyTorch

PyTorch is an **optional** dependency. The C++ core operates on
[DLPack](https://github.com/dmlc/dlpack) tensors, so `import xgrammar` and all core
functionality (grammar compilation, the grammar matcher, tokenizer info, structural tags,
JSON-schema generation, etc.) work without PyTorch installed.

When PyTorch is **not** installed:

* `allocate_token_bitmask` / `reset_token_bitmask` return and operate on **NumPy** arrays.
* `fill_next_token_bitmask` accepts any array that supports the DLPack protocol (including
  NumPy arrays), so the matcher can fill bitmasks without PyTorch.

Install the `torch` extra if you need any of the PyTorch-backed features:

```bash
python -m pip install "xgrammar[torch]"
```

This is required for:

* `allocate_token_bitmask` returning a `torch.Tensor` (convenient for inference engines that
  manage logits on GPU).
* `apply_token_bitmask_inplace` on real model logits via the CPU / Triton / Torch kernels.
* the HuggingFace `model.generate()` integration `xgrammar.contrib.hf.LogitsProcessor`.
* the PyTorch-based helpers in `xgrammar.testing`.

> When PyTorch is installed, `allocate_token_bitmask` returns a `torch.Tensor` exactly as
> before, so existing code keeps working unchanged.

## Method 2: Build XGrammar Python Package from Source

This option is useful when you want to make modification or obtain a specific version of XGrammar.

```bash
git clone --recursive https://github.com/mlc-ai/xgrammar.git && cd xgrammar
pre-commit install
# Copy cmake config. You can update the config if needed.
cp cmake/config.cmake .
# Install scikit-build-core and apache-tvm-ffi as build dependencies.
python3 -m pip install scikit-build-core apache-tvm-ffi
python3 -m pip install --no-build-isolation -e .
```

XGrammar is a library written in C++ and Python. The editable install will automatically rebuild
the package when XGrammar is imported in Python.

### Optional: Run Python Tests

```bash
# Install the test dependencies
python3 -m pip install ".[test]"

# If you have a HuggingFace token, you can run all tests including the ones that have gated models.
huggingface-cli login --token YOUR_HF_TOKEN
python3 -m pytest

# If you do not have a HuggingFace token, you can run a subset of tests that do not require gated models.
python3 -m pytest -m "not hf_token_required"
```

## Method 3: Build XGrammar C++ Library Only

XGrammar can also be built as a C++ library. This is useful for using XGrammar in C++ or Rust projects.

XGrammar uses CMake and Ninja to build the C++ library. To build only the C++ library, you can set -DXGRAMMAR_BUILD_PYTHON_BINDINGS=OFF when running CMake, or modify cmake/config.cmake manually.

```bash
git clone --recursive https://github.com/mlc-ai/xgrammar.git && cd xgrammar
# Copy cmake config. You can update the config if needed.
cp cmake/config.cmake .
mkdir build && cd build
cmake -G Ninja ..
ninja
```

### Optional: Run C++ Tests

```bash
# Run all tests
bash scripts/run_ctest.sh

# Run a subset of tests whose name contains "test_name"
bash scripts/run_ctest.sh test_name
```
