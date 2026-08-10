#!/usr/bin/env python3
"""Check Linux wheels for unversioned undefined libstdc++ symbols.

On manylinux, symbols newer than the base image's libstdc++ must be linked
statically from gcc-toolset's libstdc++_nonshared.a. If the linker instead
resolves them from another shared library's accidental exports (e.g. the
libtvm_ffi.so shipped by apache-tvm-ffi <= 0.1.12), the wheel is left with
*unversioned* undefined references that only load when some other library
happens to provide them — and fail with ImportError on older-libstdc++ hosts
(RHEL 8, Ubuntu 22.04; see issue #821). auditwheel cannot catch this class of
breakage because its policy check only inspects *versioned* symbol
requirements, so we assert it separately after wheel repair.

Usage: check_wheel_symbols.py WHEEL [WHEEL ...]
Non-Linux wheels are skipped. Exits non-zero if any Linux wheel contains an
unversioned undefined std:: symbol.
"""

import re
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

# Mangled names of std:: entities: _ZSt*, _ZNSt*/_ZNKSt*, and their vtable /
# VTT / typeinfo / guard-variable forms (_ZTVNSt..., _ZGVNSt..., ...).
STD_SYMBOL = re.compile(r"_Z(T[VTIS]|GV)?N?K?St")


def unversioned_und_std_symbols(so_path: Path) -> list:
    out = subprocess.run(
        ["readelf", "--dyn-syms", "--wide", str(so_path)],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    bad = []
    for line in out.splitlines():
        fields = line.split()
        if len(fields) < 8 or fields[6] != "UND":
            continue
        symbol = fields[7]
        # A versioned reference looks like name@GLIBCXX_3.4.21; only bare
        # (unversioned) references are the problem.
        if "@" not in symbol and STD_SYMBOL.match(symbol):
            bad.append(symbol)
    return bad


def check_wheel(wheel: Path) -> bool:
    if not re.search(r"(many|musl)linux", wheel.name):
        print(f"SKIP {wheel.name}: not a Linux wheel")
        return True
    ok = True
    with tempfile.TemporaryDirectory() as tmp, zipfile.ZipFile(wheel) as zf:
        for name in zf.namelist():
            if ".so" not in Path(name).name:
                continue
            so_path = Path(zf.extract(name, tmp))
            bad = unversioned_und_std_symbols(so_path)
            if bad:
                ok = False
                print(f"FAIL {wheel.name}: {name} has unversioned undefined std:: symbols:")
                for symbol in bad:
                    print(f"       {symbol}")
    if ok:
        print(f"OK   {wheel.name}")
    return ok


def main() -> int:
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    results = [check_wheel(Path(arg)) for arg in sys.argv[1:]]
    return 0 if all(results) else 1


if __name__ == "__main__":
    sys.exit(main())
