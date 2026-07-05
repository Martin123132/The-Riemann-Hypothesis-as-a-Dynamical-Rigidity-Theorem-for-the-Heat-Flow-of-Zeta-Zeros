# Environment Baseline

Date: 2026-07-03

Status: environment inventory. This is not a proof artifact; it records local compute resources for reproducible numerical work.

Purpose:

Record the available local compute environment for the RH proof-programme numerics.

## System

```text
OS/platform: Windows-10-10.0.19045-SP0
logical processors: 12
physical memory: 34,194,378,752 bytes, about 31.8 GiB
```

## Storage

```text
C: free about 1730.66 GB
D: free about 1195.96 GB
```

## Python

```text
Python: 3.13.13
Executable: C:\Users\ollet\AppData\Local\Programs\Python\Python313\python.exe
```

## Numeric Libraries

```text
mpmath: 1.3.0
numpy: 2.4.6
scipy: 1.17.1
sympy: 1.14.0
flint/python-flint: missing
gmpy2: missing
psutil: missing
```

## Local Additions

Installed locally under the project workspace:

```text
work/rh_compute/vendor/python-flint 0.8.0
```

Import path requirement:

```python
import sys
sys.path.insert(0, r"C:\Users\ollet\Documents\Codex\2026-07-03\l\work\rh_compute\vendor")
import flint
```

Verified:

```text
flint.arb
flint.acb
flint.arb_mat
```

## Implication

First reproduction runner should use:

```text
mpmath for high-precision moment computation
sympy for exact/symbolic Sturm-style root counting where practical
numpy/scipy only for secondary diagnostics
```

Later rigorous-certification work should consider installing:

```text
gmpy2
```

if we need faster multiprecision support. Arb/FLINT ball arithmetic is now
available locally through the workspace vendor directory.
