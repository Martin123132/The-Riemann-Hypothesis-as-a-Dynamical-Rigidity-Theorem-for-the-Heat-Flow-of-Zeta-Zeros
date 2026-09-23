# Research Snapshot: 24 September 2026

Copyright (c) 2026 TWO HANDS NETWORK LTD. All rights reserved.
Prior written permission is required under the repository's existing licence.

This publication brings the completed local research through the sealed
**Uniform Filtered First Moment** stage of 23 September 2026 into the repository.
It does not include the separate next-stage work still running at publication
time. It does not claim a proof of RH or of the remaining favourable signed bound.

## Start Here

- [Current report](latest/REPORT.md)
- [Complete current derivation](latest/DERIVATION.md)
- [Current audit and limitations](latest/AUDIT.md)
- [Current verification results](latest/VERIFICATION.json)
- [Index of recorded stages](INDEX.md)
- [Coverage and publication verification](PUBLICATION_VERIFICATION.json)

The latest deduction, under the explicitly retained reciprocal-phase estimate,
gives uniform full-filter first-moment decay for the unchanged literal rough
source. It bounds nonlinear overlap by a squared vanishing allowance without
the older first-moment reassembly dependency. The required signed alignment of
the energy-bearing nonlinear response with its original multiplier is still open.

Earlier positive, negative, conditional, finite and superseded results remain
historical evidence. Their recorded classifications have not been rescored.
This update checks publication integrity; it does not rerun every historical
source integration or independently reprove every analytic statement.

## Complete Files, Not Excerpts

The evidence bundles are ordinary ZIP files containing complete source files,
identified by SHA-256. Duplicate copies share one blob. `CATALOG.json` maps
each original logical path to its complete bytes and records original archive
identities and every member. `SOURCE_ROSTER.json` freezes the 90-stage input
selection. No mathematical file has been shortened to fit a display or upload.

The original sealed packages remain unchanged locally. These are newly assembled
publication bundles, not byte-identical copies of the original ZIP containers.
Original member bytes and hashes are retained. Session bookmarks, runtime logs,
PID files, bytecode and local runtime/environment settings are deliberately kept
local; every omission is listed with its hash and reason in the catalog. Original
ZIP containers are inventoried recursively rather than duplicated as opaque blobs.

Some historical scripts retain their original local paths and runtime dependencies.
Hash verification is portable; automatic execution of all historical research on
a new computer is not promised. The published verifier never runs research code.

## Verify or Restore

From this directory, using Python 3:

```text
python verify_snapshot.py
```

This checks the publication manifest, every ZIP member's decompressed SHA-256,
all mapped byte lengths, complete archive-member accounting, and readable reports.
The large-file check reads full binary streams, not tool previews or excerpts.

To restore all included pinned workspace files into a NEW directory:

```text
python verify_snapshot.py --restore restored --prefix workspace/
```

Restoration uses exclusive file creation and rejects path escapes. The output
also includes separately inventoried archive-member trees. Use a narrower prefix
from `CATALOG.json` to restore one package. Nothing is executed automatically.

Repository history and licence documents are preserved. Public visibility does
not grant reuse permission beyond the existing licence and applicable rights.
