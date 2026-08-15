Barcode correction
==================

``generate-permit-list`` (GPL) decides barcode corrections.  Later stages do
not search barcode neighbourhoods or reconsider collision policies: ``collate``
and the ATAC sort path apply the decisions GPL compiled.  This division keeps
the biological policy in alevin-fry and leaves RAD parsing, rewriting,
spooling, and bucket routing in libradicl.

Correction policies
-------------------

The cell-barcode policy is selected with
``--cell-bc-correction unique|frequency``.

``unique``
  Accept a non-exact observation only if all neighbouring retained sources
  identify one distinct canonical target.  Multiple Flex construction
  sequences for the same canonical sample count as one target, not as a
  collision.

``frequency``
  Score the possible canonical targets using raw, exact source counts frozen
  before any correction.  Each candidate source contributes its exact count
  plus a pseudocount of one, and contributions from construction aliases are
  summed by canonical target.  The best target must meet the confidence
  threshold.  Confidence is compared as an exact rational value, and any
  remaining equal-score tie is resolved by the packed target ordering.  The
  default confidence is 97.5% for RNA and probe/sample barcodes and 90% for
  ATAC.  It can be changed with ``--cell-bc-confidence``.

An exact permitted source always wins, regardless of neighbouring sources.
Thus correction is deterministic under whitelist order, hash iteration order,
worker completion order, and thread count.

Neighbourhoods
--------------

The neighbourhood is selected with
``--cell-bc-neighborhood hamming-1|substitution-or-shift-1``.

``hamming-1``
  One nucleotide substitution at a fixed barcode position.

``substitution-or-shift-1``
  One substitution, or the historical fixed-length one-base shift operation.
  A shift drops a base at one end and admits a base at the other.  This is not
  general Levenshtein edit distance.  ``edit-1`` remains accepted as a
  compatibility spelling for this neighbourhood.

Shift-aware Frequency correction is permitted, but alevin-fry emits a warning:
RAD stores neither barcode base qualities nor an indel-error model, so this
combination is an abundance heuristic rather than a calibrated indel
posterior.

When the user does not select a neighbourhood, the resolved defaults preserve
the existing protocol behaviour:

==============================  =============================
Workflow                        Default neighbourhood
==============================  =============================
Filtered ordinary RNA           substitution-or-shift-1
Unfiltered ordinary RNA         hamming-1
Filtered multi-sample RNA       hamming-1
Unfiltered multi-sample RNA     hamming-1
ATAC                            hamming-1
==============================  =============================

All RNA filtering methods use the same policy engine: knee, expected cells,
forced cells, an explicit filtered list, and an unfiltered external list first
select retained targets and then compile observed corrections from frozen raw
counts.

Sample and probe barcodes
-------------------------

For a multi-barcode RAD file, ``--sample-bc-correction`` accepts ``exact``,
``unique``, or ``frequency`` and defaults to ``exact``.  Non-exact modes use
``--sample-bc-neighborhood`` (Hamming-1 by default), and sample Frequency uses
``--sample-bc-confidence``.  The deprecated
``--sample-correction-mode 1-edit`` spelling resolves to Unique plus
``substitution-or-shift-1``; ``exact`` retains its previous meaning.

Flex construction sequences are resolved to their canonical sample before a
collision policy is applied.  Sample Frequency still makes one pass over the
RAD file:

1. Exact matches and observations having one structural target are routed
   immediately; observations with no target are rejected immediately.
2. Only observations with more than one canonical target are deferred as
   ``(raw sample, cell)`` pairs.
3. Each worker holds a bounded buffer.  A full buffer is sorted, run-length
   encoded, and appended to that worker's Snappy-compressed temporary spool.
4. After exact sample priors are frozen, ambiguous sample barcodes are resolved
   and the runs are replayed additively into the cell histograms.  No global
   run merge is required.

Exact and Unique sample correction create no correction spools.  GPL's
``--memory-limit`` controls the total deferred-pair buffer budget and defaults
to 512 MiB.  Values below 256 MiB produce a warning and use 256 MiB.
``--tmp-dir`` chooses the spool directory and defaults to GPL's output
directory.  Spools are removed on success and on error unwinding.

Compiled correction plan
------------------------

GPL writes ``correction_plan.bin`` in addition to the historical
``permit_map.bin``, ``sample_permit_map.bin``, and ``permit_freq.bin`` files.
The plan has a magic header and format version, declares the sample and cell
barcode lengths and resolved policies, and stores sorted observed-to-corrected
entries including retained identities.  Multi-sample cell corrections are
grouped by canonical sample.

Collation validates the plan, compiles its sorted entries into a direct fused
correction-and-bucket lookup, and performs one lookup on each RAD record.  The
record hot path contains no neighbour generation, collision-policy dispatch,
or callback.  Multi-sample output remains hierarchical: records for a sample
are contiguous, and records for each corrected cell within that sample are
collated together.

Backward compatibility is explicit:

* ordinary RNA and ATAC outputs without a plan use ``permit_map.bin``;
* multi-sample outputs without a plan use the historical Hamming correction
  and warn that GPL must be rerun to reproduce newer policies;
* a plan that is present but malformed or has an unsupported version is an
  error, rather than silently selecting a different correction path; and
* corrections whose target is absent from the output permit list remain
  excluded.

``permit_freq.bin`` contains counts aggregated by corrected target.  Existing
``sample_info.json`` fields retain their meanings; ``collatable_reads`` and
the correction/routing diagnostics are recorded separately.

Resource floors
---------------

alevin-fry uses two threads as its practical minimum.  A request for zero or
one thread, or an available-parallelism report of one, produces a warning and
continues with two threads.  Requests for two or more threads are honored.
