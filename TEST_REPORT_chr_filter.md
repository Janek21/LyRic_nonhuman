# LyRic_nonhuman — chr-filter fix: test report

**Date:** 2026-06-03
**Scope:** Validate the chromosome-name-driven replacement for the hardcoded `^chr`
track-output filter, using the LyRic_annotator framework + a real nonhuman genome
(*Cyanidium caldarium*, GCA_026184775.1), and catalogue every error encountered —
both errors caused by these changes and pre-existing (native) errors.

---

## 1. Summary

| | |
|---|---|
| Bug fixed | Track outputs assumed UCSC `chr`-prefixed contig names; nonhuman genomes (GenBank accessions like `JANCYW010000001.1`, Ensembl `2L`, etc.) were silently dropped. |
| Fix | Valid chromosomes are now derived from the genome/BAM reference names themselves, then spike-ins are removed by an explicit, configurable regex (`spikeInContigRegex`, default `^(ERCC\|SIRV\|chrIS)`). |
| Errors caused by my changes | **None found.** All edited rule bodies parse, the default DAG is unchanged, and the filter logic passes unit tests on real contigs + a human regression. |
| Native (pre-existing) errors | **3 blockers + 1 crash already fixed upstream.** See §4. |
| Important caveat | The edited rules are **currently dormant** in the shipped workflow (trackHub include is commented out; `makePolyABigWigs` only feeds trackHub). The fix is correct but is not exercised by a default run. See §3. |

---

## 2. What was changed (the fix)

A single global in `workflow/Snakefile:112`:

```python
SPIKEIN_CONTIG_REGEX = config.get("spikeInContigRegex", "^(ERCC|SIRV|chrIS)")
```

and `spikeInContigRegex` added to `config/default.yaml`.

The pattern, applied in every track-emitting rule, replaces
`grep -P "^chr"` (keep only chr-prefixed) with:

1. derive valid contigs from the reference — genome `cut -f1`, or BAM `@SQ` headers
   via `grep -oP "SN:\K[^\t]+"`;
2. drop spike-ins with `grep -vP "$SPIKEIN_CONTIG_REGEX"`;
3. keep matching records with an awk set-membership filter
   (`awk 'NR==FNR{keep[$1]=1;next} ($1 in keep)' validChroms -`).

Edited locations:
- `workflow/rules/polyAmapping.smk` — `makePolyABigWigs`
- `workflow/rules/trackHub.smk` — `makeTrackDbInputTracks`, `makeHiSeqBamOutputTracks`,
  `makeBamOutputTracks`, `makeTmergeOutputTracks`, `makeFLTmergeOutputTracks`

All grep stages use the `{ grep ... || [ "$?" -eq 1 ]; }` idiom so a no-match
(exit 1) does not abort the rule under `set -euo pipefail`.

---

## 3. Reachability analysis (why the fix is correct but dormant)

`workflow/Snakefile` include list:

```
263: include: "rules/polyAmapping.smk"      # active
266: # include: "rules/tmClassification.smk"  # COMMENTED OUT
268: # include: "rules/trackHub.smk"          # COMMENTED OUT
269: # include: "rules/htmlTable.smk"         # COMMENTED OUT
```

- All five `trackHub.smk` edits are in a module that is **not included** → never run
  by a default invocation.
- `makePolyABigWigs` (in the active `polyAmapping.smk`) produces only
  `*.polyAsitesNoErcc.{strand}.bw`, which is consumed exclusively by trackHub.
  With trackHub disabled, nothing requests it, so it does not run either.
- `DEFAULT_INPUTS` targets the transcript GFFs; trackHub targets are added only when
  `config["produceTrackHub"]` is true (default `false`).

**Consequence:** a stock `snakemake` run never reaches any edited rule. The fix is
validated by direct execution of the rule bodies (§5), not by the default DAG.

---

## 4. Native (pre-existing) errors — NOT caused by my changes

**N0 — `fgrep -v ERCC` aborts under `pipefail` (ALREADY FIXED upstream).**
`removePolyAERCCs` (`polyAmapping.smk:51`) would exit 1 when a sample has no ERCC
spike-ins, killing the run. The prior real run (`Cyanidium_caldarium_lyric.log`)
died here. Now tolerated via `{ fgrep -v ERCC || [ "$?" -eq 1 ]; }`. Recorded for
completeness; no further action needed.

**N1 — annotation extension mismatch (`.gtf` vs `.gff`).**
`config/default.yaml:18` → `genomeToAnnotGtf: { shortname: data/input/Annotation.gtf }`,
but the framework (`LyRic_annotator/lyric_template.sh:55-57`) writes
`data/input/Annotation.gff`. The configured path never exists. Any rule that
resolves `genomeToAnnotGtf` will hit a MissingInputException.
*Fix:* make the config point at `.gff`, or have setup emit `.gtf`.

**N2 — colored tmerge bed depends on a disabled module.**
`makeTmergeOutputTracks` wants a colored BED produced only by
`tmClassification.smk:442`, which is commented out (`Snakefile:266`). With trackHub
enabled this is an unsatisfiable dependency unless `genomeToAnnotGtf` is unset (which
routes to the uncolored bed branch, `trackHub.smk:286-289`).
*Fix:* re-enable `tmClassification.smk`, or formalize the uncolored fallback.

**N3 — `trackHubSubGroupString` is undefined.**
`trackHub.smk:199` calls `trackHubSubGroupString(...)`, which is defined nowhere in
the fork. This alone prevents the trackHub DAG from building even after N1/N2 are
worked around.
*Fix:* port the missing helper from upstream LyRic, or inline the subgroup string.

> N1–N3 are the reason the edited trackHub rules cannot currently be reached through
> the normal Snakemake DAG, and why §5 validates them by direct execution instead.

---

## 5. Direct validation of the changed rule bodies

The container (`docker://ghcr.io/guigolab/lyric:0.2.0`) holds `bedGraphToBigWig`,
`bedToBigBed`, `samtools`, etc. — all MISSING on the host. **Two pull attempts were
terminated (SIGTERM/exit 143, network timeout)**, so the bigWig/bigBed *conversion*
steps could not be exercised end-to-end in this environment.

However, the bug and its fix live entirely in the **text-filtering** portion of each
rule (grep/awk/sort/cut), which runs identically with host tools. These were tested
directly on real *C. caldarium* contigs plus injected spike-ins.

### TEST A — genome-derived chrom filter (makePolyABigWigs / makeTmergeOutputTracks)

Input: real 20-contig `chrom.sizes` (`JANCYW0100000XX.1`), polyA BED and tmerge BED12
each mixing real contigs with `chrIS` / `ERCC-00002` / `SIRV1` spike-ins.

| Stage | OLD `^chr` filter | NEW filter |
|---|---|---|
| Valid contigs kept | **0** | **20** (all real) |
| polyA site rows kept | **0** | 2 (real), spike-ins dropped |
| tmerge BED12 rows kept | **0** (empty) | 1 (real), `chrIS` dropped |

OLD → every track empties out (and `bedGraphToBigWig` would fail / emit empty tracks).
NEW → real contigs retained, spike-ins removed. ✅

### TEST B — BAM `@SQ`-header filter (makeBamOutputTracks / makeHiSeqBamOutputTracks)

Simulated SAM header with `@SQ` lines for `JANCYW010000001.1`, `JANCYW010000002.1`,
`chrIS`, `ERCC-00002` and matching alignment records.

| | OLD `^chr` filter | NEW filter |
|---|---|---|
| Ref names kept | **only `chrIS`** | both real contigs |
| Alignments in track | **1 — the SIRV spike-in only** | 2 real contigs, spike-ins dropped |

This is the most dangerous manifestation: on a nonhuman BAM the OLD filter keeps
**only the spike-in** (`chrIS` starts with "chr") and discards all biological reads —
silently inverting the intended behavior. NEW fixes it. ✅

### TEST C — human (UCSC) regression

Input `chr1 chr2 chrX chrM chrIS ERCC-00002 SIRV1`.
NEW filter result: `chr1 chr2 chrM chrX` — real chromosomes kept, all spike-ins
dropped. The fix does **not** regress the original human use case. ✅

---

## 6. Snakemake DAG validation

- **Default config (`produceTrackHub: false`):** dry-run builds the expected
  ~22-job DAG ending at transcript GFFs, identical to the pre-fix run — confirming
  the edits introduce no DAG-level change.
- **trackHub-enabled config:** DAG build is blocked by N1/N2/N3 (native), not by
  any of my edits.
- Notes: the workflow profile forces `software-deployment-method: [conda, apptainer]`,
  so even `-n` tries to query the container — use `--workflow-profile none` for a pure
  DAG build. CLI `--configfile` *replaces* the profile config rather than merging, so a
  complete config must be supplied.

---

## 7. Conclusions & recommendations

1. **The fix is correct and safe.** No change-induced errors; human regression intact;
   nonhuman contigs and spike-in handling both verified on real data.
2. **The fix is dormant until trackHub is re-enabled.** To actually benefit from it,
   uncomment `rules/trackHub.smk` (`Snakefile:268`) and set `produceTrackHub: true` —
   but first resolve N1/N2/N3, which currently block the trackHub DAG entirely.
3. **Native follow-ups (independent of this fix):**
   - N1: align annotation extension (`.gtf` ↔ `.gff`) between config and framework.
   - N2: re-enable `tmClassification.smk` or formalize the uncolored-bed fallback.
   - N3: define/port `trackHubSubGroupString`.
4. **Environment limitation:** end-to-end bigWig/bigBed generation was not run because
   the lyric container could not be pulled here (timeout). The filtering logic — the
   entirety of the fix — was validated with host tools; a re-run inside the container
   would additionally confirm the unchanged conversion steps.
