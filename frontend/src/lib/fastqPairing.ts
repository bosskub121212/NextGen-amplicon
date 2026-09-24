// Read-direction token handling for paired-end FASTQ filenames.
//
// The token that says "this is read 1 / read 2" only ever sits at the END of a
// FASTQ filename, immediately before the extension:
//
//     Sample_R1.fastq.gz     Sample_R1_001.fastq.gz     Sample_1.fq.gz
//
// Matching "_R1" ANYWHERE in the name — which is what this app used to do —
// breaks every dataset whose SAMPLE happens to be called R1 or R2. For
// Bacteria_R1_1.fq.gz the real token is the trailing "_1"; "_R1" is part of the
// sample name. Two separate failures came out of that:
//
//   1. the sample name collapsed to "Bacteria", so both of its files were
//      treated as forward reads of one unpaired sample, and a second sample
//      could end up sharing the same name;
//   2. worse, the mate guess rewrote "Bacteria_R1_1.fq.gz" into
//      "Bacteria_R2_2.fq.gz" (replace _R1 -> _R2, then _1. -> _2.), pairing one
//      sample's forward reads with a DIFFERENT sample's reverse reads. That
//      does not fail loudly; it produces a plausible-looking ASV table.
//
// So: anchor the token at the end of the name. Keep the old loose rule only as
// a fallback for layouts that genuinely carry it elsewhere (Trim Galore writes
// Sample_R1_val_1.fq.gz), and only when the anchored form matches nothing.

const EXT = /\.(fastq|fq)(\.gz)?$/i;

/** `_R1` / `_1` / `_R1_001` immediately before the extension. */
export const FWD_TOKEN = /_(R?)1((?:_\d{3})?)\.(fastq|fq)(\.gz)?$/i;
/** `_R2` / `_2` / `_R2_001` immediately before the extension. */
export const REV_TOKEN = /_(R?)2((?:_\d{3})?)\.(fastq|fq)(\.gz)?$/i;

const FWD_LOOSE = /_R1/i;
const REV_LOOSE = /_R2/i;

export type Pairing = {
  /** true when the anchored token was found somewhere in the file set */
  strict: boolean;
  isR1: (name: string) => boolean;
  isR2: (name: string) => boolean;
  /** filename -> sample name */
  sampleOf: (name: string) => string;
  /** forward filename -> its expected reverse filename */
  mateOf: (name: string) => string;
};

/**
 * Build the pairing rules for one set of filenames. The strict/loose decision is
 * made once for the whole set, so a single odd name cannot flip the rule for the
 * others.
 */
export function pairingFor(names: string[]): Pairing {
  const fastq = names.filter((n) => EXT.test(n));
  const hasToken = fastq.some((n) => FWD_TOKEN.test(n) || REV_TOKEN.test(n));

  // The anchored rule is right for every normal layout, but some trimmers put a
  // token at BOTH ends — Trim Galore writes S1_R1_val_1.fq.gz / S1_R2_val_2.fq.gz,
  // where only the leading one identifies the mate. So verify the anchored rule
  // against the actual file list first: if it leaves mates unfound and the loose
  // rule finds them all, the loose rule is the correct one for this set.
  const strictMate = (n: string) =>
    n.replace(
      FWD_TOKEN,
      (_m, r: string, lane: string, ext: string, gz: string) =>
        `_${r}2${lane || ""}.${ext}${gz || ""}`
    );
  // Legacy mate rewrite: swap the leading _R1 AND a trailing _1. — the pair of
  // substitutions the app used to apply unconditionally. It is only ever used
  // when the anchored rule has already been shown not to pair this file set.
  const looseMate = (n: string) =>
    n.replace(FWD_LOOSE, "_R2").replace(/_1\.(fq|fastq)/i, "_2.$1");

  let strict = hasToken;
  if (hasToken) {
    const fwd = fastq.filter((n) => FWD_TOKEN.test(n));
    const strictMisses = fwd.filter((n) => !fastq.includes(strictMate(n))).length;
    if (strictMisses > 0) {
      const looseFwd = fastq.filter((n) => FWD_LOOSE.test(n));
      const looseMisses = looseFwd.filter((n) => !fastq.includes(looseMate(n))).length;
      if (looseFwd.length > 0 && looseMisses === 0) strict = false;
    }
  }

  if (strict) {
    return {
      strict: true,
      isR1: (n) => FWD_TOKEN.test(n),
      isR2: (n) => REV_TOKEN.test(n),
      sampleOf: (n) =>
        n.replace(FWD_TOKEN, "").replace(REV_TOKEN, "").replace(EXT, ""),
      // Rebuild the SAME token shape with a 2 in it, so _R1_001.fastq.gz maps to
      // _R2_001.fastq.gz and _1.fq.gz maps to _2.fq.gz — and nothing earlier in
      // the name is touched.
      mateOf: strictMate,
    };
  }

  return {
    strict: false,
    isR1: (n) => FWD_LOOSE.test(n),
    isR2: (n) => REV_LOOSE.test(n),
    sampleOf: (n) => n.replace(/_R1.*|_R2.*/i, "").replace(EXT, ""),
    mateOf: looseMate,
  };
}

/** Sample name for single-end / ONT mode: the filename without its extension. */
export const sampleOfSingle = (name: string) => name.replace(EXT, "");

export type PairingProblem = { kind: "duplicate" | "unpaired"; samples: string[] };

/**
 * Sanity-check a detected sample list so the UI can say something before a run
 * starts, instead of the pipeline failing eight minutes in.
 */
export function checkPairing(names: string[], fileNames: string[]): PairingProblem[] {
  const problems: PairingProblem[] = [];

  const seen = new Map<string, number>();
  names.forEach((n) => seen.set(n, (seen.get(n) || 0) + 1));
  const dup = [...seen.entries()].filter(([, c]) => c > 1).map(([n]) => n);
  if (dup.length) problems.push({ kind: "duplicate", samples: dup });

  const P = pairingFor(fileNames);
  const fastq = fileNames.filter((f) => EXT.test(f));
  const unpaired = fastq
    .filter((f) => P.isR1(f))
    .filter((f) => !fastq.includes(P.mateOf(f)))
    .map((f) => P.sampleOf(f));
  if (unpaired.length) problems.push({ kind: "unpaired", samples: unpaired });

  return problems;
}
