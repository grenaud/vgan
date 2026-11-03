# Does soibean have the same damage problem as TrailMix? Yes — and it's actually worse.

**Status: the fixes described below have now been written** (`src/getLCAfromGAM.h`, `src/MCMC.h`, a note in `src/soibean.cpp`). The old lines are kept as `// OLD:` comments under explanatory headers, and
`analyse_GAM` also got the alignment-indexing fixes from TrailMix (edit-to-mapping lookup, read-order slicing, per-base quality/damage positions). They have only been syntax-checked, not built or run.

Original write-up (kept as written before the fix): nothing in the code had been touched. The goal is to explain,
in plain terms, what's mathematically wrong with how soibean currently scores a
read against a candidate species/haplotype when damage is involved, and what it
would take to bring soibean's code up to the same standard as the TrailMix fix.

## The recap: what the model is supposed to do

Like you said yourself: a C in your fragment that lines up with a T in the
reference path can happen for two reasons —

1. **A real mutation.** Rare, and it doesn't care where in the fragment it
   happens.
2. **Damage.** Much more likely in an ancient sample, and very much depends
   on *where* in the fragment the base sits (damage piles up near the ends of
   the fragment).

A sound model has to walk through this as a little chain of three events,
one after another:

1. **Mutation** — starting from the reference base, how likely is it that
   evolution actually changed it to some other base, given how long ago the
   two sequences split (the branch length)?
2. **Damage** — starting from *whatever base mutation left behind*, how
   likely is it that damage further changed it (C→T or G→A, depending on
   position in the fragment)?
3. **Sequencing error** — starting from *whatever base damage left behind*,
   how likely is it that the sequencer actually read out the base we see in
   the data?

Because we never directly observe the in-between states (the base right
after mutation, the base right after damage — we only see the final base in
the read), the correct way to score a read is to add up the probability of
*every possible path* through these three steps that ends up producing the
base we actually observed. That's the "Figure 2" model, and it's what was
fixed for TrailMix.

## What soibean does instead

Soibean is a separate, older piece of code that happens to solve a very
similar problem (how well does this read support this branch of the tree).
It was never updated to this fix because it has its *own, separate* copy of
the two pieces of code that matter here — not the ones that got fixed for
TrailMix. Think of it as a second implementation that was written around
the same time, then drifted apart over the years. So nothing TrailMix got
fixed automatically carried over.

When you read through soibean's version carefully, there isn't one problem —
there are three, stacked on top of each other.

### Problem 1: the damage/sequencing-error number is computed correctly... and then silently thrown away

Soibean's code *does* go through the trouble of combining the sequencing
quality score with the damage profile (your `--deam5p`/`--deam3p` settings)
into a number that's supposed to represent "how likely is it that damage and
then sequencing error turned the true base into the base we see". That part
is built correctly, base by base, position by position.

But then, when it's time to actually use that number for the base that was
*really* observed, the code accidentally treats a match and a mismatch
exactly the same way. It was clearly meant to say "if the candidate base
matches what we saw, do X; otherwise, do Y" — but both sides of that
decision end up doing the exact same thing. Mathematically, because the four
possible outcomes for a base always add up to 1 (something has to happen),
summing over "any of the four, I don't care which" always comes out to
`log(1) = 0`. Zero contribution. Every time. Whether the read matches or
completely disagrees with where it should sit.

In other words: soibean currently computes the damage information
correctly, and then discards it before it can influence the result. This is
the exact same category of mistake that was just found and fixed in
TrailMix (a copy/paste of an "if/else" that never actually differed) — it's
just living in a different, separate copy of the code.

### Problem 2: in place of real damage, it plugs in a number that has nothing to do with damage

Since the real damage term above contributes nothing, what's actually
driving the result in practice is a second term that was added alongside
it, meant to represent "some chance of error, so we never say something is
*impossible*". This term is computed using the classic error-rate formula —
mostly-right-answer gets `(1 - error rate)`, each wrong answer gets
`(error rate) / 3`.

The problem is what gets plugged in as the "error rate". It isn't the read's
quality score. It isn't the sample's damage profile. It's the length of the
*shortest branch anywhere in the entire reference tree* (or a flat 0.01 if
that number doesn't make sense). That's a property of the tree's geometry —
how finely resolved the reference tree happens to be — and it has nothing
to do with how damaged this particular read is, or where in the fragment
this particular base sits. It looks like an error-rate term because it's
wearing the right formula, but the number feeding it is borrowed from a
completely unrelated part of the program.

### Problem 3: no real "could have happened to any of the four bases" step at all

Even if Problem 1 and Problem 2 weren't there, soibean's damage number
(before it gets thrown away) is built assuming the base damage started from
*is* the reference base — it doesn't account for the fact that mutation
might have already changed the base before damage even got a chance to act
on it. That's the original TrailMix bug, the one you asked about at the very
start of this: damage needs to apply to *whatever mutation left behind*, not
blindly to the reference. Soibean's separate copy of the code has this same
assumption baked in, because it was written before this was worked out
properly.

### Putting it together

Right now, for a damaged ancient sample, soibean's placement decisions are
being driven almost entirely by:

- the ordinary mutation signal (computed correctly), plus
- a flat, made-up noise term borrowed from tree geometry,

with the actual damage signal (the part you specifically care about,
since it's what should make an ancient sample's C→T pattern near fragment
ends *informative* rather than *confusing*) contributing essentially
nothing. That's the opposite of what you'd want for degraded, damaged
material — exactly the regime where this matters most.

## What it would take to actually fix it

1. **Give soibean's own precompute step the same upgrade TrailMix's got.**
   Soibean has its own separate version of the "read every alignment and
   work out per-base probabilities" step (a different file from the one
   used by TrailMix). It needs the same addition: for each of the four
   possible "just-after-mutation" bases, precompute the chained probability
   that damage-then-sequencing-error turns *that* base into the base we
   actually observed. That gives four numbers instead of one, and — crucially —
   those four numbers need to actually depend on which of the four bases you
   start from, so a later step can pick the right one.

2. **Fix soibean's own version of the final combine-everything step**
   (also a separate copy from TrailMix's, used only by soibean) to do the
   proper "add up every possible path" calculation: take the mutation
   probability for each of the four possible intermediate bases, multiply
   each one by the matching damage-and-error number from step 1, and sum
   the four results. This replaces both of the broken pieces at once — the
   damage term that silently evaluates to zero, and the tree-branch-length
   number standing in for an error rate.

3. **Remove the tree-branch-length ("shortest branch in the tree") term
   from this calculation entirely.** It doesn't belong in a per-base
   likelihood at all — real sequencing error is already fully accounted for
   inside step 1's damage-and-error number once it's chained in properly.
   (This number is only used for this one purpose in soibean's code, so
   removing it from here doesn't appear to break anything else — but that's
   worth double-checking when the change is actually made.)

4. **Decide what "ancient vs. modern" means for soibean.** TrailMix has to
   juggle this per source, because it can mix an ancient sample with a
   modern contaminant in the same run. Soibean classifies one sample at a
   time and already has a single damage profile loaded from
   `--deam5p`/`--deam3p` (defaulting to "no damage" if you don't pass
   those). So soibean's version is actually simpler than TrailMix's — it
   just always uses its one loaded profile, with no per-source toggle
   needed.

## One more thing worth knowing

Because soibean keeps its own separate copy of this machinery rather than
sharing TrailMix's, fixing TrailMix did **not** fix soibean, and fixing
soibean later will not need to touch any of the TrailMix code. They're
independent fixes to independent code paths that happen to implement the
same idea.

Once this is fixed, it should be validated the same way the TrailMix fix
was: run it on a sample with a known, damaged, ancient origin and a
sample with little to no damage, and check that (a) the result doesn't get
worse for the non-damaged case, and (b) the posterior behaves better (less
noisy/more stable across MCMC chains) for the damaged case — the same kind
of before/after check that showed TrailMix's fix was actually working,
rather than just assuming the math is right because it reads correctly on
paper.
