---
name: document-method-troubleshooting
description: >-
  Capture what was learned while troubleshooting a Viash method/metric component into a
  durable two-layer record: an in-repo NOTES.md (the authoritative how-it-works / why-the-
  setup-is-this-way / where-it-breaks reference, committed with the code) plus a memory file
  (a short pointer + a dated iteration log). Use after debugging a component (build/deps/
  container/CLI/runtime issue), when the user says "document this", "write it up", "remember
  the segger saga", "add to the notes", or asks to keep a running log of troubleshooting
  iterations. Invoke again after EACH further iteration to append the log.
---

# document-method-troubleshooting

Troubleshooting any component (segger, baysor, an atera loader, a metric, …) produces
hard-won knowledge — why the Docker base or a version pin is load-bearing, a subtle
logic/format/output-contract bug and its fix, what dead ends were ruled out, what still isn't
confirmed. That knowledge is worthless if it evaporates when the session ends or lives only
in a commit message nobody re-reads. This skill writes it down in **two layers that do not
duplicate each other**:

| Layer | Where | Lives with | Contains |
|-------|-------|-----------|----------|
| **NOTES.md** | `src/<namespace>/<component>/NOTES.md` | the code, committed to git | the durable detail: how the component works step-by-step, why the setup/pins are the way they are, arguments, risk points/gotchas, wiring |
| **Memory** | auto-memory dir, `<component>-<stage>.md` (type `project`) | the session (recalled automatically) | a short pointer to NOTES.md + the current one-paragraph state + a **dated iteration log** |

**The rule that keeps them from drifting: deep detail goes in NOTES.md; the memory is a
pointer plus a running log.** If you're tempted to paste the Docker recipe rationale into the
memory, it belongs in NOTES.md. If you're tempted to put "tried X on 2026-07-23, still fails"
into NOTES.md, it belongs in the memory iteration log.

## When invoked

1. **Identify the component.** From the conversation or the user's words, resolve it to a
   `src/<namespace>/<component>/` dir (e.g. `methods_transcript_assignment/segger`). Confirm
   the dir exists.
2. **Check what already exists** — a `NOTES.md` in that dir, and a memory file for it
   (search the auto-memory dir + `MEMORY.md` index). This decides create-vs-append below.
3. **Do the smaller of two jobs:**
   - **First time / understanding changed** → write or revise `NOTES.md`, then create/refresh
     the memory pointer.
   - **Just another iteration** (tried something, got a result) → append one line to the
     memory's iteration log; touch NOTES.md only if the *understanding* changed (a new root
     cause, a new pin, a ruled-out dead end), not merely because an attempt was made.

## NOTES.md — the in-repo reference

Committed alongside the code so it survives context loss and helps the next person (human or
Claude). Write it as prose a maintainer would want before touching the component. **Scale it
to the component** — a method that needed one small fix gets a few honest paragraphs; a
fought-over one (many pins, dead ends, external CLI) gets the full treatment. Use the
sections that apply and drop the rest — don't pad a simple component to fit a big template.
Candidate sections:

- **What this component is** — its stage/API (`src/api/comp_*` for methods; `comp_metric_*`
  for metrics), what's unusual about it (GPU-only? R vs Python? wraps an external CLI?
  adapter?), links (docs/repo/DOI).
- **script.py / script.R — step by step** — the real control flow, with `file:line` anchors,
  and the *why* behind non-obvious choices.
- **Arguments** — table of args, defaults, and what they map to upstream.
- **Setup / Docker (often the highest-value section when it applies)** — why the base image,
  why each version pin, any ORDER constraints, and the history: what was tried before and why
  it failed (e.g. "NGC base → UCX collision → SIGABRT, unfixable by pip pins → rebased on
  RAPIDS"). This is what stops someone re-breaking it. If the component just merges a base
  setup and worked first try, say that in a line — no drama to record.
- **Wiring** — where it's registered (for a method: workflow config dep list + default
  methods string, `main.nf` fan-out; metrics wire via the metrics workflow), the run/test
  script.
- **Risk points / gotchas** — commit-pinned deps and why, output-contract assumptions, CI
  coverage gaps, anything fragile.

Keep `file:line` references honest — verify them against current code before writing them as
fact (they rot). State what is **validated** vs **not yet confirmed** explicitly.

## Memory — the pointer + iteration log

Follow the repo's memory conventions (see the memory instructions in the system context):
one file, frontmatter with `type: project`, body, and a one-line pointer added to `MEMORY.md`.
Structure the body as:

- **One paragraph of current state** — what the component is, and the non-obvious facts
  (load-bearing pins, commit pin, "not yet run end-to-end on GPU", CI gaps). Convert relative
  dates to absolute.
- **A line that says NOTES.md is the authoritative detail** and this memory is the index/log,
  so future-you reads NOTES.md for depth and doesn't try to keep both fully detailed.
- **`**How to apply:**`** — note that the log is appended after each iteration.
- **Links** — `[[other-memory]]` cross-refs to related components/lessons.
- **`## Iteration log`** — reverse-or-forward-chronological dated lines. **This is the part
  that grows.**

### Iteration log entry format

One line per troubleshooting round: `- YYYY-MM-DD — <what was tried> → <result / current
blocker>`. Be concrete and include dead ends (a ruled-out cause is as valuable as a fix):

```
- 2026-07-23 — memory created; state captured. Builds/imports validated on CPU (rmm mocked);
  awaiting first confirmed GPU end-to-end run.
- 2026-07-24 — ran run_test_segger_nebius.sh on Nebius GPU → failed at <stage> with <error>;
  suspect <hypothesis>. Ruled out stale image (check-component: revision current).
```

When an iteration changes a durable fact (new pin, confirmed root cause, "now runs
end-to-end"), **also** update the state paragraph and, if it's about how the component works
or breaks, NOTES.md — then log that you did.

## Output

Tell the user, briefly: which files you wrote/appended, and (if this was an iteration) the
one-line log entry you added. Don't reprint the whole NOTES.md.
