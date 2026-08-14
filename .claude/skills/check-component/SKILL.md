---
name: check-component
description: >-
  Verify a Viash component's change is fully deployed before (re)running a benchmark:
  committed to origin/main, regenerated in build/main, and rebuilt as a fresh container
  on ghcr (build_main tag). Use when a component edit doesn't seem to take effect on the
  cluster / in a Nextflow or Seqera run, when a run errors with symptoms of an old image
  (missing deps, wrong base image), or to confirm an image is not stale before launching.
---

# check-component

A component change only reaches a benchmark run after **three** independent things
happen. A failure at any stage means the cluster keeps using the **old container**, so
the fix silently has no effect:

1. **Source** — the change is committed and on `origin/main`.
2. **Codegen** — CI's `Build` workflow regenerates the `build/main` deploy branch
   (the inlined `target/nextflow/.../main.nf` that run scripts launch via
   `--revision build/main`).
3. **Container** — CI builds the component's Docker image and **pushes it to `ghcr.io`**
   at the `build_main` tag. This is the step that most often lags or fails (heavy images,
   source compiles, runner disk/time limits). Codegen redeploying does **not** imply the
   container was rebuilt — they are separate steps.

The decisive signal is the image's `org.opencontainers.image.revision` label: it records
the git commit the container was built from. Comparing it to `origin/main` HEAD tells you
whether the pushed image is current — read straight from the registry, without pulling the
(multi-GB) image.

## How to run

```bash
.claude/skills/check-component/check_component.sh <component>
```

Accepts a bare component name (`segger`), a `namespace/component` path
(`methods_transcript_assignment/segger`), or a path to `config.vsh.yaml`. Optional 2nd arg
overrides the image tag (default `build_main`; for a `build/<branch>` deploy the tag is
`build_<branch>`).

The script prints the three checks and a `VERDICT`. Report the verdict to the user and, when
something is stale, which stage failed and what to do about it.

## Interpreting results

- **All three OK** → a run with `--pull-latest` will use the current build. Safe to launch.
- **[1] source not on main** → the fix is uncommitted or on another branch (e.g. `fixes`).
  Commit/merge to `main` first; the push to `main` is what triggers CI's `Build`.
- **[2] build/main behind** → the `Build` workflow hasn't redeployed the latest `main` yet.
  It usually runs within ~1 min of the merge; if not, check the Actions "Build" run.
- **[3] container stale** (revision label ≠ `origin/main` HEAD) → the container push for the
  latest commit hasn't succeeded. Either the `Build` run is still compiling the image, or the
  image build **failed** (common for heavy GPU/compile images). Check the `Build` run for this
  component; if it failed, rebuild + push the image manually:

  ```bash
  # from a linux/amd64 host with docker + `docker login ghcr.io`
  viash ns build -q <component> --setup push
  ```

  (macOS/Apple-Silicon can't practically build amd64 CUDA images — use a Linux build host or
  a build pod on the cluster.)

- **[3] container `created-vs-commit: STALE`** (revision label **matches** HEAD, but the image
  was **built before** the commit it is labelled with) → a **false-OK trap** the revision label
  alone misses. The `revision` label was re-stamped onto old layers without a real rebuild —
  Docker layer-cache reuse, a re-pushed/re-tagged manifest, or a heavy build that failed and
  left the previous image in place. The deployed image does **not** contain the change even
  though it claims the right SHA. Fix is the same as a plain stale container, but you must
  **bust the cache** so the changed layers actually re-execute (build on a fresh builder, e.g. a
  clean Docker-in-Docker pod, rather than one with warm cache). Tell: `image built` timestamp is
  earlier than the revision commit's time.

## Notes

- Requires network access to the registry and `git fetch` access to `origin`. The script does
  a `git fetch origin` itself.
- The container check queries `ghcr.io` with an anonymous pull token, which works for public
  packages. For a private package or a non-ghcr registry the container stage reports
  `NOT FOUND`; verify via `docker manifest inspect` with credentials instead.
- Project identity (organization, name, registry) is read from `_viash.yaml`; the image repo
  path is `<registry>/<organization>/<name>/<namespace>/<component>`.
