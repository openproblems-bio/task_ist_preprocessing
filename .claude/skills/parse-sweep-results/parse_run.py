#!/usr/bin/env python3
"""Analyse a Seqera run's task list into a per-dataset x pipeline-combo end-to-end
success/failure table for a task_ist_preprocessing benchmark run.

Works for ANY run, not just single-method parameter sweeps:
  * one method swept over parameters   -> combos = param variants
  * several methods at one stage        -> combos = the methods (param "default")
  * mixed / multiple varying stages     -> combos = the tuple of varying tokens

Reads the task JSON (list of {status, tag, process, exit}) fetched by parse_run.sh.
Task tags encode the pipeline config as
    <dataset>/<stage1_token>/<stage2_token>/.../metric_<name>
controls appear as <dataset>/control_<name>/metric_<name>.

"succeeded in the end" == the combo has a COMPLETED terminal-metric task.

Two Seqera facts this handles:
  * task TAGS are TRUNCATED, so the trailing "/metric_<name>" is chopped for long combos
    -> metric tasks are detected by PROCESS name (never truncated), discovered from the
       short tags that still show a metric_ segment.
  * the swept/ varying token sits near the FRONT of the tag, which survives truncation.
"""
import argparse, json, os, sys
from collections import defaultdict

OK = {"COMPLETED", "CACHED"}
# stage-token prefix -> results-tree folder
STAGE_FOLDER = {
    "segm": "segmentation", "ass": "transcript_assignment",
    "cta": "cell_type_annotation", "corr": "expression_correction",
    "aggr": "count_aggregation", "qc": "qc_filter", "cell": "volume",
    "norm": "normalization", "gene": "gene_efficiency",
}


def chain_of(tag):
    return tag.split("/")[1:]


def combo_chain(tag):
    """tokens after the dataset, with any metric_ segment removed."""
    return [s for s in chain_of(tag) if not s.startswith("metric_")]


def is_control(ch):
    return any(s.startswith("control_") for s in ch)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tasks", required=True)
    ap.add_argument("--id", required=True)
    ap.add_argument("--base", required=True)
    ap.add_argument("--stage", default=None)
    ap.add_argument("--method", default=None)
    ap.add_argument("--overall", default="?")
    ap.add_argument("--run-name", default="?")
    args = ap.parse_args()

    tasks = json.load(open(args.tasks))
    if not tasks:
        print("No tasks for run", args.id, "(still queued, or wrong id/workspace)."); sys.exit(1)

    # ---- metric tasks by PROCESS (truncation-proof) ----
    metric_procs = {t.get("process", "") for t in tasks
                    if any(s.startswith("metric_") for s in chain_of(t["tag"]))}

    def variance(drop_last):
        """distinct-token sets per chain index over non-metric, non-control tasks."""
        idx = defaultdict(set)
        for t in tasks:
            if t.get("process", "") in metric_procs:
                continue
            ch = combo_chain(t["tag"])
            if is_control(ch):
                continue
            toks = ch[:-1] if (drop_last and len(ch) > 1) else ch
            for i, tok in enumerate(toks):
                idx[i].add(tok)
        return idx

    # Drop each chain's LAST token first (it may be a truncated tail); if that finds no
    # variation (e.g. everything failed at the first stage), retry without dropping.
    idx_vals = variance(drop_last=True)
    V = sorted(i for i, s in idx_vals.items() if len(s) > 1)
    if not V:
        idx_vals = variance(drop_last=False)
        V = sorted(i for i, s in idx_vals.items() if len(s) > 1)
    multi = len(V) > 1

    def key_of(ch):
        if not V:
            return ("(all-default)",)
        return tuple(ch[i] if i < len(ch) else "?" for i in V)

    # ---- aggregate per (dataset, combo) ----
    attempted = defaultdict(lambda: {"metric_ok": False, "metric_any": False})
    controls = defaultdict(lambda: {"metric_ok": False, "metric_any": False})
    failures = []                    # raw non-OK: (ds, chain, step, status, exit, taskId)
    failinfo = defaultdict(list)     # (ds, combo-key) -> [(taskId, step, exit)] of its non-OK tasks
    for t in tasks:
        tag, st, proc = t["tag"], t.get("status", "?"), t.get("process", "")
        tid = t.get("taskId")
        ds = tag.split("/")[0]
        ch = combo_chain(tag)
        is_metric = proc in metric_procs
        step = proc.split(":")[-1].replace("_process", "")   # e.g. stardist, basic_qc_filter
        if st not in OK:
            failures.append((ds, "/".join(chain_of(tag)), step, st, t.get("exit"), tid))
        ctrl = next((s for s in ch if s.startswith("control_")), None)
        if ctrl:
            c = controls[(ds, ctrl)]
            if is_metric:
                c["metric_any"] = True; c["metric_ok"] |= st in OK
            continue
        k = key_of(ch)
        if "?" in k:
            continue   # auxiliary task (e.g. extract_uns_metadata) with no combo token
        a = attempted[(ds, k)]
        if is_metric:
            a["metric_any"] = True; a["metric_ok"] |= st in OK
        if st not in OK:
            failinfo[(ds, k)].append((tid, step, t.get("exit")))

    datasets = sorted({d for (d, _) in attempted} | {d for (d, _) in controls})
    combos = sorted({k for (_, k) in attempted})

    # ---- label combos + parse method/param/value (single varying axis) ----
    def stage_of(tok):
        return STAGE_FOLDER.get(tok.split("_")[0], tok.split("_")[0])

    single = (len(V) == 1)
    parsed = {}   # key-tuple -> (stage, method, param, value, label)
    if single:
        toks = [k[0] for k in combos]
        prefix = toks[0].split("_")[0] if toks else ""
        def rem(tok):
            return tok[len(prefix) + 1:] if tok.startswith(prefix + "_") else tok
        rems = [rem(t) for t in toks]
        # a remainder is a "base method" if no OTHER remainder is a proper prefix of it
        base = [r for r in rems if not any(o != r and r.startswith(o + "_") for o in rems)]
        def parse(tok):
            r = rem(tok)
            cand = [m for m in base if r == m or r.startswith(m + "_")]
            method = max(cand, key=len) if cand else r
            if r == method:
                return stage_of(tok), method, "default", "-"
            rest = r[len(method) + 1:]
            p, v = rest.rsplit("_", 1) if "_" in rest else (rest, "-")
            return stage_of(tok), method, p, v
        for k in combos:
            s, m, p, v = parse(k[0])
            label = m if p == "default" else f"{m} {p}={v}"
            parsed[k] = (s, m, p, v, label)
        methods_with_params = sorted({parsed[k][1] for k in combos if parsed[k][2] != "default"})
        stage_folder = args.stage or (parsed[combos[0]][0] if combos else "misc")
        fname_method = args.method or (methods_with_params[0] if len(methods_with_params) == 1
                                       else stage_folder)
    else:
        for k in combos:
            label = " | ".join(k)
            parsed[k] = ("misc", "-", "-", "-", label)
        # infer stage folder from the first varying index's tokens
        first_tok = combos[0][0] if combos and combos[0][0] not in ("(all-default)", "?") else ""
        stage_folder = args.stage or (stage_of(first_tok) if first_tok else "misc")
        fname_method = args.method or stage_folder

    # ---- write CSV ----
    out_dir = os.path.join(args.base, stage_folder)
    os.makedirs(out_dir, exist_ok=True)
    csv_path = os.path.join(out_dir, f"{fname_method}_run_{args.id}_results.csv")
    rows = []
    for (ds, k) in sorted(attempted):
        s, m, p, v, _ = parsed[k]
        a = attempted[(ds, k)]
        rows.append((ds, " | ".join(k), s, m, p, v,
                     "yes" if a["metric_any"] else "no",
                     "SUCCEEDED" if a["metric_ok"] else "FAILED"))
    for (ds, ctrl) in sorted(controls):
        c = controls[(ds, ctrl)]
        rows.append((ds, ctrl, "control", ctrl[len("control_"):], "control", "-",
                     "yes" if c["metric_any"] else "no",
                     "SUCCEEDED" if c["metric_ok"] else "FAILED"))
    with open(csv_path, "w") as f:
        f.write("dataset,combo,stage,method,param,value,reached_metric,end_to_end_status\n")
        for r in rows:
            f.write(",".join(str(x) for x in r) + "\n")

    # ---- classify failures: method-implementation error vs param-combo error ----
    # Heuristics (see SKILL.md "Classifying failures"):
    #   * where it dies: FAILED in the combo's OWN method step -> suspect the method/impl;
    #     FAILED in a DOWNSTREAM step -> the method ran, a param produced degenerate output.
    #   * failure ratio: if ~all of a method's variants fail on a dataset, it's method/dataset
    #     level (not one param); if only this variant fails, it's param-specific.
    ds_method_total = defaultdict(int); ds_method_fail = defaultdict(int)
    for (ds, k) in attempted:
        m = parsed[k][1]
        ds_method_total[(ds, m)] += 1
        if not attempted[(ds, k)]["metric_ok"]:
            ds_method_fail[(ds, m)] += 1

    def classify(ds, k):
        m = parsed[k][1]
        info = failinfo.get((ds, k), [])
        steps = sorted({s for (_, s, _) in info})
        tid = info[0][0] if info else None
        ex = info[0][2] if info else "?"
        own = m in steps                                   # failed in its own method step
        tot, fl = ds_method_total[(ds, m)], ds_method_fail[(ds, m)]
        if not steps:
            cls = "no-failing-task (metric missing; upstream/aux or truncation)"
        elif not own:
            cls = f"PARAM-COMBO -> degenerate output; crashes downstream in {'/'.join(steps)}"
        elif tot > 1 and fl >= max(2, tot - 1):
            cls = (f"METHOD/DATASET-level: {fl}/{tot} of {m}'s combos fail here "
                   f"-> impl bug or data-incompat (not one param)")
        else:
            cls = f"PARAM-SPECIFIC at method step ({'/'.join(steps)}) -> bad param or fragile handling"
        return steps, tid, ex, cls.replace(",", ";")

    fail_combos = [(ds, k) for (ds, k) in sorted(attempted) if not attempted[(ds, k)]["metric_ok"]]

    # ---- TRIAGE: flag for troubleshooting ONLY methods that fail on their DEFAULT variant ----
    # A default-variant failure = the method is broken out-of-the-box on that data (a real
    # implementation / robustness issue worth troubleshooting). A failure only on a NON-default
    # param = the param choice / a fragile edge — recorded, but NOT flagged. (In a multi-method
    # run every method's sole variant IS its default, so any such failure flags that method.)
    method_default_fail = set()
    for (ds, k) in attempted:
        _, m, p, _, _ = parsed[k]
        if p == "default" and not attempted[(ds, k)]["metric_ok"]:
            method_default_fail.add(m)
    methods_with_fail = {parsed[k][1] for (ds, k) in fail_combos}
    flagged_methods = sorted(method_default_fail)
    param_only_methods = sorted(methods_with_fail - method_default_fail)

    fail_csv = os.path.join(out_dir, f"{fname_method}_run_{args.id}_failures.csv")
    ft_path = os.path.join(out_dir, f".failed_tasks_{args.id}.tsv")   # hand-off for parse_run.sh --diagnose
    if fail_combos:
        with open(fail_csv, "w") as f, open(ft_path, "w") as ft:
            f.write("dataset,combo,method,param,value,is_default,failing_step,exit,task_id,troubleshoot,likely_class\n")
            for (ds, k) in fail_combos:
                s, m, p, v, _ = parsed[k]
                steps, tid, ex, cls = classify(ds, k)
                flag_method = m in method_default_fail          # method-level -> CSV `troubleshoot`
                flag_diag = (p == "default")                    # this row is a default-variant failure
                f.write(",".join(str(x) for x in
                        (ds, " | ".join(k), m, p, v, "yes" if p == "default" else "no",
                         "/".join(steps) or "?", ex, tid, "yes" if flag_method else "no", cls)) + "\n")
                if tid is not None:   # --diagnose default fetches only default-variant failures (flag_diag)
                    ft.write(f"{tid}\t{ds}\t{' | '.join(k)}\t{'/'.join(steps)}\t{ex}\t{1 if flag_diag else 0}\n")

    # ---- report ----
    n_ok = sum(1 for r in rows if r[7] == "SUCCEEDED")
    print(f"Run {args.id} ({args.run_name})  overall={args.overall}")
    axis = "none" if not V else (f"single (chain index {V[0]})" if single else f"MULTIPLE {V}")
    print(f"datasets: {len(datasets)}  |  varying axis: {axis}  |  combos/dataset: {len(combos)}")
    if single:
        print(f"methods: {sorted({parsed[k][1] for k in combos})}")
    print(f"tasks: {len(tasks)}  |  non-OK tasks: {len(failures)}")
    print(f"rows (incl. controls): {len(rows)}  ->  SUCCEEDED {n_ok}  FAILED {len(rows) - n_ok}")
    print(f"CSV: {csv_path}\n")

    # matrix: combos (rows) x datasets (cols)
    labels = {k: parsed[k][4] for k in combos}
    ctrl_ks = sorted({c for (_, c) in controls})
    width = max([len(labels[k]) for k in combos] + [len("control: " + c[8:]) for c in ctrl_ks] + [12])
    ds_short = [d.replace("_test_combined", "") for d in datasets]
    print("  " + "combo".ljust(width) + "   " + "  ".join(s[:22].center(22) for s in ds_short))
    def cell(present, ok):
        return ("PASS" if ok else "FAIL") if present else " -  "
    for k in combos:
        line = "  " + labels[k].ljust(width)
        for ds in datasets:
            a = attempted.get((ds, k))
            line += "   " + cell(a is not None, a and a["metric_ok"]).center(22)
        print(line)
    for c in ctrl_ks:
        line = "  " + ("control: " + c[8:]).ljust(width)
        for ds in datasets:
            cc = controls.get((ds, c))
            line += "   " + cell(cc is not None, cc and cc["metric_ok"]).center(22)
        print(line)

    if fail_combos:
        print(f"\nFAILURE DIAGNOSIS ({len(fail_combos)} failed combos; {len(failures)} non-OK task records):")
        if flagged_methods:
            print(f"  >> TROUBLESHOOT (fail on DEFAULT variant): {', '.join(flagged_methods)}")
        if param_only_methods:
            print(f"  .. param-only failures (NOT flagged): {', '.join(param_only_methods)}")
        print()
        for (ds, k) in fail_combos:
            _, m, _, _, lab = parsed[k]
            steps, tid, ex, cls = classify(ds, k)
            mark = "[FLAG]" if m in method_default_fail else "  -   "
            print(f"  {mark} {ds.replace('_test_combined','')}  |  {lab}")
            print(f"         step={'/'.join(steps) or '?'}  exit={ex}  taskId={tid}  ->  {cls}")
        print(f"\n  Failures CSV: {fail_csv}")
        print("  Retry note: exit-1 combos get 2 attempts, OOM(137/139) up to 3, then 'ignore' —")
        print("  so a FAILED combo failed on ALL attempts. Re-run with --diagnose for tracebacks")
        print("  (flagged/fail-on-default only; add --all-failures to include param-only ones).")
        print("  Tip: an impl-error verdict may just be an UNDEPLOYED fix — cross-check with")
        print("  check-component before blaming the code.")
    if multi:
        print("\nNOTE: multiple varying axes detected — combos are labelled by the full "
              "tuple of varying tokens; deep tokens may be affected by tag truncation.")
    print(f"\nOUT_DIR={out_dir}")
    if fail_combos:
        print(f"FAILURES_CSV={fail_csv}")
        print(f"FAILED_TASKS={ft_path}")
    print(f"WROTE_CSV={csv_path}")


if __name__ == "__main__":
    main()
