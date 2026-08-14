#!/usr/bin/env python3
"""Validate a run_benchmark `method_params` sweep block.

Usage:
    validate_sweep.py <method_params.yaml> [component_name]

The file is the one referenced by `method_parameters_yaml` — a mapping shaped as:

    parameters:
      <component>:
        default: { arg: val, ... }
        sweep:   { arg: [v1, v2, ...], ... }

For each method it reports the default args, each swept parameter's value count, and the
TOTAL number of variants the workflow will expand — which is `1 (default) + Σ len(sweep list)`
because the benchmark varies ONE parameter at a time (a "star" around the default, not a
grid). It warns when a sweep list repeats the default value (a redundant duplicate variant)
and when the sweep is heavy for a quick test run.

Exit status is non-zero if the file has no usable `parameters:` block, so it can gate a script.
"""
import sys
import yaml


def check_method(comp, spec):
    spec = spec or {}
    default = spec.get("default", {}) or {}
    sweep = spec.get("sweep", {}) or {}
    print(f"\n== {comp} ==")
    print(f"  default: {default or '(none)'}")
    n_variants = 1  # the default variant
    for par, vals in sweep.items():
        if not isinstance(vals, list):
            vals = [vals]
        n_variants += len(vals)
        # bool must be checked before the general == (True == 1 in Python)
        dups = [
            v for v in vals
            if par in default
            and type(v) is type(default[par])
            and v == default[par]
        ]
        warn = ""
        if dups:
            warn = f"   ** repeats default ({default[par]}) -> redundant variant, drop it"
        elif isinstance(default.get(par), bool):
            warn = "   (boolean: one non-default value is all there is)"
        print(f"  sweep {par}: {len(vals)} value(s) {vals}{warn}")
    print(f"  -> variants: 1 default + {n_variants - 1} sweep = {n_variants}")
    return n_variants


def main(argv):
    if len(argv) < 2:
        print("usage: validate_sweep.py <method_params.yaml> [component]", file=sys.stderr)
        return 2
    with open(argv[1]) as fh:
        doc = yaml.safe_load(fh)
    only = argv[2] if len(argv) > 2 else None

    params = doc.get("parameters", doc) if isinstance(doc, dict) else None
    if not isinstance(params, dict) or not params:
        print("no usable `parameters:` block found", file=sys.stderr)
        return 1

    total = 0
    seen = False
    for comp, spec in params.items():
        if only and comp != only:
            continue
        seen = True
        total += check_method(comp, spec)
    if only and not seen:
        print(f"component '{only}' not found in file", file=sys.stderr)
        return 1

    print(f"\nTOTAL variants across method(s): {total}")
    if total > 12:
        print("NOTE: >12 variants -> heavy for a quick test run; consider trimming lists.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
