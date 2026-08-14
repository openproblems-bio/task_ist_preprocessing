#!/usr/bin/env python3
"""Inspect a container image's build provenance from a registry without pulling it.

Usage: ghcr_inspect.py <registry> <image-repo-path> <tag> <expected-git-sha>

Prints key/value lines consumed by check_component.sh:
  STATUS ok|missing
  DETAIL <code>            (only when missing)
  REVISION <git-sha>       org.opencontainers.image.revision label
  CREATED <iso8601>
  CREATED_EPOCH <int>      CREATED as unix seconds (empty if unparseable)
  BASE_CREATED <iso8601>   created time of the first (base) layer
  LAYERS <n>
  MATCH yes|no             REVISION == expected-git-sha
"""
import json
import re
import subprocess
import sys
from datetime import datetime, timezone


def to_epoch(iso):
    """ISO8601 (with optional fractional seconds and trailing Z) -> unix seconds.
    Fractional seconds are dropped; registry timestamps are UTC. '' if unparseable."""
    m = re.match(r"(\d{4})-(\d{2})-(\d{2})T(\d{2}):(\d{2}):(\d{2})", iso or "")
    if not m:
        return ""
    y, mo, d, h, mi, s = (int(x) for x in m.groups())
    return str(int(datetime(y, mo, d, h, mi, s, tzinfo=timezone.utc).timestamp()))

ACCEPT = (
    "application/vnd.oci.image.index.v1+json,"
    "application/vnd.oci.image.manifest.v1+json,"
    "application/vnd.docker.distribution.manifest.list.v2+json,"
    "application/vnd.docker.distribution.manifest.v2+json"
)


def curl(url, accept=None, token=None, follow=False):
    cmd = ["curl", "-s"]
    if follow:
        cmd.append("-L")
    if token:
        cmd += ["-H", f"Authorization: Bearer {token}"]
    if accept:
        cmd += ["-H", f"Accept: {accept}"]
    cmd.append(url)
    return subprocess.run(cmd, capture_output=True, text=True).stdout


def jload(s):
    try:
        return json.loads(s)
    except Exception:
        return {}


def main():
    if len(sys.argv) < 5:
        print("STATUS missing\nDETAIL bad-args")
        return
    reg, img, tag, main_sha = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
    base = "https://ghcr.io" if reg.strip() == "ghcr.io" else f"https://{reg}"

    tok = jload(curl(f"{base}/token?scope=repository:{img}:pull&service={reg}")).get("token", "")
    man = jload(curl(f"{base}/v2/{img}/manifests/{tag}", ACCEPT, tok))
    if not man or man.get("errors"):
        code = man.get("errors", [{}])[0].get("code") if man else "NO_RESPONSE"
        print(f"STATUS missing\nDETAIL {code}")
        return

    if "manifests" in man:  # multi-arch index -> prefer linux/amd64
        child = None
        for m in man["manifests"]:
            p = m.get("platform", {})
            if p.get("os") == "linux" and p.get("architecture") == "amd64":
                child = m["digest"]
                break
        child = child or man["manifests"][0]["digest"]
        man = jload(curl(f"{base}/v2/{img}/manifests/{child}", ACCEPT, tok))

    cfg_digest = (man.get("config") or {}).get("digest", "")
    if not cfg_digest:
        print("STATUS missing\nDETAIL no-config-digest")
        return

    cfg = jload(curl(f"{base}/v2/{img}/blobs/{cfg_digest}", token=tok, follow=True))
    labels = (cfg.get("config") or {}).get("Labels") or {}
    rev = labels.get("org.opencontainers.image.revision", "")
    if not rev:  # fall back to scanning build history for the LABEL line
        for h in cfg.get("history", []):
            cb = h.get("created_by", "")
            if "org.opencontainers.image.revision=" in cb:
                rev = cb.split("org.opencontainers.image.revision=", 1)[1].split()[0].strip('"')
                break

    hist = cfg.get("history", [])
    base_created = hist[0].get("created", "") if hist else ""
    match = "yes" if (rev and main_sha and rev == main_sha) else "no"
    created = cfg.get("created", "")
    print("STATUS ok")
    print(f"REVISION {rev}")
    print(f"CREATED {created}")
    print(f"CREATED_EPOCH {to_epoch(created)}")
    print(f"BASE_CREATED {base_created}")
    print(f"LAYERS {len(man.get('layers', []))}")
    print(f"MATCH {match}")


if __name__ == "__main__":
    main()
