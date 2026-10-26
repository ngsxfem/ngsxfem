"""Plan and prepare release candidates / releases for .github/workflows/pypi.yml.

    release_plan.py plan      decide what the CI run does (writes GITHUB_OUTPUT)
    release_plan.py prepare   create/update the release-candidate branch (bumps the NGSolve pin)

Flow:
  1. The ngsolve-release-watcher dispatches pypi.yml on master with the new
     NGSolve version.  If master is not pinned to it yet, `prepare` creates the
     branch `release-candidate` (master + one commit that only bumps the pin).
  2. Builds of that branch are published as release candidates 2.1.<ngsolve>rcN.
     A pull request release-candidate -> master collects them.
  3. Merging it changes the pin on master: that push is built and published as
     the release 2.1.<ngsolve> (and tagged).

Everything is decided from the pin in pyproject.toml, the tags and PyPI, so a
failed run can simply be repeated.
"""
import os
import re
import subprocess
import sys
import tomllib
from urllib.error import HTTPError
from urllib.request import urlopen

from packaging.version import InvalidVersion, Version

RC_BRANCH = "release-candidate"
DRYRUN_BRANCH = "release-candidate-dryrun"
MASTER = "master"


def git(*args, check=True):
    result = subprocess.run(["git", *args], capture_output=True, text=True)
    if check and result.returncode:
        raise RuntimeError(f"git {' '.join(args)}: {result.stderr.strip()}")
    return result.stdout.strip() if result.returncode == 0 else None


def ngsolve_pin(pyproject_text):
    """The exact NGSolve version `dependencies` is pinned to (or None)."""
    project = tomllib.loads(pyproject_text)
    for dep in project["project"].get("dependencies", []):
        match = re.fullmatch(r"\s*ngsolve\s*==\s*([^\s;]+)\s*", dep)
        if match:
            return match[1]


def package_name(pyproject_text):
    return tomllib.loads(pyproject_text)["project"]["name"]


def release_series(ngsolve_version, tags):
    """`2.1.2608` for NGSolve 6.2.2608: major.minor of the latest stable tag + NGSolve's micro."""
    stable = []
    for tag in tags:
        try:
            v = Version(tag.removeprefix("v"))
        except InvalidVersion:
            continue
        if not v.is_prerelease and not v.is_devrelease and not v.is_postrelease:
            stable.append(v)
    latest = max(stable)
    return f"{latest.major}.{latest.minor}.{Version(ngsolve_version).micro}"


def pypi_versions(name):
    import json
    try:
        with urlopen(f"https://pypi.org/pypi/{name}/json", timeout=30) as r:
            releases = json.load(r)["releases"]
    except HTTPError as error:
        if error.code != 404:
            raise
        return []
    return [Version(v) for v, files in releases.items()
            if any(not f.get("yanked", False) for f in files)]


def next_rc(series, published, tags):
    """Number of the next release candidate of `series` (PyPI and git tags)."""
    numbers = [0]
    for v in published:
        if v.base_version == series and v.pre and v.pre[0] == "rc":
            numbers.append(v.pre[1])
    for tag in tags:
        match = re.fullmatch(rf"v{re.escape(series)}rc(\d+)", tag)
        if match:
            numbers.append(int(match[1]))
    return max(numbers) + 1


def plan(event, ref_name, head_ref, input_ngsolve, dry_run, pin_now, pin_before,
         published, tags, tags_at_head):
    """Pure decision function.  Returns the dict written to GITHUB_OUTPUT."""
    out = {"mode": "none", "skip": "false", "prepare": "false", "publish": "false",
           "ngsolve_version": "", "version": "", "reason": ""}

    def done(reason, skip=False):
        out["reason"] = reason
        out["skip"] = str(skip).lower()
        if skip:
            out["prepare"] = out["publish"] = "false"
        return out

    if event == "pull_request":
        if head_ref in (RC_BRANCH, DRYRUN_BRANCH):
            return done("built by the push to the release-candidate branch", skip=True)
        return done("pull request: build and test only")
    if event == "push" and ref_name.startswith("v") and ref_name[1:2].isdigit():
        # manually pushed tag: build as tagged (version from the tag, pins from the commit)
        out.update(mode="tag", publish="true")
        return done(f"tag {ref_name}")

    on_rc_branch = ref_name in (RC_BRANCH, DRYRUN_BRANCH)
    dispatch = event == "workflow_dispatch" and bool(input_ngsolve)
    if dispatch and not on_rc_branch:
        ngsolve = input_ngsolve
        if pin_now != ngsolve:
            out["prepare"] = "true"
        mode = "rc" if pin_now != ngsolve else "release"
    elif on_rc_branch and event in ("push", "workflow_dispatch"):
        ngsolve, mode = pin_now, "rc"
    elif event == "push" and ref_name == MASTER:
        if pin_before is None or pin_before == pin_now:
            return done("master: NGSolve pin unchanged, build and test only")
        ngsolve, mode = pin_now, "release"
    else:
        return done("build and test only")

    series = release_series(ngsolve, tags)
    if any(v.base_version == series and not v.is_prerelease for v in published):
        return done(f"{series} is already released on PyPI", skip=True)
    out.update(mode=mode, ngsolve_version=ngsolve)
    if mode == "rc":
        rc_tags = [t for t in tags_at_head if re.fullmatch(rf"v{re.escape(series)}rc\d+", t)]
        if rc_tags and out["prepare"] == "false":
            return done(f"this commit was already built as {rc_tags[0]}", skip=True)
        out["version"] = f"{series}rc{next_rc(series, published, tags)}"
    else:
        out["version"] = series
    out["publish"] = str(not dry_run).lower()
    out["reason"] = f"{mode} {out['version']} for NGSolve {ngsolve}" + (" (dry run)" if dry_run else "")
    return out


def bump_pin(text, ngsolve):
    return re.sub(r'(?<=ngsolve==)[0-9][^"\',;\s]*', ngsolve, text)


def bump_notes(text, series, ngsolve):
    text = re.sub(r"^## \d+\.\d+\.\d+ \(unreleased\)", f"## {series} (unreleased)", text,
                  count=1, flags=re.M)
    return re.sub(r"^Requires NGSolve [^\s]+?\.$", f"Requires NGSolve {ngsolve}.", text,
                  count=1, flags=re.M)


def write_output(values):
    for key, value in values.items():
        print(f"{key}={value}")
    if os.getenv("GITHUB_OUTPUT"):
        with open(os.environ["GITHUB_OUTPUT"], "a") as f:
            f.writelines(f"{k}={v}\n" for k, v in values.items())


def cmd_plan():
    env = os.environ
    pyproject = open("pyproject.toml").read()
    before = env.get("BEFORE_SHA", "")
    pin_before = None
    if before and set(before) != {"0"}:
        old = git("show", f"{before}:pyproject.toml", check=False)
        pin_before = ngsolve_pin(old) if old else None
    tags = git("tag", "--list", "v*").split()
    out = plan(
        event=env["EVENT_NAME"], ref_name=env["REF_NAME"], head_ref=env.get("HEAD_REF", ""),
        input_ngsolve=env.get("INPUT_NGSOLVE", ""),
        dry_run=env.get("INPUT_DRY_RUN", "false").lower() == "true",
        pin_now=ngsolve_pin(pyproject), pin_before=pin_before,
        published=pypi_versions(package_name(pyproject)), tags=tags,
        tags_at_head=git("tag", "--points-at", "HEAD").split())
    out["branch"] = DRYRUN_BRANCH if env.get("INPUT_DRY_RUN", "").lower() == "true" else RC_BRANCH
    write_output(out)
    if env.get("GITHUB_STEP_SUMMARY"):
        with open(env["GITHUB_STEP_SUMMARY"], "a") as f:
            f.write("### Release plan\n\n" + "\n".join(f"- `{k}`: {v}" for k, v in out.items()) + "\n")


def cmd_prepare():
    """Branch from master (= the checked out commit), bump the pin, push; print the commit."""
    ngsolve, branch = os.environ["NGSOLVE_VERSION"], os.environ["RC_BRANCH"]
    git("config", "user.name", "github-actions[bot]")
    git("config", "user.email", "41898282+github-actions[bot]@users.noreply.github.com")
    remote = git("ls-remote", "--heads", "origin", branch)
    if remote:
        git("fetch", "origin", branch)
        existing = git("show", "FETCH_HEAD:pyproject.toml")
        if ngsolve_pin(existing) == ngsolve:
            print(f"{branch} already pins NGSolve {ngsolve}: keeping it")
            write_output({"ref": git("rev-parse", "FETCH_HEAD")})
            return
        print(f"{branch} pins NGSolve {ngsolve_pin(existing)}: recreating it from master")
    git("checkout", "-B", branch)  # from the checked out master commit
    pyproject = open("pyproject.toml").read()
    series = release_series(ngsolve, git("tag", "--list", "v*").split())
    open("pyproject.toml", "w").write(bump_pin(pyproject, ngsolve))
    if os.path.exists("RELEASE_NOTES.md"):
        notes = open("RELEASE_NOTES.md").read()
        open("RELEASE_NOTES.md", "w").write(bump_notes(notes, series, ngsolve))
    git("commit", "-a", "-m", f"release candidate: build against NGSolve {ngsolve}\n\n"
        f"Only bumps the NGSolve pin (pyproject.toml) and the release notes header to {series}.\n"
        f"Created by .github/workflows/pypi.yml; merging this branch into {MASTER} releases {series}.")
    git("push", "--force", "origin", f"HEAD:refs/heads/{branch}")  # the branch belongs to the automation
    write_output({"ref": git("rev-parse", "HEAD")})


if __name__ == "__main__":
    {"plan": cmd_plan, "prepare": cmd_prepare}[sys.argv[1]]()
