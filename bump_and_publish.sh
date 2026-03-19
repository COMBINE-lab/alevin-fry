#!/usr/bin/env bash
set -euo pipefail

die() {
    echo "error: $*" >&2
    exit 1
}

usage() {
    cat <<'EOF'
Usage:
  ./bump_and_publish.sh <version> [--publish] [--dry-run]
  ./bump_and_publish.sh [--publish] [--dry-run] <version>

Options:
  --publish  Publish the bumped version to crates.io after committing
  --dry-run  Show the mutating commands without changing files, commits, tags, pushes, or publishing
  -h, --help Show this help message
EOF
}

print_cmd() {
    printf '+'
    printf ' %q' "$@"
    printf '\n'
}

run() {
    print_cmd "$@"
    if [[ "$DRY_RUN" == true ]]; then
        return 0
    fi
    "$@"
}

require_cmd() {
    command -v "$1" >/dev/null 2>&1 || die "required command not found: $1"
}

sanitize_manifest() {
    local input="$1"
    local output="$2"
    if command -v cargo-sanitize >/dev/null 2>&1; then
        cargo-sanitize -i "$input" -o "$output"
    else
        cargo sanitize -i "$input" -o "$output"
    fi
}

run_with_sanitized_manifest() {
    local backup
    backup="$(mktemp "${TMPDIR:-/tmp}/alevin-fry-Cargo.toml.XXXXXX")"
    cp "$ROOT_CARGO" "$backup"

    (
        set -euo pipefail
        trap 'mv "$backup" "$ROOT_CARGO"' EXIT
        sanitize_manifest "$backup" "$ROOT_CARGO"
        print_cmd "$@"
        "$@"
    )
}

VERSION=""
PUBLISH=false
DRY_RUN=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --publish)
            PUBLISH=true
            ;;
        --dry-run)
            DRY_RUN=true
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        -*)
            die "unknown option: $1"
            ;;
        *)
            if [[ -n "$VERSION" ]]; then
                die "version specified more than once"
            fi
            VERSION="$1"
            ;;
    esac
    shift
done

[[ -n "$VERSION" ]] || {
    usage
    exit 1
}

if ! [[ "$VERSION" =~ ^[0-9]+\.[0-9]+\.[0-9]+([+-][0-9A-Za-z.-]+)*$ ]]; then
    die "version must look like X.Y.Z, optionally with prerelease/build suffixes"
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

ROOT_CARGO="Cargo.toml"
LOCKFILE="Cargo.lock"
TAG="v${VERSION}"

[[ -f "$ROOT_CARGO" ]] || die "not found: $ROOT_CARGO"
[[ -f "$LOCKFILE" ]] || die "not found: $LOCKFILE"

require_cmd git
require_cmd cargo
require_cmd mktemp
if ! command -v cargo-sanitize >/dev/null 2>&1 && ! cargo sanitize --help >/dev/null 2>&1; then
    die "cargo-sanitize is required; install it with: cargo install cargo-sanitize"
fi

CURRENT_VERSION="$(sed -n 's/^version = "\(.*\)"/\1/p' "$ROOT_CARGO" | head -1)"
[[ -n "$CURRENT_VERSION" ]] || die "could not determine current crate version from $ROOT_CARGO"

if [[ "$CURRENT_VERSION" == "$VERSION" ]]; then
    die "crate version is already $VERSION"
fi

if [[ -n "$(git status --porcelain)" ]]; then
    die "working tree is not clean; commit or stash existing changes first"
fi

CURRENT_BRANCH="$(git branch --show-current)"
[[ -n "$CURRENT_BRANCH" ]] || die "HEAD is detached; check out a branch before releasing"

if ! UPSTREAM_BRANCH="$(git rev-parse --abbrev-ref --symbolic-full-name '@{upstream}' 2>/dev/null)"; then
    die "current branch does not have an upstream configured"
fi

LOCAL_HEAD="$(git rev-parse HEAD)"
UPSTREAM_HEAD="$(git rev-parse "$UPSTREAM_BRANCH")"
MERGE_BASE="$(git merge-base HEAD "$UPSTREAM_BRANCH")"

if [[ "$LOCAL_HEAD" != "$UPSTREAM_HEAD" || "$LOCAL_HEAD" != "$MERGE_BASE" ]]; then
    die "branch $CURRENT_BRANCH is not in sync with $UPSTREAM_BRANCH"
fi

if git rev-parse "$TAG" >/dev/null 2>&1; then
    die "tag $TAG already exists locally"
fi

if git ls-remote --exit-code --tags origin "refs/tags/$TAG" >/dev/null 2>&1; then
    die "tag $TAG already exists on origin"
fi

echo "Current version           : $CURRENT_VERSION"
echo "New version               : $VERSION"
echo "Tag                       : $TAG"
echo "Current branch            : $CURRENT_BRANCH"
echo "Upstream branch           : $UPSTREAM_BRANCH"
if [[ "$PUBLISH" == true ]]; then
    echo "Publish to crates.io      : yes"
else
    echo "Publish to crates.io      : no"
fi
if [[ "$DRY_RUN" == true ]]; then
    echo "Dry-run                   : yes"
else
    echo "Dry-run                   : no"
fi
echo

echo "Preflight checks"
print_cmd cargo check --locked -q
cargo check --locked -q

echo "Checking publishable sanitized manifest"
run_with_sanitized_manifest cargo publish --dry-run --locked --allow-dirty

echo
echo "Updating $ROOT_CARGO"
echo "  version: $CURRENT_VERSION -> $VERSION"

if [[ "$DRY_RUN" == false ]]; then
    sed -i.bak "1,/^version = /s/^version = \".*\"/version = \"${VERSION}\"/" "$ROOT_CARGO"
    rm -f "${ROOT_CARGO}.bak"
fi

UPDATED_VERSION="$(sed -n 's/^version = "\(.*\)"/\1/p' "$ROOT_CARGO" | head -1)"
if [[ "$DRY_RUN" == false ]]; then
    [[ "$UPDATED_VERSION" == "$VERSION" ]] || die "crate version update failed"
else
    echo "Dry-run: would rewrite $ROOT_CARGO"
fi

run cargo check -q
run cargo check --locked -q
run git add "$ROOT_CARGO" "$LOCKFILE"
run git commit -m "chore(release): bump alevin-fry to v${VERSION}"

if [[ "$PUBLISH" == true ]]; then
    if [[ "$DRY_RUN" == true ]]; then
        print_cmd cargo publish --locked --allow-dirty
    else
        run_with_sanitized_manifest cargo publish --locked --allow-dirty
    fi
fi

run git tag -a "$TAG" -m "Release ${VERSION}"
run git push origin HEAD
run git push origin "$TAG"

if [[ "$DRY_RUN" == true ]]; then
    echo
    echo "Dry-run complete"
else
    echo
    echo "Release bump complete for v${VERSION}"
fi
