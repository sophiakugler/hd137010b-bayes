#!/bin/sh
# Make notebooks painless in git. Run once per clone:
#
#     ./scripts/setup-notebooks.sh <your-folder>
#
# Why this is per-person and not a committed .gitattributes:
#
# A clean filter runs whenever git compares the working tree with the index.
# Point it at a notebook that was committed WITH its outputs and git reports a
# difference that exists only inside the filter -- so the file shows as modified
# forever, `git merge` and `git switch` abort with "local changes would be
# overwritten", and `git add -A` silently commits the stripped version, deleting
# someone else's outputs. So each of us strips only our own notebooks and leaves
# everyone else's byte-exact.
#
# Nothing this writes is committed: it touches .git/config and .git/info/attributes.

set -e

FOLDERS="denis ian Nicolas sophia"
mine=$1

if [ -z "$mine" ]; then
    echo "usage: $0 <your-folder>   (one of: $FOLDERS)" >&2
    exit 1
fi

root=$(git rev-parse --show-toplevel)
nbstripout=$(command -v nbstripout || true)

if [ -z "$nbstripout" ]; then
    echo "nbstripout not on PATH -- install it first (uv tool install nbstripout)" >&2
    exit 1
fi

cat > "$root/.git/nbstrip.sh" <<EOF
#!/bin/sh
# Git clean filter for Jupyter notebooks (local to this clone; see .git/config).
# Anything that is not a valid notebook is passed through byte-for-byte, so
# 'git add' can never fail on it.
in=\$(mktemp) || exit 1
out=\$(mktemp) || exit 1
trap 'rm -f "\$in" "\$out"' EXIT
cat > "\$in"
if "$nbstripout" --extra-keys 'metadata.kernelspec metadata.language_info.version' < "\$in" > "\$out" 2>/dev/null; then
    cat "\$out"
else
    cat "\$in"
fi
EOF
chmod +x "$root/.git/nbstrip.sh"

git config filter.nbstripout.clean "$root/.git/nbstrip.sh"
git config filter.nbstripout.smudge cat
git config filter.nbstripout.required true
git config diff.ipynb.textconv "$nbstripout -t --extra-keys 'metadata.kernelspec metadata.language_info.version'"

{
    echo "# Written by scripts/setup-notebooks.sh -- local to this clone, not committed."
    echo "# Strip my own notebooks; leave everyone else's exactly as they committed them."
    echo "*.ipynb        filter=nbstripout diff=ipynb"
    echo ""
    for folder in $FOLDERS; do
        [ "$folder" = "$mine" ] && continue
        printf '%-14s -filter -diff\n' "$folder/**"
    done
} > "$root/.git/info/attributes"

echo "notebooks under $mine/ (and anywhere outside $FOLDERS) are stripped on commit;"
echo "every other folder is left byte-exact. Check with:  git status --short -- '*.ipynb'"
