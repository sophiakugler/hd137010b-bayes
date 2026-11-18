#!/bin/sh
# Install the notebook filter for this clone. Run once, from anywhere in the repo:
#
#     ./scripts/setup-notebooks.sh
#
# .gitattributes routes every *.ipynb through it, so outputs and execution counts
# are stripped on the way into git. Notebooks then diff as code and stop
# conflicting inside base64 image blobs. Your working copies keep their outputs
# until you check them out again; the repository never stores them.
#
# Without this, the .gitattributes line is inert and your notebooks would be
# committed with their outputs, which is what everyone else's filter then strips
# back out -- so please run it before committing a notebook.
#
# Nothing it writes is committed: it touches .git/config and .git/nbstrip.sh.

set -e

root=$(git rev-parse --show-toplevel)
nbstripout=$(command -v nbstripout || true)

if [ -z "$nbstripout" ]; then
    echo "nbstripout is not on PATH -- install it first:  uv tool install nbstripout" >&2
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

echo "notebook filter installed. 'git status' should be quiet about notebooks you"
echo "have only re-run; check with:  git status --short -- '*.ipynb'"
