## Creation Commands
```bash
echo "3.11" > .python-version
uv add sphinx myst-parser sphinx-book-theme sphinx-autobuild
. .venv/bin/activate
sphinx-quickstart sphinx
no separate source/build directories
rmdir sphinx/_build
# replace index.rst with index.md (markdown is easier for all involved)
```

## Building Documentation
```bash
. .venv/bin/activate
sphinx-build --write-all src/_sphinx/ docs/
sphinx-autobuild src/_sphinx/ docs/ # live updates for development
```

## Admonitions

```{important}
Public data access is currently unavailable. If you are a reviewer, a link will have been provided to you.
```
Includes: attention, caution, danger, error, hint, important, note, tip, admonition