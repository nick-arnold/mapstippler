# Publishing mapstippler

Everything is built and validated. What remains is yours to do: the name is
reserved on PyPI by whoever uploads first, and the accounts have to be yours.

## Before anything: two decisions in `pyproject.toml`

1. **License.** It is currently MIT (permissive: anyone may use it, including
   commercially). The `LICENSE` file says "Copyright (c) 2026 Nicholas D.
   Arnold" — change the name or year if that is wrong. If you would rather have
   GPL, replace `LICENSE` and change the `license` field. Decide now: changing
   it after people have installed the package is awkward.
2. **URLs.** `USERNAME` in the three `[project.urls]` entries is a placeholder.
   Replace it with your GitHub username.

Also check `authors` — add your email if you want it public on PyPI.

## 1. GitHub

```
cd mapstippler_proj
git init
git add .
git commit -m "mapstippler 0.4.0: blue-noise and graduated dot maps"
git branch -M main
git remote add origin https://github.com/USERNAME/mapstippler.git
git push -u origin main
```

`.gitignore` already excludes build artifacts, `__pycache__`, and the GPL
reference headers (those stay out of the repo deliberately — see
`reference/compare_reference.py` for how to fetch them).

## 2. PyPI account

1. Register at https://pypi.org/account/register/
2. Enable two-factor authentication (required for uploading).
3. Create an API token: Account settings → API tokens → Add API token. Scope it
   to "Entire account" for the first upload; afterwards you can make a
   project-scoped token for `mapstippler` and delete the broad one.
4. Save the token somewhere safe. It is shown once.

## 3. Test upload first

TestPyPI is a separate copy of PyPI for rehearsals. Register separately at
https://test.pypi.org/account/register/ and make a token there too.

```
python -m twine upload --repository testpypi dist/*
```

Username is `__token__`, password is the token including the `pypi-` prefix.

Then check the page renders (README, links, classifiers) and that it installs:

```
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ mapstippler
```

The extra index is needed because TestPyPI does not carry numpy, scipy and the
rest.

## 4. Real upload

```
python -m twine upload dist/*
```

Then, from anywhere:

```
pip install mapstippler
```

## Releasing a new version later

1. Bump `version` in `pyproject.toml`.
2. `rm -rf dist build *.egg-info`
3. `python -m build`
4. `python -m twine check dist/*`
5. `python -m twine upload dist/*`

A version number can never be reused on PyPI, even if you delete the release, so
bump it for every upload.

## What is already done

- Package renamed throughout; imports are `import mapstippler`.
- Console script: `mapstippler advise ...` works without `python -m`.
- All dependencies required (numpy, scipy, pandas, geopandas, shapely, pyarrow,
  matplotlib), so one install gets the whole tool.
- Wheel and sdist build cleanly; `twine check` passes both.
- Verified by installing the wheel into a clean virtual environment and running
  a full dot map end to end.
- 46 tests pass under the new name.
