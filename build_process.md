# Build Process
Package build and upload to Pypi is handled by `.github/workflows/publish.yml`. This is triggered by publication of a new tagged release on the [Github releases page](https://github.com/wtclarke/spec2nii/releases). Upload access to Pypi is via stored secret.

[Conda-forge](https://github.com/conda-forge/spec2nii-feedstock) then picks up the new Pypi package automatically.

## Old manual build process
1. Commit, push and run CI tests
2. git tag -m VX.X.X X.X.X
3. git push github master --tags
4. git clean -fdxn -e prototyping -e *.code-workspace
5. git clean -fdx -e prototyping -e *.code-workspace
6. python setup.py sdist
7. python setup.py bdist_wheel
8. python -m twine upload dist/*
