# Manual build process
1. Commit, push and run CI tests
2. git tag -m VX.X.X X.X.X
3. git push github master --tags
4. git clean -fdxn -e prototyping -e *.code-workspace
5. git clean -fdx -e prototyping -e *.code-workspace
6. python setup.py sdist
7. python setup.py bdist_wheel
8. python -m twine upload dist/*
