import pathlib

import nox


_py_versions = range(9, 13)


@nox.session(python=[f"3.{v}" for v in _py_versions])
def test(session):
    py_version = session.python.replace(".", "")
    session.install(".[test]")
    session.install(".")
    session.chdir("tests")
    session.run(
        "pytest",
        "-x",
    )
