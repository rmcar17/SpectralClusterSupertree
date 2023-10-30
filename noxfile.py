import pathlib

import nox


_py_versions = range(7, 11)


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


@nox.session(python=["3.10"])
def testcov(session):
    session.install(*dependencies)
    session.install(".[test]")
    session.chdir("tests")
    session.run(
        "pytest",
        "--cov-report",
        "html",
        "--cov",
        "myproj",
    )
