import pathlib

import nox


_py_versions = range(10, 12)


@nox.session(python=[f"3.{v}" for v in _py_versions])
def test(session):
    py_version = session.python.replace(".", "")
    session.install(".[test]")
    session.install(".")
    session.chdir("tests")
    session.run(
        "pytest",
        "-s",
        "-x",
        "--cov-report",
        f"lcov:lcov-{session.python}.info",
        "--cov"
    )
