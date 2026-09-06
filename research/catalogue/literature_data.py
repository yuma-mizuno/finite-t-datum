"""Exact finite-type matrices from Mizuno, arXiv:1912.05710, Tables 1 and 2.

The labels retain the table row order. These are reference examples used for
polynomial comparisons; classification completeness is verified separately.
"""
from __future__ import annotations
from dataclasses import dataclass
import sympy as sp
z = sp.symbols("z")


@dataclass(frozen=True)
class Datum:
    label: str
    a_plus: sp.Matrix
    a_minus: sp.Matrix


def m(rows: list[list[sp.Expr]]) -> sp.Matrix:
    return sp.Matrix(rows)


DATA = (
    Datum(
        "R2-1",
        m([[1 + z**2, -z], [-z, 1 + z**2]]),
        sp.diag(1 + z**2, 1 + z**2),
    ),
    Datum(
        "R2-2",
        m([[1 + z**2, -z], [-z, 1 + z**2]]),
        sp.diag(1 - z + z**2, 1 - z + z**2),
    ),
    Datum(
        "R2-3",
        m([[1 + z**2, -z], [-z - z**5, 1 + z**6]]),
        m([[1 + z**2, 0], [-z**3, 1 + z**6]]),
    ),
    Datum(
        "R2-4",
        m([[1 + z**2, -z], [-z - z**2, 1 + z**3]]),
        sp.diag(1 - z + z**2, 1 + z**3),
    ),
    Datum(
        "R2-5",
        m([[1 + z**2, -z], [-z - z**5 - z**9, 1 + z**10]]),
        m([[1 + z**2, 0], [-z**3 - z**7, 1 + z**10]]),
    ),
    Datum(
        "R3-1",
        m(
            [
                [1 + z**2, -z, 0],
                [-z, 1 + z**2, -z],
                [0, -z, 1 + z**2],
            ]
        ),
        sp.diag(1 + z**2, 1 + z**2, 1 + z**2),
    ),
    Datum(
        "R3-2",
        m(
            [
                [1 + z**2, -z, 0],
                [-z, 1 + z**2, -z],
                [0, -z, 1 + z**2],
            ]
        ),
        sp.diag(1 - z + z**2, 1 - z + z**2, 1 - z + z**2),
    ),
    Datum(
        "R3-3",
        m(
            [
                [1 + z**2, -z, 0],
                [-z, 1 + z**2, -z],
                [0, -z - z**2, 1 + z**3],
            ]
        ),
        sp.diag(1 - z + z**2, 1 - z + z**2, 1 + z**3),
    ),
    Datum(
        "R3-4",
        m(
            [
                [1 + z**2, 0, -z],
                [-z**3, 1 + z**6, 0],
                [-z - z**7, -z**2 - z**6, 1 + z**8],
            ]
        ),
        m(
            [
                [1 + z**2, -z, 0],
                [-z - z**5, 1 + z**6, 0],
                [0, 0, 1 + z**8],
            ]
        ),
    ),
    Datum(
        "R3-5",
        m(
            [
                [1 + z**2, -z, 0],
                [-z, 1 + z**2, -z],
                [0, -z, 1 - z + z**2],
            ]
        ),
        sp.diag(1 + z**2, 1 + z**2, 1 + z**2),
    ),
    Datum(
        "R3-6",
        m(
            [
                [1 + z**2, -z, 0],
                [-z - z**5, 1 + z**6, -z**3],
                [0, -z**3, 1 + z**6],
            ]
        ),
        m(
            [
                [1 + z**2, 0, 0],
                [-z**3, 1 + z**6, 0],
                [0, 0, 1 + z**6],
            ]
        ),
    ),
    Datum(
        "R3-7",
        m(
            [
                [1 - z + z**2, -z, 0],
                [-z, 1 + z**2, 0],
                [0, 0, 1 + z**5],
            ]
        ),
        m(
            [
                [1 + z**2, 0, 0],
                [0, 1 + z**2, -z],
                [-z**2 - z**3, -z - z**4, 1 + z**5],
            ]
        ),
    ),
)

