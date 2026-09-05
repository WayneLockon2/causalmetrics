#!/usr/bin/env python3
"""Out-of-fold nuisance predictions for AIPW.

Reads a CSV of (y, d, X), produces stratified K-fold out-of-fold predictions
of p(X) = P(D = 1 | X), mu0(X) = E[Y | D = 0, X], and mu1(X) = E[Y | D = 1, X],
and writes a new CSV with the original columns plus `p_hat`, `mu0_hat`,
`mu1_hat`, and `fold_id`. The R-side estimator consumes this CSV directly.

Example
-------
    python nuisance_predictor.py \\
        --input data.csv --output data_aug.csv \\
        --y y --d d --x x1,x2,x3 \\
        --folds 5 --seed 1
"""
from __future__ import annotations

import argparse
import sys

import numpy as np
import pandas as pd

from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Cross-fit nuisance predictor for AIPW."
    )
    p.add_argument("--input", required=True, help="Input CSV path.")
    p.add_argument("--output", required=True, help="Output CSV path.")
    p.add_argument("--y", required=True, help="Outcome column name.")
    p.add_argument("--d", required=True, help="Treatment column name (binary 0/1).")
    p.add_argument(
        "--x",
        required=True,
        help="Comma-separated covariate column names (e.g. 'x1,x2,x3').",
    )
    p.add_argument("--folds", type=int, default=5, help="Number of cross-fit folds.")
    p.add_argument("--seed", type=int, default=1, help="Random seed for folds and learners.")
    return p.parse_args()


def cross_fit(
    X: np.ndarray,
    y: np.ndarray,
    d: np.ndarray,
    n_folds: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return out-of-fold p_hat, mu0_hat, mu1_hat, and fold_id (1-indexed)."""
    n = len(y)
    p_hat = np.full(n, np.nan)
    mu0_hat = np.full(n, np.nan)
    mu1_hat = np.full(n, np.nan)
    fold_id = np.zeros(n, dtype=int)

    skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=seed)

    for k, (train_idx, test_idx) in enumerate(skf.split(X, d), start=1):
        fold_id[test_idx] = k

        X_tr, X_te = X[train_idx], X[test_idx]
        y_tr = y[train_idx]
        d_tr = d[train_idx]

        # Propensity. Logistic regression is well-calibrated by construction;
        # AIPW depends on 1 / p, so calibration matters more than raw accuracy.
        # Swap in a calibrated tree-based classifier if you prefer.
        clf = LogisticRegression(max_iter=1000)
        clf.fit(X_tr, d_tr)
        p_hat[test_idx] = clf.predict_proba(X_te)[:, 1]

        # mu_0(X): fit only on controls in the training fold.
        m0 = d_tr == 0
        if m0.sum() < 2:
            raise RuntimeError(f"Fold {k} training set has fewer than 2 controls.")
        reg0 = HistGradientBoostingRegressor(random_state=seed)
        reg0.fit(X_tr[m0], y_tr[m0])
        mu0_hat[test_idx] = reg0.predict(X_te)

        # mu_1(X): fit only on treated in the training fold.
        m1 = d_tr == 1
        if m1.sum() < 2:
            raise RuntimeError(f"Fold {k} training set has fewer than 2 treated.")
        reg1 = HistGradientBoostingRegressor(random_state=seed)
        reg1.fit(X_tr[m1], y_tr[m1])
        mu1_hat[test_idx] = reg1.predict(X_te)

    return p_hat, mu0_hat, mu1_hat, fold_id


def main() -> None:
    args = parse_args()
    x_cols = [c.strip() for c in args.x.split(",") if c.strip()]
    if not x_cols:
        sys.exit("No covariate columns supplied via --x.")

    df = pd.read_csv(args.input)
    required = [args.y, args.d, *x_cols]
    missing = [c for c in required if c not in df.columns]
    if missing:
        sys.exit(f"Missing columns in input: {missing}")

    before = len(df)
    df = df.dropna(subset=required).reset_index(drop=True)
    if len(df) < before:
        print(
            f"Dropped {before - len(df)} rows with missing required values.",
            file=sys.stderr,
        )

    y = df[args.y].to_numpy(dtype=float)
    d = df[args.d].to_numpy(dtype=int)
    if not set(np.unique(d)).issubset({0, 1}):
        sys.exit("Treatment column must be binary 0/1.")
    X = df[x_cols].to_numpy(dtype=float)

    p_hat, mu0_hat, mu1_hat, fold_id = cross_fit(
        X, y, d, n_folds=args.folds, seed=args.seed
    )

    out = df.copy()
    out["p_hat"] = p_hat
    out["mu0_hat"] = mu0_hat
    out["mu1_hat"] = mu1_hat
    out["fold_id"] = fold_id
    out.to_csv(args.output, index=False)
    print(
        f"Wrote {len(out)} rows with nuisance predictions to {args.output}.",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
