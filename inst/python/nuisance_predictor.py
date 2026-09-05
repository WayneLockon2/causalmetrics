#!/usr/bin/env python3
"""Out-of-fold nuisance predictions for causalmetrics estimators.

Reads a CSV of (y, d, X) and writes the same rows plus out-of-fold nuisance
predictions and a `fold_id` column. Two models are supported:

  --model irm   (default) interactive regression model / AIPW, for est_aipw()
                and est_dml(model = "irm"): writes p_hat = P(D = 1 | X),
                mu0_hat = E[Y | D = 0, X], mu1_hat = E[Y | D = 1, X].
  --model plr   partially linear regression, for est_dml(model = "plr"):
                writes l_hat = E[Y | X] and m_hat = E[D | X].

The treatment must be binary 0/1 for `irm`; for `plr` it may be continuous.
Binary treatments are predicted with logistic regression (calibrated
probabilities matter because the scores divide by them); continuous targets
use the regression learner chosen with --learner (default: histogram gradient
boosting). Folds are stratified by the treatment when it is binary.

Examples
--------
    python nuisance_predictor.py --model irm \\
        --input data.csv --output data_aug.csv \\
        --y y --d d --x x1,x2,x3 --folds 5 --seed 1

    python nuisance_predictor.py --model plr --learner hgb \\
        --input data.csv --output data_aug.csv \\
        --y y --d d --x x1,x2,x3 --folds 5 --seed 1
"""
from __future__ import annotations

import argparse
import sys

import numpy as np
import pandas as pd

from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.model_selection import KFold, StratifiedKFold


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Cross-fit nuisance predictor for causalmetrics.")
    p.add_argument("--input", required=True, help="Input CSV path.")
    p.add_argument("--output", required=True, help="Output CSV path.")
    p.add_argument("--y", required=True, help="Outcome column name.")
    p.add_argument("--d", required=True, help="Treatment column name.")
    p.add_argument(
        "--x",
        required=True,
        help="Comma-separated covariate column names (e.g. 'x1,x2,x3').",
    )
    p.add_argument(
        "--model",
        choices=["irm", "plr"],
        default="irm",
        help="irm: p_hat/mu0_hat/mu1_hat (default); plr: l_hat/m_hat.",
    )
    p.add_argument(
        "--learner",
        choices=["hgb", "linear"],
        default="hgb",
        help="Regression learner for continuous targets (default hgb).",
    )
    p.add_argument("--folds", type=int, default=5, help="Number of cross-fit folds.")
    p.add_argument("--seed", type=int, default=1, help="Random seed for folds and learners.")
    return p.parse_args()


def make_regressor(learner: str, seed: int):
    if learner == "linear":
        return LinearRegression()
    return HistGradientBoostingRegressor(random_state=seed)


def is_binary(v: np.ndarray) -> bool:
    return set(np.unique(v)).issubset({0, 1})


def splitter(d: np.ndarray, n_folds: int, seed: int):
    """Stratified folds for a binary treatment, plain folds otherwise."""
    if is_binary(d):
        return StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=seed).split(d, d.astype(int))
    return KFold(n_splits=n_folds, shuffle=True, random_state=seed).split(d)


def cross_fit_irm(X, y, d, n_folds, seed, learner):
    """Out-of-fold p_hat, mu0_hat, mu1_hat, and fold_id (1-indexed)."""
    n = len(y)
    p_hat = np.full(n, np.nan)
    mu0_hat = np.full(n, np.nan)
    mu1_hat = np.full(n, np.nan)
    fold_id = np.zeros(n, dtype=int)

    for k, (train_idx, test_idx) in enumerate(splitter(d, n_folds, seed), start=1):
        fold_id[test_idx] = k
        X_tr, X_te = X[train_idx], X[test_idx]
        y_tr, d_tr = y[train_idx], d[train_idx].astype(int)

        # Propensity: logistic regression is calibrated by construction.
        clf = LogisticRegression(max_iter=1000)
        clf.fit(X_tr, d_tr)
        p_hat[test_idx] = clf.predict_proba(X_te)[:, 1]

        # Arm-specific outcome regressions (T-learner convention).
        for arm, store in ((0, mu0_hat), (1, mu1_hat)):
            m = d_tr == arm
            if m.sum() < 2:
                raise RuntimeError(f"Fold {k} training set has fewer than 2 observations with D = {arm}.")
            reg = make_regressor(learner, seed)
            reg.fit(X_tr[m], y_tr[m])
            store[test_idx] = reg.predict(X_te)

    return p_hat, mu0_hat, mu1_hat, fold_id


def cross_fit_plr(X, y, d, n_folds, seed, learner):
    """Out-of-fold l_hat = E[Y | X], m_hat = E[D | X], and fold_id (1-indexed)."""
    n = len(y)
    l_hat = np.full(n, np.nan)
    m_hat = np.full(n, np.nan)
    fold_id = np.zeros(n, dtype=int)
    binary_d = is_binary(d)

    for k, (train_idx, test_idx) in enumerate(splitter(d, n_folds, seed), start=1):
        fold_id[test_idx] = k
        X_tr, X_te = X[train_idx], X[test_idx]

        reg_y = make_regressor(learner, seed)
        reg_y.fit(X_tr, y[train_idx])
        l_hat[test_idx] = reg_y.predict(X_te)

        if binary_d:
            clf = LogisticRegression(max_iter=1000)
            clf.fit(X_tr, d[train_idx].astype(int))
            m_hat[test_idx] = clf.predict_proba(X_te)[:, 1]
        else:
            reg_d = make_regressor(learner, seed)
            reg_d.fit(X_tr, d[train_idx])
            m_hat[test_idx] = reg_d.predict(X_te)

    return l_hat, m_hat, fold_id


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
        print(f"Dropped {before - len(df)} rows with missing required values.", file=sys.stderr)

    y = df[args.y].to_numpy(dtype=float)
    d = df[args.d].to_numpy(dtype=float)
    X = df[x_cols].to_numpy(dtype=float)

    out = df.copy()
    if args.model == "irm":
        if not is_binary(d):
            sys.exit("Treatment column must be binary 0/1 for --model irm.")
        p_hat, mu0_hat, mu1_hat, fold_id = cross_fit_irm(X, y, d, args.folds, args.seed, args.learner)
        out["p_hat"] = p_hat
        out["mu0_hat"] = mu0_hat
        out["mu1_hat"] = mu1_hat
    else:
        l_hat, m_hat, fold_id = cross_fit_plr(X, y, d, args.folds, args.seed, args.learner)
        out["l_hat"] = l_hat
        out["m_hat"] = m_hat
    out["fold_id"] = fold_id
    out.to_csv(args.output, index=False)
    print(f"Wrote {len(out)} rows with {args.model} nuisance predictions to {args.output}.", file=sys.stderr)


if __name__ == "__main__":
    main()
