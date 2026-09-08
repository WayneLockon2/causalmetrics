#!/usr/bin/env python3
"""First-stage networks for structured double machine learning (PyTorch).

Companion to causalmetrics::est_dml_structured(). Fits, out of fold, either

  --model structured   the generalized sigmoid model of Ye et al. (2025):
                       E[Y | X = x, T = t] = c * sigmoid(theta_0(x) + theta(x)' t),
                       where theta_0(x) and theta(x) come from ReLU networks and c is
                       a single trained scalar. Writes theta_0, theta_1..theta_m, theta_c.
  --model pure         a plain ReLU network on (x, t) (the PDL benchmark). Writes the
                       predicted outcome under every target combination given by
                       --targets (columns pred_<combo>) and under the realized
                       treatment (pred_obs).

Two ways to define training and prediction rows:

  --fold kfold                  cross-fitting: for every fold value, train on the other
                                folds (optionally restricted to the treatment cells in
                                --train-cells) and predict the rows of that fold.
  --train-col split --train-value train
                                one split: train on the rows whose column equals the
                                value and predict every row.

--group <col> repeats the whole procedure separately for each value of <col> (used for
simulation replications). The output CSV has one row per input row that received a
prediction, identified by `.row` (1-based position in the input file).

Architectures for --model structured:
  --arch two_branch   two networks (one for theta_0, one for theta_1..theta_m), the
                      empirical architecture of the paper (Figure 6).
  --arch single       one network with m + 1 outputs (theta_0, theta_1..theta_m), the
                      architecture of the paper's synthetic experiments (Appendix D).

Example
-------
    python dedl_nuisance.py --model structured --input data.csv --output theta.csv \
        --y y --t is_DP,is_LP,is_FYP --x-all-others --fold kfold \
        --train-cells 000,001,010,100,111 --hidden 20 --depth 3 --epochs 50 \
        --batch 128 --lr 0.001 --seed 1 --history history.csv
"""
from __future__ import annotations

import argparse
import sys
import time

import numpy as np
import pandas as pd
import torch
import torch.nn as nn


def mlp(p: int, out: int, hidden: int, depth: int) -> nn.Sequential:
    layers: list[nn.Module] = []
    d_in = p
    for _ in range(depth):
        layers += [nn.Linear(d_in, hidden), nn.ReLU()]
        d_in = hidden
    layers.append(nn.Linear(d_in, out))
    return nn.Sequential(*layers)


class StructuredNet(nn.Module):
    """c * sigmoid(theta_0(x) + theta(x)' t)."""

    def __init__(self, p: int, m: int, hidden: int, depth: int, arch: str, c_init: float):
        super().__init__()
        self.arch = arch
        self.m = m
        if arch == "two_branch":
            self.net_a = mlp(p, 1, hidden, depth)
            self.net_b = mlp(p, m, hidden, depth)
        else:
            self.net = mlp(p, m + 1, hidden, depth)
        self.c = nn.Parameter(torch.tensor(float(c_init)))

    def theta(self, x: torch.Tensor) -> torch.Tensor:
        if self.arch == "two_branch":
            return torch.cat([self.net_a(x), self.net_b(x)], dim=1)
        return self.net(x)

    def forward(self, x: torch.Tensor, t: torch.Tensor) -> torch.Tensor:
        th = self.theta(x)
        u = th[:, 0] + (th[:, 1:] * t).sum(dim=1)
        return self.c * torch.sigmoid(u)


class PureNet(nn.Module):
    def __init__(self, p: int, m: int, hidden: int, depth: int):
        super().__init__()
        self.net = mlp(p + m, 1, hidden, depth)

    def forward(self, x: torch.Tensor, t: torch.Tensor) -> torch.Tensor:
        return self.net(torch.cat([x, t], dim=1)).squeeze(1)


def parse_cells(spec: str | None, m: int) -> list[tuple[int, ...]] | None:
    if spec is None or spec == "":
        return None
    out = []
    for s in spec.split(","):
        s = s.strip()
        if len(s) != m or any(ch not in "01" for ch in s):
            raise SystemExit(f"bad treatment cell '{s}' for m = {m}")
        out.append(tuple(int(ch) for ch in s))
    return out


def all_cells(m: int) -> list[tuple[int, ...]]:
    return [tuple((i >> (m - 1 - k)) & 1 for k in range(m)) for i in range(2 ** m)]


def cell_labels(t: np.ndarray) -> np.ndarray:
    return np.array(["".join(str(int(v)) for v in row) for row in t])


def train_one(model: nn.Module, x: torch.Tensor, t: torch.Tensor, y: torch.Tensor,
              epochs: int, batch: int, lr: float, val_split: float, rng: torch.Generator,
              epoch_hook=None, betas=(0.9, 0.99), weight_decay: float = 0.0):
    n = x.shape[0]
    perm = torch.randperm(n, generator=rng)
    n_val = int(round(val_split * n))
    val_idx, tr_idx = perm[:n_val], perm[n_val:]
    xt, tt, yt = x[tr_idx], t[tr_idx], y[tr_idx]
    opt = torch.optim.Adam(model.parameters(), lr=lr, betas=betas, eps=1e-7, weight_decay=weight_decay)
    history = []
    n_tr = xt.shape[0]
    bs = n_tr if batch <= 0 else batch
    for ep in range(epochs):
        model.train()
        order = torch.randperm(n_tr, generator=rng)
        for i in range(0, n_tr, bs):
            b = order[i:i + bs]
            opt.zero_grad()
            loss = ((model(xt[b], tt[b]) - yt[b]) ** 2).mean()
            loss.backward()
            opt.step()
        model.eval()
        with torch.no_grad():
            tr_mse = float(((model(xt, tt) - yt) ** 2).mean())
            val_mse = float(((model(x[val_idx], t[val_idx]) - y[val_idx]) ** 2).mean()) if n_val > 0 else float("nan")
        history.append((ep + 1, tr_mse, val_mse))
        if epoch_hook is not None:
            epoch_hook(ep + 1, model)
    return history


def fit_with_restarts(make_model, x, t, y, args, rng, epoch_hook=None):
    best = None
    for r in range(max(1, args.restarts)):
        model = make_model()
        hist = train_one(model, x, t, y, args.epochs, args.batch, args.lr, args.val_split, rng,
                         epoch_hook=epoch_hook if r == 0 else None, weight_decay=args.weight_decay)
        final_val = hist[-1][2] if not np.isnan(hist[-1][2]) else hist[-1][1]
        if best is None or final_val < best[2]:
            best = (model, hist, final_val, r + 1)
        if args.val_threshold is None or final_val <= args.val_threshold:
            break
    return best


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", choices=["structured", "pure"], default="structured")
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--y", required=True)
    ap.add_argument("--t", required=True, help="comma-separated treatment indicator columns")
    ap.add_argument("--x", default=None, help="comma-separated covariate columns")
    ap.add_argument("--x-all-others", action="store_true", help="use every column not otherwise named")
    ap.add_argument("--fold", default=None, help="fold column for cross-fitting")
    ap.add_argument("--train-cells", default=None, help="restrict training rows to these cells, e.g. 000,001,010,100,111")
    ap.add_argument("--train-col", default=None)
    ap.add_argument("--train-value", default="train")
    ap.add_argument("--group", default=None, help="repeat separately per value of this column")
    ap.add_argument("--targets", default=None, help="cells to predict for --model pure (default all 2^m)")
    ap.add_argument("--arch", choices=["two_branch", "single"], default="two_branch")
    ap.add_argument("--hidden", type=int, default=20)
    ap.add_argument("--depth", type=int, default=3)
    ap.add_argument("--epochs", type=int, default=50)
    ap.add_argument("--batch", type=int, default=128, help="0 = full batch")
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--weight-decay", type=float, default=0.0)
    ap.add_argument("--val-split", type=float, default=0.1)
    ap.add_argument("--restarts", type=int, default=1, help="retrain from scratch until the validation MSE is below --val-threshold, at most this many times")
    ap.add_argument("--val-threshold", type=float, default=None)
    ap.add_argument("--c-init", default="1", help="initial value of c: a number or 'max_y'")
    ap.add_argument("--standardize", action="store_true", help="z-score the covariates on the training rows")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--threads", type=int, default=0)
    ap.add_argument("--history", default=None, help="CSV of per-epoch training and validation MSE")
    ap.add_argument("--epoch-theta", default=None, help="(structured, single fold) gzipped CSV of theta for the prediction rows after every epoch")
    args = ap.parse_args(argv)

    if args.threads > 0:
        torch.set_num_threads(args.threads)
    df = pd.read_csv(args.input)
    tvar = [c.strip() for c in args.t.split(",")]
    m = len(tvar)
    named = set(tvar + [args.y] + [c for c in [args.fold, args.train_col, args.group] if c])
    if args.x_all_others:
        xvar = [c for c in df.columns if c not in named and c != ".row"]
    elif args.x:
        xvar = [c.strip() for c in args.x.split(",")]
    else:
        raise SystemExit("give --x or --x-all-others")
    missing = [c for c in xvar + tvar + [args.y] if c not in df.columns]
    if missing:
        raise SystemExit(f"missing columns: {missing}")
    if args.fold is None and args.train_col is None:
        raise SystemExit("give --fold (cross-fitting) or --train-col/--train-value (one split)")

    train_cells = parse_cells(args.train_cells, m)
    targets = parse_cells(args.targets, m) if args.targets else all_cells(m)
    df[".row"] = np.arange(1, len(df) + 1)
    groups = [None] if args.group is None else list(pd.unique(df[args.group]))

    out_frames, hist_rows, epoch_frames = [], [], []
    t_start = time.time()
    for g in groups:
        sub = df if g is None else df[df[args.group] == g]
        cells = cell_labels(sub[tvar].values)
        if args.fold is not None:
            folds = list(pd.unique(sub[args.fold]))
            tasks = []
            for f in folds:
                is_pred = (sub[args.fold] == f).values
                is_train = ~is_pred
                if train_cells is not None:
                    keep = np.isin(cells, ["".join(map(str, c)) for c in train_cells])
                    is_train = is_train & keep
                tasks.append((f, is_train, is_pred))
        else:
            is_train = (sub[args.train_col].astype(str) == str(args.train_value)).values
            if train_cells is not None:
                is_train = is_train & np.isin(cells, ["".join(map(str, c)) for c in train_cells])
            tasks = [("all", is_train, np.ones(len(sub), dtype=bool))]

        for f, is_train, is_pred in tasks:
            seed = args.seed + (0 if g is None else int(abs(hash(str(g))) % 10_000)) + (0 if f == "all" else int(abs(hash(str(f))) % 100))
            torch.manual_seed(seed)
            rng = torch.Generator().manual_seed(seed)
            X = sub[xvar].values.astype(np.float32)
            if args.standardize:
                mu, sd = X[is_train].mean(0), X[is_train].std(0)
                sd[sd == 0] = 1.0
                X = (X - mu) / sd
            X = torch.tensor(X)
            T = torch.tensor(sub[tvar].values.astype(np.float32))
            Y = torch.tensor(sub[args.y].values.astype(np.float32))
            c_init = float(Y[is_train].max()) if args.c_init == "max_y" else float(args.c_init)
            p = X.shape[1]

            if args.model == "structured":
                make = lambda: StructuredNet(p, m, args.hidden, args.depth, args.arch, c_init)  # noqa: E731
            else:
                make = lambda: PureNet(p, m, args.hidden, args.depth)  # noqa: E731

            pred_rows = np.where(is_pred)[0]
            hook = None
            if args.epoch_theta is not None and args.model == "structured":
                def hook(ep, model, _rows=pred_rows, _f=f, _g=g):
                    model.eval()
                    with torch.no_grad():
                        th = model.theta(X[_rows]).numpy()
                    fr = pd.DataFrame(th, columns=["theta_0"] + [f"theta_{k + 1}" for k in range(m)])
                    fr.insert(0, ".row", sub[".row"].values[_rows])
                    fr.insert(0, "epoch", ep)
                    fr["theta_c"] = float(model.c.detach())
                    if _g is not None:
                        fr.insert(0, "group", _g)
                    epoch_frames.append(fr)

            model, hist, final_val, n_used = fit_with_restarts(
                make, X[is_train], T[is_train], Y[is_train], args, rng, epoch_hook=hook)
            for ep, tr, va in hist:
                hist_rows.append({"group": g, "fold": f, "epoch": ep, "train_mse": tr, "val_mse": va, "restarts": n_used})

            model.eval()
            with torch.no_grad():
                if args.model == "structured":
                    th = model.theta(X[pred_rows]).numpy()
                    fr = pd.DataFrame(th, columns=["theta_0"] + [f"theta_{k + 1}" for k in range(m)])
                    fr["theta_c"] = float(model.c.detach())
                    fr["pred_obs"] = model(X[pred_rows], T[pred_rows]).numpy()
                else:
                    fr = pd.DataFrame({"pred_obs": model(X[pred_rows], T[pred_rows]).numpy()})
                    for cell in targets:
                        tt = torch.tensor(np.tile(np.array(cell, dtype=np.float32), (len(pred_rows), 1)))
                        fr["pred_" + "".join(map(str, cell))] = model(X[pred_rows], tt).numpy()
            fr.insert(0, "fold", f)
            fr.insert(0, ".row", sub[".row"].values[pred_rows])
            if g is not None:
                fr.insert(0, "group", g)
            out_frames.append(fr)
            print(f"group={g} fold={f}: n_train={int(is_train.sum())} n_pred={len(pred_rows)} "
                  f"final train_mse={hist[-1][1]:.5f} val_mse={final_val:.5f} restarts={n_used} "
                  f"[{time.time() - t_start:.0f}s]", file=sys.stderr)

    out = pd.concat(out_frames, ignore_index=True).sort_values(".row")
    out.to_csv(args.output, index=False)
    if args.history:
        pd.DataFrame(hist_rows).to_csv(args.history, index=False)
    if args.epoch_theta and epoch_frames:
        pd.concat(epoch_frames, ignore_index=True).to_csv(args.epoch_theta, index=False, compression="gzip")
    return 0


if __name__ == "__main__":
    sys.exit(main())
