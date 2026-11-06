#!/usr/bin/env python3
"""Fit six-predictor Maavara forests using strictly positive training responses."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split

from importlib.machinery import SourceFileLoader

base = SourceFileLoader("maavara", "scripts/01a_refit_maavara_with_without_uparea.py").load_module()
out = Path("data/processed/50_positive_only_no_uparea")
out.mkdir(parents=True, exist_ok=True)
raw = Path("data/raw/dryad_data")
training = pd.read_csv(raw / "appling_data_06202024.csv")
scaling = pd.read_csv(raw / "scaling_dataset_Qmonth_06082024.csv")
metrics = []
counts = []
for process, spec in base.PROCESS_SPECS.items():
    x, y = base.training_data(training, spec, "without_uparea")
    keep = y > 0
    counts.append(dict(process=process, n_all=len(y), n_positive=int(keep.sum()), n_excluded=int((~keep).sum())))
    x, y = x.loc[keep], y[keep]
    x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=base.TEST_SIZE, random_state=base.RANDOM_STATE)
    model = RandomForestRegressor(**base.RF_PARAMETERS)
    model.fit(x_train, y_train)
    for dataset, obs, pred_x in (("train", y_train, x_train), ("test", y_test, x_test)):
        pred = model.predict(pred_x)
        metrics.append(dict(process=process, variant="positive_only_no_uparea", dataset=dataset,
                            n=len(obs), n_predictors=x.shape[1], oob_score=model.oob_score_,
                            **base.calculate_metrics(obs, pred)))
    for month in base.MONTHS:
        monthly_x = base.scaling_data(scaling, spec, "without_uparea", month)
        pred, _ = base.predict_with_missing_rows(model, monthly_x)
        if np.any(np.isfinite(pred) & (pred < 0)):
            raise AssertionError("Negative prediction from positive-only forest")
        frame = pd.DataFrame({"COMID": scaling.COMID, f"{process}_pred_meanQ_{month}": pred})
        base.atomic_csv(frame, out / process / f"scaling_{process}_pred_{month}.csv")
    print(f"{process}: trained on {keep.sum():,}/{len(keep):,} positive rows", flush=True)
base.atomic_csv(pd.DataFrame(metrics), out / "model_metrics.csv")
base.atomic_csv(pd.DataFrame(counts), out / "training_counts.csv")
(out / "run_manifest.json").write_text(json.dumps({"design": "six predictors; strictly positive responses only",
    "rf_parameters": base.RF_PARAMETERS, "random_state": base.RANDOM_STATE,
    "test_size": base.TEST_SIZE}, indent=2))
