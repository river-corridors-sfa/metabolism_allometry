#!/usr/bin/env python3
"""Run the approved ER training-label sensitivity without touching production outputs."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import platform
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
import sklearn
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import train_test_split


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "reports/no_uparea_migration/er_qc_sensitivity"
SOURCE = ROOT / "data/raw/dryad_data/appling_data_06202024.csv"
SCALING = ROOT / "data/raw/dryad_data/scaling_dataset_Qmonth_06082024.csv"
ORIGINAL = ROOT / "scripts/01a_refit_maavara_with_without_uparea.py"
EXPECTED = {"numpy": "1.26.4", "pandas": "2.2.1", "sklearn": "1.4.2", "scipy": "1.13.0"}

spec = importlib.util.spec_from_file_location("maavara_refit", ORIGINAL)
assert spec and spec.loader
refit = importlib.util.module_from_spec(spec)
import sys
sys.modules[spec.name] = refit
spec.loader.exec_module(refit)


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def save(data: pd.DataFrame, name: str) -> None:
    path = OUTPUT / name
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite existing result: {path}")
    data.to_csv(path, index=False)


def scores(actual: np.ndarray, predicted: np.ndarray) -> dict:
    return {
        "n": len(actual),
        "r2": r2_score(actual, predicted) if len(actual) > 1 and np.var(actual) > 0 else np.nan,
        "rmse": mean_squared_error(actual, predicted, squared=False),
        "mae": mean_absolute_error(actual, predicted),
        "bias": float(np.mean(predicted - actual)),
        "negative_n": int(np.sum(predicted < 0)),
        "zero_n": int(np.sum(predicted == 0)),
    }


def main() -> None:
    versions = {"numpy": np.__version__, "pandas": pd.__version__, "sklearn": sklearn.__version__, "scipy": scipy.__version__}
    if versions != EXPECTED:
        raise RuntimeError(f"Pinned environment required: observed {versions}; expected {EXPECTED}")
    OUTPUT.mkdir(parents=True, exist_ok=True)
    print("Reading training and scaling sources", flush=True)
    data = pd.read_csv(SOURCE)
    scaling = pd.read_csv(SCALING)
    er_spec = refit.PROCESS_SPECS["ER"]
    if data.ER.isna().any() or scaling.COMID.duplicated().any():
        raise ValueError("Unexpected missing ER or duplicate COMID")
    train_idx, test_idx = train_test_split(np.arange(len(data)), test_size=refit.TEST_SIZE,
                                          random_state=refit.RANDOM_STATE)
    valid_test_idx = test_idx[data.ER.to_numpy()[test_idx] <= 0]
    assert len(train_idx) == 343250 and len(valid_test_idx) == 134336
    y_test = -data.ER.to_numpy()[valid_test_idx]
    save(pd.DataFrame({"source_row": valid_test_idx, "observed_er_magnitude": y_test,
                       "site_name": data.site_name.to_numpy()[valid_test_idx],
                       "uparea": data.uparea.to_numpy()[valid_test_idx],
                       "discharge": data.discharge.to_numpy()[valid_test_idx],
                       "doy": data.doy.to_numpy()[valid_test_idx],
                       "longitude": data.longitude.to_numpy()[valid_test_idx]}), "heldout_common.csv")
    metadata = {"created_utc": datetime.now(timezone.utc).isoformat(),
                "python": platform.python_version(), "versions": versions,
                "source_sha256": digest(SOURCE), "scaling_sha256": digest(SCALING),
                "original_script_sha256": digest(ORIGINAL),
                "rf_parameters": refit.RF_PARAMETERS, "split_seed": refit.RANDOM_STATE,
                "test_size": refit.TEST_SIZE, "valid_test_n": len(valid_test_idx)}
    (OUTPUT / "run_manifest.json").write_text(json.dumps(metadata, indent=2))
    metric_rows = []
    support_rows = []
    variants = [
        ("with_uparea", "with_uparea"),
        ("without_uparea", "without_uparea"),
        ("strict_exclusion", "without_uparea"),
        ("near_zero", "without_uparea"),
    ]
    for label, feature_variant in variants:
        print(f"Fitting {label}: {datetime.now().isoformat()}", flush=True)
        features = refit.feature_names(er_spec, feature_variant)
        x = data.loc[:, features]
        raw_train = data.ER.to_numpy()[train_idx]
        if label == "strict_exclusion":
            keep = raw_train <= 0
        elif label == "near_zero":
            keep = raw_train <= 0.5
        else:
            keep = np.ones(len(train_idx), dtype=bool)
        fit_idx = train_idx[keep]
        y_train = -data.ER.to_numpy()[fit_idx]
        if label == "near_zero":
            y_train = np.maximum(y_train, 0)
        assert np.isfinite(x.loc[fit_idx].to_numpy(dtype=float)).all()
        assert np.isfinite(y_train).all()
        if label in ("strict_exclusion", "near_zero"):
            assert (y_train >= 0).all()
        for feature in features:
            values = x.loc[fit_idx, feature].to_numpy(dtype=float)
            support_rows.append({"variant": label, "feature": feature,
                                 "n": len(values), "min": np.min(values),
                                 "p01": np.quantile(values, .01), "median": np.median(values),
                                 "p99": np.quantile(values, .99), "max": np.max(values)})
        model = RandomForestRegressor(**refit.RF_PARAMETERS)
        model.fit(x.loc[fit_idx], y_train)
        predicted = model.predict(x.loc[valid_test_idx])
        metric_rows.append({"variant": label, "train_n": len(fit_idx),
                            "train_sites": data.site_name.iloc[fit_idx].nunique(),
                            "oob_r2": model.oob_score_, **scores(y_test, predicted)})
        save(pd.DataFrame({"source_row": valid_test_idx, "prediction": predicted}),
             f"heldout_{label}.csv")
        # Existing baseline network predictions are reused; only QC variants need network inference.
        if label in ("strict_exclusion", "near_zero"):
            for month in refit.MONTHS:
                monthly = refit.scaling_data(scaling, er_spec, feature_variant, month)
                network_pred, missing = refit.predict_with_missing_rows(model, monthly)
                assert missing == 5
                save(pd.DataFrame({"COMID": scaling.COMID,
                                   f"ER_pred_meanQ_{month}": network_pred}),
                     f"conus/{label}/scaling_ER_pred_{month}.csv")
                print(f"  Predicted {month}", flush=True)
        del model, predicted, x
    save(pd.DataFrame(metric_rows), "heldout_metrics.csv")
    save(pd.DataFrame(support_rows), "training_predictor_support.csv")
    print("Done", flush=True)


if __name__ == "__main__":
    main()
