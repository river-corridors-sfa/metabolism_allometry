#!/usr/bin/env python3
"""Refit the final Maavara et al. GPP and ER models with/without uparea.

This reproduces the final seven-predictor random forests in the archived Dryad
notebooks, then fits a paired six-predictor model whose only difference is the
removal of `uparea` (total upstream watershed area).

For each process and model variant, the script:
  1. uses the notebook's 70/30 split and random seed;
  2. fits the notebook's RandomForestRegressor specification;
  3. predicts all CONUS reaches for all 12 mean-flow months;
  4. compares with-uparea predictions to the downloaded Dryad predictions; and
  5. writes training metrics, feature importance, predictions, and a manifest.

Run from the repository root with the pinned environment documented in
`requirements-maavara-models.txt`.
"""

from __future__ import annotations

import argparse
import json
import platform
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import scipy
import sklearn
from scipy.stats import pearsonr
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import train_test_split


EXPECTED_VERSIONS = {
    "numpy": "1.26.4",
    "pandas": "2.2.1",
    "scikit-learn": "1.4.2",
    "scipy": "1.13.0",
}

RANDOM_STATE = 24
TEST_SIZE = 0.30
RF_PARAMETERS = {
    "n_estimators": 500,
    "max_features": 5,
    "max_depth": 30,
    "max_samples": 0.4,
    "n_jobs": -1,
    "random_state": RANDOM_STATE,
    "oob_score": True,
}

MONTHS = {
    "Jan": {"number": 1, "name": "January", "short": "jan", "doy": 16},
    "Feb": {"number": 2, "name": "February", "short": "feb", "doy": 46},
    "Mar": {"number": 3, "name": "March", "short": "mar", "doy": 75},
    "Apr": {"number": 4, "name": "April", "short": "apr", "doy": 105},
    "May": {"number": 5, "name": "May", "short": "may", "doy": 136},
    "Jun": {"number": 6, "name": "June", "short": "jun", "doy": 166},
    "Jul": {"number": 7, "name": "July", "short": "jul", "doy": 197},
    "Aug": {"number": 8, "name": "August", "short": "aug", "doy": 228},
    "Sep": {"number": 9, "name": "September", "short": "sep", "doy": 258},
    "Oct": {"number": 10, "name": "October", "short": "oct", "doy": 289},
    "Nov": {"number": 11, "name": "November", "short": "nov", "doy": 319},
    "Dec": {"number": 12, "name": "December", "short": "dec", "doy": 350},
}


@dataclass(frozen=True)
class ProcessSpec:
    name: str
    response: str
    landcover_training: str
    landcover_scaling: str
    base_features: tuple[str, ...]


PROCESS_SPECS = {
    "GPP": ProcessSpec(
        name="GPP",
        response="GPP",
        landcover_training="uparea_prop_shrubs",
        landcover_scaling="up_shrubs",
        base_features=(
            "doy",
            "temp.water",
            "discharge",
            "shortwave",
            "uparea",
            "TCC_mean",
            "uparea_prop_shrubs",
        ),
    ),
    "ER": ProcessSpec(
        name="ER",
        response="ER",
        landcover_training="uparea_prop_wet",
        landcover_scaling="up_wet",
        base_features=(
            "doy",
            "temp.water",
            "discharge",
            "shortwave",
            "uparea",
            "TCC_mean",
            "uparea_prop_wet",
        ),
    ),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path("data/raw/dryad_data"),
        help="Directory containing the Maavara et al. Dryad files.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/processed/17_maavara_refits"),
        help="Output directory for refits and CONUS predictions.",
    )
    parser.add_argument(
        "--processes",
        nargs="+",
        choices=tuple(PROCESS_SPECS),
        default=list(PROCESS_SPECS),
    )
    parser.add_argument(
        "--variants",
        nargs="+",
        choices=("with_uparea", "without_uparea"),
        default=["with_uparea", "without_uparea"],
    )
    parser.add_argument(
        "--months",
        nargs="+",
        choices=tuple(MONTHS),
        default=list(MONTHS),
    )
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=-1,
        help="Random-forest workers. The notebook value is -1 (all cores).",
    )
    parser.add_argument(
        "--save-models",
        action="store_true",
        help="Save fitted forests with joblib; these files can be very large.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace existing prediction outputs.",
    )
    parser.add_argument(
        "--check-only",
        action="store_true",
        help="Validate versions, inputs, and feature construction without fitting.",
    )
    parser.add_argument(
        "--allow-version-mismatch",
        action="store_true",
        help="Run despite package differences; not appropriate for replication.",
    )
    return parser.parse_args()


def current_versions() -> dict[str, str]:
    return {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scikit-learn": sklearn.__version__,
        "scipy": scipy.__version__,
        "joblib": version("joblib"),
    }


def check_versions(allow_mismatch: bool) -> dict[str, str]:
    observed = current_versions()
    mismatches = {
        package: {"expected": expected, "observed": observed[package]}
        for package, expected in EXPECTED_VERSIONS.items()
        if observed[package] != expected
    }
    if mismatches and not allow_mismatch:
        details = ", ".join(
            f"{name}: expected {values['expected']}, got {values['observed']}"
            for name, values in mismatches.items()
        )
        raise RuntimeError(
            "Package versions do not match the Maavara et al. environment. "
            f"{details}. Use the pinned requirements file."
        )
    if mismatches:
        print(f"WARNING: running with version mismatches: {mismatches}")
    return observed


def require_columns(data: pd.DataFrame, columns: list[str], label: str) -> None:
    missing = sorted(set(columns) - set(data.columns))
    if missing:
        raise ValueError(f"{label} is missing required columns: {missing}")


def ensure_finite(data: pd.DataFrame, label: str) -> None:
    values = data.to_numpy(dtype=float)
    if not np.isfinite(values).all():
        bad = data.columns[~np.isfinite(values).all(axis=0)].tolist()
        raise ValueError(f"{label} contains non-finite values in: {bad}")


def feature_names(spec: ProcessSpec, variant: str) -> list[str]:
    features = list(spec.base_features)
    if variant == "without_uparea":
        features.remove("uparea")
    return features


def training_data(
    source: pd.DataFrame,
    spec: ProcessSpec,
    variant: str,
) -> tuple[pd.DataFrame, np.ndarray]:
    features = feature_names(spec, variant)
    require_columns(source, [spec.response, *features], "training data")

    x = source.loc[:, features].copy()
    y = source[spec.response].to_numpy(copy=True)
    if spec.name == "ER":
        # StreamMetabolizer ER is negative; the published model predicts its
        # positive magnitude.
        y = y * -1

    ensure_finite(x, f"{spec.name} {variant} training predictors")
    if not np.isfinite(y).all():
        raise ValueError(f"{spec.name} training response contains non-finite values")
    return x, y


def scaling_data(
    source: pd.DataFrame,
    spec: ProcessSpec,
    variant: str,
    month_abbr: str,
) -> pd.DataFrame:
    month = MONTHS[month_abbr]
    required = [
        "COMID",
        "uparea",
        "TCC_mean",
        spec.landcover_scaling,
        f"ML_river_temp_{month['number']}",
        f"shortwave_{month['short']}_mean",
        f"Q_avg_{month['name']}",
    ]
    require_columns(source, required, "CONUS scaling data")

    out = pd.DataFrame(
        {
            "doy": month["doy"],
            "temp.water": source[f"ML_river_temp_{month['number']}"],
            "discharge": source[f"Q_avg_{month['name']}"],
            "shortwave": source[f"shortwave_{month['short']}_mean"],
            "uparea": source["uparea"],
            "TCC_mean": source["TCC_mean"],
            spec.landcover_training: source[spec.landcover_scaling],
        }
    )
    out = out.loc[:, feature_names(spec, variant)]
    return out


def predict_with_missing_rows(
    model: RandomForestRegressor,
    predictors: pd.DataFrame,
) -> tuple[np.ndarray, int]:
    """Predict complete rows and retain unavailable reaches as NA.

    The archived Dryad scaling CSV has five incomplete reaches in every month
    (three missing water temperatures and two missing shortwave values), while
    the downloaded predictions are complete. The notebooks contain no
    imputation step, so the original values cannot be recovered from the
    archive. Keeping those five results as NA is the only non-invented choice,
    and the same row rule is applied to both paired model variants.
    """
    values = predictors.to_numpy(dtype=float)
    complete = np.isfinite(values).all(axis=1)
    predictions = np.full(len(predictors), np.nan, dtype=float)
    predictions[complete] = model.predict(predictors.loc[complete])
    return predictions, int((~complete).sum())


def calculate_metrics(observed: np.ndarray, predicted: np.ndarray) -> dict[str, float]:
    pearson_r = float(pearsonr(observed, predicted).statistic)
    return {
        "pearson_r": pearson_r,
        "pearson_r_squared": pearson_r**2,
        "r2_score": float(r2_score(observed, predicted)),
        "rmse": float(mean_squared_error(observed, predicted, squared=False)),
        "mae": float(mean_absolute_error(observed, predicted)),
        "bias": float(np.mean(predicted - observed)),
    }


def atomic_csv(data: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    data.to_csv(temporary, index=False)
    temporary.replace(path)


def compare_with_downloaded(
    predictions: pd.DataFrame,
    data_dir: Path,
    spec: ProcessSpec,
    month_abbr: str,
) -> dict[str, float | int | str]:
    downloaded_path = (
        data_dir
        / f"scaling_{spec.name}_pred"
        / f"scaling_{spec.name}_pred_{month_abbr}.csv"
    )
    downloaded = pd.read_csv(downloaded_path)
    require_columns(downloaded, ["COMID"], str(downloaded_path))

    downloaded_value = next(column for column in downloaded if column != "COMID")
    model_value = next(column for column in predictions if column != "COMID")
    comparison = predictions.merge(
        downloaded.rename(columns={downloaded_value: "downloaded"}),
        on="COMID",
        how="inner",
        validate="one_to_one",
    ).rename(columns={model_value: "refit"})
    n_incomplete_refit = int(comparison["refit"].isna().sum())
    comparison = comparison.dropna(subset=["downloaded", "refit"])

    metrics = calculate_metrics(
        comparison["downloaded"].to_numpy(),
        comparison["refit"].to_numpy(),
    )
    difference = comparison["refit"] - comparison["downloaded"]
    return {
        "process": spec.name,
        "month": month_abbr,
        "downloaded_file": str(downloaded_path),
        "n_downloaded": len(downloaded),
        "n_refit": len(predictions),
        "n_matched": len(comparison),
        "n_incomplete_refit": n_incomplete_refit,
        **metrics,
        "max_absolute_difference": float(difference.abs().max()),
    }


def planned_prediction_paths(args: argparse.Namespace) -> list[Path]:
    return [
        args.output_dir
        / variant
        / process
        / f"scaling_{process}_pred_{month}.csv"
        for variant in args.variants
        for process in args.processes
        for month in args.months
    ]


def main() -> None:
    args = parse_args()
    versions = check_versions(args.allow_version_mismatch)

    training_path = args.data_dir / "appling_data_06202024.csv"
    scaling_path = args.data_dir / "scaling_dataset_Qmonth_06082024.csv"
    for path in (training_path, scaling_path):
        if not path.exists():
            raise FileNotFoundError(path)

    existing = [path for path in planned_prediction_paths(args) if path.exists()]
    if existing and not args.overwrite and not args.check_only:
        raise FileExistsError(
            f"{len(existing)} planned outputs already exist. "
            "Use --overwrite to replace them."
        )

    print(f"Reading training data: {training_path}")
    training = pd.read_csv(training_path)
    print(f"Reading CONUS scaling data: {scaling_path}")
    scaling = pd.read_csv(scaling_path)
    require_columns(scaling, ["COMID"], "CONUS scaling data")
    if scaling["COMID"].duplicated().any():
        raise ValueError("CONUS scaling data contains duplicated COMIDs")

    # Validate every requested design before committing to long model fits.
    for process in args.processes:
        spec = PROCESS_SPECS[process]
        for variant in args.variants:
            x, y = training_data(training, spec, variant)
            for month in args.months:
                monthly_x = scaling_data(scaling, spec, variant, month)
                if list(monthly_x.columns) != list(x.columns):
                    raise ValueError(
                        f"Feature mismatch for {process} {variant} {month}: "
                        f"training={list(x.columns)}, scaling={list(monthly_x.columns)}"
                    )
                n_incomplete = int(
                    (~np.isfinite(monthly_x.to_numpy(dtype=float)).all(axis=1)).sum()
                )
                if n_incomplete:
                    print(
                        f"  NOTE: {process} {variant} {month} has "
                        f"{n_incomplete} incomplete scaling rows; outputs will be NA."
                    )
            print(
                f"Validated {process} {variant}: {len(x):,} training rows, "
                f"{len(x.columns)} predictors"
            )

    if args.check_only:
        print("All version, input, and feature checks passed; no models were fit.")
        return

    args.output_dir.mkdir(parents=True, exist_ok=True)
    metric_rows: list[dict[str, object]] = []
    importance_rows: list[dict[str, object]] = []
    replication_rows: list[dict[str, object]] = []

    for process in args.processes:
        spec = PROCESS_SPECS[process]
        for variant in args.variants:
            print(f"\nFitting {process} {variant} at {datetime.now().isoformat()}")
            x, y = training_data(training, spec, variant)
            x_train, x_test, y_train, y_test = train_test_split(
                x,
                y,
                test_size=TEST_SIZE,
                random_state=RANDOM_STATE,
            )

            parameters = {**RF_PARAMETERS, "n_jobs": args.n_jobs}
            model = RandomForestRegressor(**parameters)
            model.fit(x_train, np.ravel(y_train))

            for dataset, observed, predictors in (
                ("train", y_train, x_train),
                ("test", y_test, x_test),
            ):
                predicted = model.predict(predictors)
                metric_rows.append(
                    {
                        "process": process,
                        "variant": variant,
                        "dataset": dataset,
                        "n": len(observed),
                        "n_predictors": len(x.columns),
                        "oob_score": float(model.oob_score_),
                        **calculate_metrics(observed, predicted),
                    }
                )

            importance_rows.extend(
                {
                    "process": process,
                    "variant": variant,
                    "feature": feature,
                    "importance": float(importance),
                }
                for feature, importance in zip(
                    model.feature_names_in_, model.feature_importances_
                )
            )

            if args.save_models:
                model_path = args.output_dir / "models" / f"{process}_{variant}.joblib"
                model_path.parent.mkdir(parents=True, exist_ok=True)
                joblib.dump(model, model_path, compress=3)

            for month in args.months:
                print(f"  Predicting {month}")
                monthly_x = scaling_data(scaling, spec, variant, month)
                prediction_column = f"{process}_pred_meanQ_{month}"
                predicted, n_incomplete = predict_with_missing_rows(model, monthly_x)
                predictions = pd.DataFrame(
                    {
                        "COMID": scaling["COMID"],
                        prediction_column: predicted,
                    }
                )
                if n_incomplete:
                    print(f"    Retained {n_incomplete} incomplete reaches as NA")
                output_path = (
                    args.output_dir
                    / variant
                    / process
                    / f"scaling_{process}_pred_{month}.csv"
                )
                atomic_csv(predictions, output_path)

                if variant == "with_uparea":
                    replication_rows.append(
                        compare_with_downloaded(
                            predictions, args.data_dir, spec, month
                        )
                    )

            # Forests are large; release each before fitting the next pair.
            del model, x, y, x_train, x_test, y_train, y_test

    atomic_csv(pd.DataFrame(metric_rows), args.output_dir / "model_metrics.csv")
    atomic_csv(
        pd.DataFrame(importance_rows),
        args.output_dir / "feature_importance.csv",
    )
    atomic_csv(
        pd.DataFrame(replication_rows),
        args.output_dir / "replication_against_downloaded.csv",
    )

    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": str(Path(__file__).resolve()),
        "training_file": str(training_path.resolve()),
        "scaling_file": str(scaling_path.resolve()),
        "processes": args.processes,
        "variants": args.variants,
        "months": args.months,
        "random_forest_parameters": {**RF_PARAMETERS, "n_jobs": args.n_jobs},
        "test_size": TEST_SIZE,
        "versions": versions,
        "features": {
            process: {
                variant: feature_names(PROCESS_SPECS[process], variant)
                for variant in args.variants
            }
            for process in args.processes
        },
    }
    with (args.output_dir / "run_manifest.json").open("w", encoding="utf-8") as file:
        json.dump(manifest, file, indent=2)

    print(f"\nFinished. Outputs: {args.output_dir}")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise
