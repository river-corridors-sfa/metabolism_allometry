"""Read-only ER training-label and paired-prediction diagnostics."""

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("data/raw/dryad_data")
REFITS = Path("data/processed/17_maavara_refits")
OUT = Path("reports/no_uparea_migration/step1_er_label_and_shift")
OUT.mkdir(parents=True, exist_ok=True)

training = pd.read_csv(
    ROOT / "appling_data_06202024.csv",
    usecols=["site_name", "month", "ER", "ER.lower", "ER.upper", "ER.Rhat",
             "ER.n_eff", "GPP", "uparea", "discharge", "temp.water", "shortwave",
             "TCC_mean", "uparea_prop_wet"],
)
training["raw_er_positive"] = training["ER"] > 0
training["wrong_sign_interval"] = training["ER.lower"] > 0
training["interval_crosses_zero"] = (training["ER.lower"] <= 0) & (
    training["ER.upper"] >= 0
)
quality = training.groupby("raw_er_positive", as_index=False).agg(
    rows=("ER", "size"),
    sites=("site_name", "nunique"),
    median_raw_er=("ER", "median"),
    median_rhat=("ER.Rhat", "median"),
    median_n_eff=("ER.n_eff", "median"),
    interval_entirely_positive=("wrong_sign_interval", "sum"),
    interval_crosses_zero=("interval_crosses_zero", "sum"),
    median_uparea=("uparea", "median"),
    median_discharge=("discharge", "median"),
    median_water_temp=("temp.water", "median"),
    median_shortwave=("shortwave", "median"),
    median_canopy=("TCC_mean", "median"),
    median_wetland_fraction=("uparea_prop_wet", "median"),
)
quality.to_csv(OUT / "training_label_quality.csv", index=False)
training.groupby("month", as_index=False).agg(
    rows=("ER", "size"), wrong_sign_rows=("raw_er_positive", "sum"),
).to_csv(OUT / "training_wrong_sign_by_month.csv", index=False)
training.groupby("site_name", as_index=False).agg(
    rows=("ER", "size"), wrong_sign_rows=("raw_er_positive", "sum"),
).to_csv(OUT / "training_wrong_sign_by_site.csv", index=False)

scaling = pd.read_csv(ROOT / "scaling_dataset_Qmonth_06082024.csv")
scaling["COMID"] = scaling["COMID"].astype(str)
months = [
    ("Jan", "January"), ("Feb", "February"), ("Mar", "March"),
    ("Apr", "April"), ("May", "May"), ("Jun", "June"),
    ("Jul", "July"), ("Aug", "August"), ("Sep", "September"),
    ("Oct", "October"), ("Nov", "November"), ("Dec", "December"),
]
monthly = []
by_sign = []
for number, (short, long) in enumerate(months, start=1):
    cols = ["COMID", "uparea", "TCC_mean", "up_wet", f"Q_avg_{long}",
            f"ML_river_temp_{number}", f"shortwave_{short.lower()}_mean"]
    features = scaling[cols].rename(columns={
        f"Q_avg_{long}": "discharge", f"ML_river_temp_{number}": "water_temp",
        f"shortwave_{short.lower()}_mean": "shortwave",
    })
    parts = [features]
    for variant in ("with_uparea", "without_uparea"):
        p = REFITS / variant / "ER" / f"scaling_ER_pred_{short}.csv"
        prediction = pd.read_csv(p, dtype={"COMID": str}).rename(
            columns={f"ER_pred_meanQ_{short}": variant}
        )
        parts.append(prediction)
    original = pd.read_csv(
        ROOT / "scaling_ER_pred" / f"scaling_ER_pred_{short}.csv",
        dtype={"COMID": str},
    )
    original = original.rename(columns={f"ER_pred_meanQ_{short}": "downloaded"})
    parts.append(original)
    x = parts[0]
    for part in parts[1:]:
        x = x.merge(part, on="COMID", how="outer", validate="one_to_one")
    if len(x) != 175116:
        raise ValueError(f"Unexpected COMID count for {short}: {len(x)}")
    x["sign_group"] = np.select(
        [
            (x.with_uparea >= 0) & (x.without_uparea < 0),
            (x.with_uparea < 0) & (x.without_uparea >= 0),
            (x.with_uparea < 0) & (x.without_uparea < 0),
            (x.with_uparea >= 0) & (x.without_uparea >= 0),
        ],
        ["new_negative_only", "with_negative_only", "negative_both", "nonnegative_both"],
        default="missing",
    )
    x["month"] = short
    monthly.append(x[["downloaded", "with_uparea", "without_uparea"]].agg(
        ["count", "mean", "median"]
    ).reset_index().rename(columns={"index": "statistic"}).assign(month=short))
    by_sign.append(x.groupby(["month", "sign_group"], as_index=False).agg(
        rows=("COMID", "size"),
        median_downloaded=("downloaded", "median"),
        median_with=("with_uparea", "median"),
        median_without=("without_uparea", "median"),
        median_uparea=("uparea", "median"),
        median_discharge=("discharge", "median"),
        median_water_temp=("water_temp", "median"),
        median_shortwave=("shortwave", "median"),
        median_canopy=("TCC_mean", "median"),
        median_wetland_fraction=("up_wet", "median"),
    ))
pd.concat(monthly, ignore_index=True).to_csv(OUT / "prediction_distribution_by_month.csv", index=False)
pd.concat(by_sign, ignore_index=True).to_csv(OUT / "prediction_sign_groups_by_month.csv", index=False)
print("Saved ER label and prediction diagnostics to", OUT)
