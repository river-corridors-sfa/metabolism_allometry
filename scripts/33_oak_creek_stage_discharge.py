"""Plot one year of paired USGS 12492900 gage height and discharge observations."""

import argparse
import json
from datetime import date, timedelta
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


SITE = "12492900"
ROOT = Path(__file__).resolve().parents[1]


def rating_curve(height, coefficient, zero_stage, exponent):
    return coefficient * np.maximum(height - zero_stage, 0) ** exponent


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--end-date", type=date.fromisoformat, default=date.today())
    args = parser.parse_args()
    start_date = args.end_date - timedelta(days=365)

    url = "https://waterservices.usgs.gov/nwis/iv/?" + urlencode(
        {
            "format": "json",
            "sites": SITE,
            "startDT": start_date.isoformat(),
            "endDT": args.end_date.isoformat(),
            "parameterCd": "00060,00065",
            "siteStatus": "all",
        }
    )
    with urlopen(url, timeout=90) as response:
        data = json.load(response)

    series = {}
    for item in data["value"]["timeSeries"]:
        code = item["variable"]["variableCode"][0]["value"]
        if code in {"00060", "00065"}:
            observations = item["values"][0]["value"]
            series[code] = pd.DataFrame(
                [(row["dateTime"], row["value"], ",".join(row["qualifiers"]))
                 for row in observations],
                columns=["date_time", "value", "qualifier"],
            ).set_index("date_time")

    if set(series) != {"00060", "00065"}:
        raise ValueError("USGS did not return both discharge and gage height")

    paired = series["00065"].join(
        series["00060"], how="inner", lsuffix="_height", rsuffix="_discharge"
    ).reset_index()
    paired = paired.rename(
        columns={"value_height": "gage_height_ft", "value_discharge": "discharge_cfs"}
    )
    paired[["gage_height_ft", "discharge_cfs"]] = paired[
        ["gage_height_ft", "discharge_cfs"]
    ].apply(pd.to_numeric, errors="coerce")
    paired = paired.dropna(subset=["gage_height_ft", "discharge_cfs"])
    paired = paired[(paired.gage_height_ft > 0) & (paired.discharge_cfs >= 0)]
    if len(paired) < 10:
        raise ValueError("Too few paired observations for a fit")

    x = paired.gage_height_ft.to_numpy()
    y = paired.discharge_cfs.to_numpy()
    coefficients, _ = curve_fit(
        rating_curve,
        x,
        y,
        p0=[25, min(x) - 0.3, 2.7],
        bounds=([0, 0, 0], [1e5, min(x) - 0.001, 10]),
        maxfev=20000,
    )
    fitted = rating_curve(x, *coefficients)
    r_squared = 1 - np.sum((y - fitted) ** 2) / np.sum((y - y.mean()) ** 2)

    data_dir = ROOT / "data" / "raw" / "usgs"
    data_dir.mkdir(parents=True, exist_ok=True)
    stem = f"oak_creek_stage_discharge_{start_date}_to_{args.end_date}"
    paired.to_csv(data_dir / f"{stem}.csv", index=False)

    plt.rcParams.update({"font.size": 12, "axes.spines.top": False,
                         "axes.spines.right": False})
    fig, ax = plt.subplots(figsize=(9, 6.3), constrained_layout=True)
    ax.scatter(x, y, s=9, color="#187A8C", alpha=0.10, edgecolors="none",
               rasterized=True)
    grid = np.linspace(x.min(), x.max(), 500)
    ax.plot(grid, rating_curve(grid, *coefficients), color="#A35C43", linewidth=2.6)
    ax.set(xlabel="Gage height (ft)", ylabel="Discharge (ft3/s)")
    ax.grid(alpha=0.22)
    ax.set_axisbelow(True)
    equation = (
        f"Q = {coefficients[0]:.2f}(H - {coefficients[1]:.3f})"
        f"^{coefficients[2]:.3f}\nR2 = {r_squared:.3f}"
    )
    ax.text(0.03, 0.97, equation, transform=ax.transAxes, va="top", ha="left",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85,
                  "pad": 6})
    output = ROOT / "figures" / f"{stem}.png"
    fig.savefig(output, dpi=220)
    print(f"{output}\n{url}\npaired observations: {len(paired)}\n{equation}")


if __name__ == "__main__":
    main()
