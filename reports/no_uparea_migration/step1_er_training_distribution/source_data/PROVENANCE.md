# Appling source assessment files

Downloaded 2026-10-30 from the [USGS Appling metabolism data release](https://doi.org/10.5066/F70864KX). Files are copied here so the adjacent Rmd can reproduce its joins and figures. They are source data and were not edited.

| File | USGS release item | Direct file URL | SHA-256 |
|---|---|---|---|
| `diagnostics.zip` | [7. Model diagnostics](https://data.usgs.gov/datacatalog/data/USGS%3A59eb9bafe4b0026a55ffe382) | [Download](https://www.sciencebase.gov/catalog/file/get/59eb9bafe4b0026a55ffe382?f=__disk__d7%2Ff0%2F11%2Fd7f011dd4f44524d439ae5fa168e0592c0317503) | `c882ade8a8d45225ca44d07679d9dfdf47081f517622150e5ecda4ac3bac99ad` |
| `site_data.tsv` | [1. Site data](https://catalog.data.gov/dataset/metabolism-estimates-for-356-u-s-rivers-2007-2017-1-site-data) | [Download](https://www.sciencebase.gov/catalog/file/get/59bff64be4b091459a5e098b?f=__disk__01%2F7a%2Fb0%2F017ab098602c222607bff7e3eb2b957063e032d9) | `9228e92d57bdf1b2d6675e53841da772684321951b177a85273ea4e3a117c880` |
| `config.zip` | [5. Model configurations](https://catalog.data.gov/dataset/metabolism-estimates-for-356-u-s-rivers-2007-2017-5-model-configurations) | [Download](https://www.sciencebase.gov/catalog/file/get/59eb9b86e4b0026a55ffe379?f=__disk__ec%2F81%2F0d%2Fec810d105483c71596a9359613d38d24068b6023) | `3d4564807f0f8092235a3bbc3fbc639c0b213bc0cacccf29f0961e3923b62811` |

The original Maavara training table remains at `data/raw/dryad_data/appling_data_06202024.csv`. It has 355 sites. The Appling release has diagnostics for 356 sites and 433 models. Its model confidence is per model; the report uses `site_min_confidence` because the Maavara training table lacks a daily model/resolution key.
