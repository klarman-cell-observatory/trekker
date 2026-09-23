## Trekker-terra

A Terra/WDL implementation of Trekker from TakaraBio.

Specific documentation and input samplesheets to the pipeline running found on their official website: https://www.takarabio.com/products/next-generation-sequencing/bioinformatics-tools/trekker-bioinformatics-solutions

Analysis Pipeline: https://www.takarabio.com/documents/User%20Guides/Trekker%20Primary%20Analysis%20Pipeline%20User%20Manual-071426.pdf

For WDL specific modifications: provide the cloud paths in the samplesheet to pull from GCP Buckets.

Additionally, this pipeline will attempt to auto-convert filtered Cellbender output .h5s to inputs that the pipeline is expecting (barcodes/ matrix/features), provide the path ending in specific Cellbender file.

e.g. gs://path/to/file/cellbender.h5
