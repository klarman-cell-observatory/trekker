#!/usr/bin/env python3

import gzip
import os
import sys

import h5py
import numpy as np
from scipy.io import mmwrite
from scipy.sparse import csc_matrix


def decode_array(values):
    """Convert HDF5 byte/string arrays to normal Python strings."""
    return [
        value.decode("utf-8") if isinstance(value, bytes) else str(value)
        for value in values
    ]


def main():

    if len(sys.argv) != 3:
        print(
            "Usage: cellbender_to_10x.py "
            "<input.h5> <output_directory>",
            file=sys.stderr,
        )
        sys.exit(1)

    input_h5 = sys.argv[1]
    output_dir = sys.argv[2]

    os.makedirs(output_dir, exist_ok=True)

    print(f"Input:  {input_h5}")
    print(f"Output: {output_dir}")

    with h5py.File(input_h5, "r") as h5:

        print("H5 structure:")
        h5.visititems(
            lambda name, obj: print(
                f"  {name}"
                if isinstance(obj, h5py.Dataset)
                else f"  {name}/"
            )
        )

        # CellBender outputs are generally 10x-style HDF5 files
        # containing a matrix group.
        if "matrix" not in h5:
            raise RuntimeError(
                "Could not find 'matrix' group in CellBender H5 file."
            )

        matrix = h5["matrix"]

        data = matrix["data"][:]
        indices = matrix["indices"][:]
        indptr = matrix["indptr"][:]
        shape = tuple(matrix["shape"][:])

        # 10x HDF5 stores matrices as CSC.
        counts = csc_matrix(
            (data, indices, indptr),
            shape=shape,
        )

        # --------------------------------------------------------
        # Barcodes
        # --------------------------------------------------------

        if "barcodes" not in matrix:
            raise RuntimeError(
                "Could not find matrix/barcodes in CellBender H5."
            )

        barcodes = decode_array(matrix["barcodes"][:])

        # --------------------------------------------------------
        # Features
        # --------------------------------------------------------

        if "features" not in matrix:
            raise RuntimeError(
                "Could not find matrix/features in CellBender H5."
            )

        features = matrix["features"]

        feature_ids = decode_array(features["id"][:])
        feature_names = decode_array(features["name"][:])

        # 10x Matrix Market features.tsv.gz normally contains:
        #
        # feature_id <TAB> feature_name <TAB> feature_type
        #
        # Use feature_type when available.
        if "feature_type" in features:
            feature_types = decode_array(features["feature_type"][:])
        else:
            feature_types = ["Gene Expression"] * len(feature_ids)

    # ------------------------------------------------------------
    # Write matrix.mtx.gz
    # ------------------------------------------------------------

    matrix_path = os.path.join(
        output_dir,
        "matrix.mtx.gz",
    )

    with gzip.open(matrix_path, "wb") as output:
        mmwrite(output, counts)

    # ------------------------------------------------------------
    # Write barcodes.tsv.gz
    # ------------------------------------------------------------

    barcode_path = os.path.join(
        output_dir,
        "barcodes.tsv.gz",
    )

    with gzip.open(
        barcode_path,
        "wt",
        encoding="utf-8",
    ) as output:

        for barcode in barcodes:
            output.write(barcode + "\n")

    # ------------------------------------------------------------
    # Write features.tsv.gz
    # ------------------------------------------------------------

    feature_path = os.path.join(
        output_dir,
        "features.tsv.gz",
    )

    with gzip.open(
        feature_path,
        "wt",
        encoding="utf-8",
    ) as output:

        for feature_id, feature_name, feature_type in zip(
            feature_ids,
            feature_names,
            feature_types,
        ):

            output.write(
                f"{feature_id}\t"
                f"{feature_name}\t"
                f"{feature_type}\n"
            )

    print()
    print("Conversion complete:")
    print(f"  {matrix_path}")
    print(f"  {barcode_path}")
    print(f"  {feature_path}")

    print()
    print(f"Matrix shape: {counts.shape}")
    print(f"Non-zero entries: {counts.nnz}")


if __name__ == "__main__":
    main()