#!/usr/bin/env python3

import gzip
import os
import sys

import h5py
from scipy.io import mmwrite
from scipy.sparse import csc_matrix


def decode_array(values):
    return [
        value.decode("utf-8") if isinstance(value, bytes) else str(value)
        for value in values
    ]


def main():
    if len(sys.argv) != 3:
        print(
            "Usage: cellbender_to_10x.py <input.h5> <output_directory>",
            file=sys.stderr,
        )
        sys.exit(1)

    input_h5 = sys.argv[1]
    output_dir = sys.argv[2]

    os.makedirs(output_dir, exist_ok=True)

    print(f"Input CellBender H5: {input_h5}")
    print(f"Output directory:    {output_dir}")

    with h5py.File(input_h5, "r") as h5:

        if "matrix" not in h5:
            raise RuntimeError(
                "Could not find a 'matrix' group in the H5 file. "
                "This converter expects a 10x-style H5 structure."
            )

        matrix = h5["matrix"]

        required_matrix_fields = [
            "data",
            "indices",
            "indptr",
            "shape",
            "barcodes",
            "features",
        ]

        for field in required_matrix_fields:
            if field not in matrix:
                raise RuntimeError(
                    f"Could not find matrix/{field} in {input_h5}"
                )

        data = matrix["data"][:]
        indices = matrix["indices"][:]
        indptr = matrix["indptr"][:]
        shape = tuple(matrix["shape"][:])

        print(f"Original matrix shape: {shape}")

        counts = csc_matrix(
            (data, indices, indptr),
            shape=shape,
        )

        barcodes = decode_array(matrix["barcodes"][:])

        features = matrix["features"]

        if "id" not in features:
            raise RuntimeError(
                "Could not find matrix/features/id in the H5 file."
            )

        feature_ids = decode_array(features["id"][:])

        if "name" in features:
            feature_names = decode_array(features["name"][:])
        else:
            feature_names = feature_ids

        if "feature_type" in features:
            feature_types = decode_array(features["feature_type"][:])
        else:
            feature_types = ["Gene Expression"] * len(feature_ids)

    if counts.shape[1] != len(barcodes):
        raise RuntimeError(
            f"Matrix has {counts.shape[1]} columns, but "
            f"{len(barcodes)} barcodes were found."
        )

    if counts.shape[0] != len(feature_ids):
        raise RuntimeError(
            f"Matrix has {counts.shape[0]} rows, but "
            f"{len(feature_ids)} features were found."
        )

    if not (
        len(feature_ids)
        == len(feature_names)
        == len(feature_types)
    ):
        raise RuntimeError(
            "Feature metadata arrays have inconsistent lengths."
        )

    # ------------------------------------------------------------
    # Remove barcodes with zero corrected UMI counts.
    #
    # CellBender can produce barcodes whose corrected expression
    # vector contains no counts. Trekker calculates log_umi from
    # the total UMI count, and log(0) produces -Inf, which Trekker
    # rejects.
    # ------------------------------------------------------------

    cell_umi_totals = counts.sum(axis=0).A1

    zero_umi_mask = cell_umi_totals == 0
    zero_umi_count = int(zero_umi_mask.sum())

    print("")
    print("CellBender corrected UMI filtering:")
    print(f"  Original cells:       {counts.shape[1]}")
    print(f"  Zero-UMI cells:       {zero_umi_count}")

    if zero_umi_count > 0:

        zero_umi_barcodes = [
            barcode
            for barcode, is_zero in zip(barcodes, zero_umi_mask)
            if is_zero
        ]

        print("  Removing zero-UMI barcodes:")

        for barcode in zero_umi_barcodes:
            print(f"    {barcode}")

        keep_mask = ~zero_umi_mask

        counts = counts[:, keep_mask]
        barcodes = [
            barcode
            for barcode, keep in zip(barcodes, keep_mask)
            if keep
        ]

    print(f"  Final cells:          {counts.shape[1]}")
    print(f"  Minimum UMI:          {counts.sum(axis=0).A1.min()}")
    print(f"  Maximum UMI:          {counts.sum(axis=0).A1.max()}")

    # ------------------------------------------------------------
    # Write 10x Matrix Market output.
    # ------------------------------------------------------------

    matrix_path = os.path.join(output_dir, "matrix.mtx.gz")

    print("")
    print(f"Writing: {matrix_path}")

    with gzip.open(matrix_path, "wb") as output:
        mmwrite(output, counts)

    barcode_path = os.path.join(output_dir, "barcodes.tsv.gz")

    print(f"Writing: {barcode_path}")

    with gzip.open(
        barcode_path,
        "wt",
        encoding="utf-8",
    ) as output:
        for barcode in barcodes:
            output.write(f"{barcode}\n")

    feature_path = os.path.join(output_dir, "features.tsv.gz")

    print(f"Writing: {feature_path}")
ß
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

    print("")
    print("CellBender -> 10x conversion complete.")
    print(f"  Matrix:   {matrix_path}")
    print(f"  Barcodes: {barcode_path}")
    print(f"  Features: {feature_path}")
    print(f"  Shape:    {counts.shape}")
    print(f"  Nonzero:  {counts.nnz}")


if __name__ == "__main__":
    main()