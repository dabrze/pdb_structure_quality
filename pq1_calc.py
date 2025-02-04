import logging
import os
import numpy as np
import pandas as pd
import requests
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer

logging.basicConfig(
    level=logging.INFO,
    format="{asctime} {levelname}: {message}",
    datefmt="%Y-%m-%d %H:%M:%S",
    style="{",
)

PCA_COLUMNS = ["Clashscore", "Ramachandran outliers [%]", "Rotamer outliers [%]"]


def query_pdbj(sql, data_format, dest_file):
    pdbj_rest_url = "https://pdbj.org/rest/mine2_sql"
    params = {
        "q": sql,
        "format": data_format,
    }

    response = requests.get(pdbj_rest_url, params)
    response.raise_for_status()

    directory = os.path.dirname(dest_file)
    if not os.path.exists(directory):
        os.makedirs(directory)

    with open(dest_file, "wb") as handle:
        for block in response.iter_content(2048):
            handle.write(block)


def unify_names_based_on_issn(df):
    journal_catalog = {}
    grouped_df = df.groupby(["ISSN", "Journal"])

    for name, group in grouped_df:
        issn = name[0]
        journal_name = name[1]

        if issn not in journal_catalog:
            journal_catalog[issn] = journal_name
        else:
            df.loc[df.loc[:, "ISSN"] == issn, "Journal"] = journal_catalog[issn]


def clean_quality_data(df):
    cleaned_df = df.copy()
    logging.info("Initial number of structures: %d" % cleaned_df.shape[0])

    duplicates = cleaned_df.duplicated("PDB code", keep="first")
    logging.info(
        "Found {0} duplicates: {1}".format(
            np.sum(duplicates), cleaned_df.loc[duplicates, "PDB code"].values
        )
    )
    cleaned_df = cleaned_df.drop_duplicates("PDB code")
    cleaned_df = cleaned_df.set_index("PDB code")

    cleaned_df.loc[:, "Year"] = (
        cleaned_df.loc[:, "Deposition date"].apply(str).str.slice(stop=4).apply(int)
    )
    cleaned_df.loc[:, "Deposition date"] = pd.to_datetime(
        cleaned_df.loc[:, "Deposition date"], format="%Y-%m-%d", errors="ignore"
    ).dt.date
    cleaned_df.loc[:, "Contains protein"] = (
        cleaned_df.loc[:, "Protein residues"] != 0
    ) | ~(cleaned_df.loc[:, "Ramachandran outliers [%]"].isna())
    cleaned_df.loc[:, "R [%]"] = cleaned_df.loc[:, "R"] * 100
    cleaned_df.loc[:, "R free [%]"] = cleaned_df.loc[:, "R free"] * 100
    cleaned_df.loc[:, "Contains protein"] = cleaned_df.loc[:, "Protein residues"] != 0
    cleaned_df = cleaned_df.drop(
        ["Protein residues", "Nucleic acids residues", "R", "R free"], axis=1
    )

    cleaned_df.loc[cleaned_df.Journal.isna(), "Journal"] = "TO BE PUBLISHED"
    cleaned_df.loc[:, "Journal"] = cleaned_df.loc[:, "Journal"].str.upper()
    cleaned_df.loc[:, "Journal"] = cleaned_df.loc[:, "Journal"].str.replace(
        ".", " ", regex=False
    )
    cleaned_df.loc[:, "Journal"] = cleaned_df.loc[:, "Journal"].str.replace(
        ",", " ", regex=False
    )
    cleaned_df.loc[:, "Journal"] = cleaned_df.loc[:, "Journal"].str.replace(
        "   ", " ", regex=False
    )
    cleaned_df.loc[:, "Journal"] = cleaned_df.loc[:, "Journal"].str.replace(
        "  ", " ", regex=False
    )
    cleaned_df.loc[:, "Journal"] = cleaned_df.loc[:, "Journal"].str.strip()
    unify_names_based_on_issn(cleaned_df)

    logging.info("Final number of structures: %d" % cleaned_df.shape[0])
    return cleaned_df


def clip_imputed_values_at_zero(df):
    measure_columns = [
        "Clashscore",
        "Ramachandran outliers [%]",
        "Rotamer outliers [%]",
        "RSRZ outliers [%]",
        "R [%]",
        "R free [%]",
    ]

    for measure_column in measure_columns:
        df.loc[df.loc[:, measure_column] < 0, measure_column] = 0

    return df


def fit_scaler_and_pca(PCA_COLUMNS, protein_df):
    logging.info("Fitting scaler and PCA...")
    logging.info("Looking for quality outliers...")
    protein_df_filtered_columns = protein_df.loc[:, PCA_COLUMNS]
    outlier_filter = protein_df_filtered_columns[
        (protein_df_filtered_columns.loc[:, "Rotamer outliers [%]"] > 50)
        | (protein_df_filtered_columns.loc[:, "Ramachandran outliers [%]"] > 45)
        | (protein_df_filtered_columns.loc[:, "Clashscore"] > 250)
    ].index
    logging.info("Number of outliers: {0}".format(len(outlier_filter)))
    logging.info(
        "Outliers: {0}".format(
            sorted(list(protein_df_filtered_columns.loc[outlier_filter, :].index))
        )
    )
    logging.info("Preparing data for PCA...")
    df_filtered = protein_df_filtered_columns.drop(outlier_filter, axis=0)
    scaler = StandardScaler()
    scaler.fit(df_filtered)

    df_rescaled = scaler.transform(df_filtered)
    df_rescaled = pd.DataFrame(df_rescaled, columns=PCA_COLUMNS)

    logging.info("Fitting PCA...")
    pca = PCA(svd_solver="full")
    pca.fit(df_rescaled.dropna())

    return scaler, pca


def get_data_from_pdbj():
    data_file = "data/pdb_quality.csv"
    with open("data/pdbj_query.txt", "r") as file:
        sql_query = file.read()

    logging.info("Getting data from PDBj...")
    query_pdbj(sql_query, "csv", data_file)
    quality_df_raw = pd.read_csv(data_file, na_values="", keep_default_na=False)

    return quality_df_raw


def clean_data(quality_df_raw):
    logging.info("Cleaning data...")
    full_quality_df = clean_quality_data(quality_df_raw)

    logging.info("Mapping journal names...")
    mapping_df = pd.read_csv("data/journal_mapping.csv", index_col=False)
    mapping_dict = mapping_df.set_index("From")["To"].to_dict()
    full_quality_df = full_quality_df.replace(mapping_dict).copy()
    full_quality_df = full_quality_df.drop(
        ["Journal CSD Id", "Publication year"], axis=1
    )

    return full_quality_df


def impute_missing_quality_data(full_quality_df):
    logging.info("Imputing missing data...")
    best_imputer = IterativeImputer(random_state=23, max_iter=10)

    tmp_df = full_quality_df.drop(["Deposition date", "Journal", "ISSN"], axis=1)
    full_quality_imputed_df = pd.DataFrame(
        best_imputer.fit_transform(tmp_df),
        columns=tmp_df.columns,
        index=full_quality_df.index,
    )
    full_quality_imputed_df.loc[:, "Deposition date"] = full_quality_df.loc[
        :, "Deposition date"
    ]
    full_quality_imputed_df.loc[:, "Journal"] = full_quality_df.loc[:, "Journal"]
    full_quality_imputed_df.loc[:, "ISSN"] = full_quality_df.loc[:, "ISSN"]
    full_quality_imputed_df.loc[:, "Contains protein"] = full_quality_df.loc[
        :, "Contains protein"
    ].astype("bool")
    full_quality_imputed_df.loc[
        (~full_quality_imputed_df.loc[:, "Contains protein"]),
        "Ramachandran outliers [%]",
    ] = np.nan
    full_quality_imputed_df.loc[
        (~full_quality_imputed_df.loc[:, "Contains protein"]), "Rotamer outliers [%]"
    ] = np.nan

    logging.info("Clipping imputed values at zero...")
    full_quality_imputed_df = clip_imputed_values_at_zero(full_quality_imputed_df)

    return full_quality_imputed_df


def add_Q1(df, scaler, pca, is_proteins):
    protein_metrics = list(PCA_COLUMNS)
    protein_metrics.extend(["R free [%]", "RSRZ outliers [%]"])
    nucleic_metrics = ["R free [%]", "RSRZ outliers [%]", "Clashscore"]

    if is_proteins:
        full_df = df.copy().dropna(subset=protein_metrics)
        pca_df = pd.DataFrame(
            data=pca.transform(scaler.transform(full_df.loc[:, PCA_COLUMNS])),
            columns=["PC1", "PC2", "PC3"],
            index=full_df.index,
        )
        pc1_p = pca_df.loc[:, "PC1"].rank(pct=True, ascending=False)
    else:
        full_df = df.copy().dropna(subset=nucleic_metrics)
        pc1_p = full_df.loc[:, "Clashscore"].rank(pct=True, ascending=False)

    rfree_p = full_df.loc[:, "R free [%]"].rank(pct=True, ascending=False)
    zscore_p = full_df.loc[:, "RSRZ outliers [%]"].rank(pct=True, ascending=False)

    full_df.loc[:, "Q1"] = (pc1_p + rfree_p + zscore_p) / 3.0
    full_df.loc[:, "Q1 percentile"] = full_df.loc[:, "Q1"].rank(pct=True)

    return full_df


if __name__ == "__main__":
    quality_df_raw = get_data_from_pdbj()
    quality_df = clean_data(quality_df_raw)
    quality_imputed_df = impute_missing_quality_data(quality_df)

    protein_df = quality_df.loc[quality_df.loc[:, "Contains protein"], :]
    nucleic_df = quality_df.loc[~quality_df.loc[:, "Contains protein"], :]
    protein_imputed_df = quality_imputed_df.loc[
        quality_imputed_df.loc[:, "Contains protein"], :
    ]
    nucleic_imputed_df = quality_imputed_df.loc[
        ~quality_imputed_df.loc[:, "Contains protein"], :
    ]

    scaler, pca = fit_scaler_and_pca(PCA_COLUMNS, protein_df)

    logging.info("Calculating Q1 percentiles...")
    protein_imputed_q1_df = add_Q1(protein_imputed_df, scaler, pca, is_proteins=True)
    nucleic_imputed_q1_df = add_Q1(nucleic_imputed_df, None, None, is_proteins=False)

    logging.info("Saving results to files...")
    protein_imputed_q1_df.to_csv("protein_imputed_q1.csv")
    nucleic_imputed_q1_df.to_csv("nucleic_imputed_q1.csv")
