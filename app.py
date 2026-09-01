import sys
import csv
from padelpy import from_smiles


# BACE1 QSAR BATCH PREDICTION TOOL

INPUT_CSV = "./assets/Inactive Molecule.csv"
ID_COLUMN = "pub chem SID"
SMILES_COLUMN = "Smile "
OUTPUT_CSV = "./assets/BACE1_Active_predictions.csv"

CHUNK_SIZE = 50
REQUIRED_DESCRIPTORS = ["SIC1", "SpMin4_Bhm", "SC-5"]

def compute_pic50(sic1, spmin4_bhm, sc5):
    pIC50 = (
            -15.13463
            + (13.84656 * sic1)
            + (7.33758 * spmin4_bhm)
            + (1.90156 * sc5)
    )
    IC50_nM = (10 ** (-pIC50)) * 1000000000
    return pIC50, IC50_nM


def extract_descriptors(desc_dict):
    """Pull required descriptors out of one PaDEL result dict."""
    missing = [d for d in REQUIRED_DESCRIPTORS if d not in desc_dict]
    if missing:
        raise KeyError("missing descriptor(s): {}".format(missing))
    return (
        float(desc_dict["SIC1"]),
        float(desc_dict["SpMin4_Bhm"]),
        float(desc_dict["SC-5"]),
    )


def read_input(path, id_col, smiles_col):
    rows = []
    with open(path, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        if id_col not in reader.fieldnames or smiles_col not in reader.fieldnames:
            print("ERROR: expected columns '{}' and '{}'.".format(id_col, smiles_col))
            print("Found columns: {}".format(reader.fieldnames))
            sys.exit(1)
        for row in reader:
            mol_id = row[id_col].strip()
            smi = row[smiles_col].strip()
            if smi:
                rows.append((mol_id, smi))
    return rows


def process_chunk(chunk):
    """  
    Try to run a whole chunk through PaDEL in one call (fast).    If that fails, fall back to one-by-one so we can isolate and    skip whichever SMILES is bad, instead of losing the whole chunk.    """
    results = []
    ids = [mol_id for mol_id, _ in chunk]
    smiles_list = [smi for _, smi in chunk]

    try:
        desc_list = from_smiles(smiles_list, fingerprints=False)
        if not isinstance(desc_list, list):
            desc_list = [desc_list]
        if len(desc_list) != len(chunk):
            raise ValueError("descriptor count mismatch - falling back")

        for mol_id, smi, desc in zip(ids, smiles_list, desc_list):
            results.append(_score_one(mol_id, smi, desc))
        return results

    except Exception:
        # fall back to per-molecule processing for this chunk
        for mol_id, smi in chunk:
            try:
                desc = from_smiles(smi, fingerprints=False)
                if isinstance(desc, list):
                    desc = desc[0]
                results.append(_score_one(mol_id, smi, desc))
            except Exception as e:
                results.append({
                    "ID": mol_id, "SMILES": smi, "SIC1": "", "SpMin4_Bhm": "",
                    "SC-5": "", "pIC50": "", "IC50_nM": "",
                    "status": "FAILED: {}".format(e),
                })
        return results


def _score_one(mol_id, smi, desc):
    try:
        sic1, spmin4_bhm, sc5 = extract_descriptors(desc)
        pIC50, IC50_nM = compute_pic50(sic1, spmin4_bhm, sc5)
        return {
            "ID": mol_id, "SMILES": smi,
            "SIC1": sic1, "SpMin4_Bhm": spmin4_bhm, "SC-5": sc5,
            "pIC50": round(pIC50, 4), "IC50_nM": round(IC50_nM, 4),
            "status": "OK",
        }
    except Exception as e:
        return {
            "ID": mol_id, "SMILES": smi, "SIC1": "", "SpMin4_Bhm": "",
            "SC-5": "", "pIC50": "", "IC50_nM": "",
            "status": "FAILED: {}".format(e),
        }


def main():
    print("==========================================")
    print("     BACE1 QSAR BATCH PREDICTION TOOL")
    print("==========================================")

    rows = read_input(INPUT_CSV, ID_COLUMN, SMILES_COLUMN)
    total = len(rows)
    print("\nLoaded {} molecules from {}".format(total, INPUT_CSV))

    all_results = []
    for start in range(0, total, CHUNK_SIZE):
        chunk = rows[start:start + CHUNK_SIZE]
        end = min(start + CHUNK_SIZE, total)
        print("Processing {}-{} of {}...".format(start + 1, end, total))
        all_results.extend(process_chunk(chunk))

    with open(OUTPUT_CSV, "w", newline="", encoding="utf-8") as f:
        fieldnames = ["ID", "SMILES", "SIC1", "SpMin4_Bhm", "SC-5",
                      "pIC50", "IC50_nM", "status"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_results)

    n_ok = sum(1 for r in all_results if r["status"] == "OK")
    n_fail = total - n_ok
    print("\n==========================================")
    print("Done. {} succeeded, {} failed.".format(n_ok, n_fail))
    print("Results written to {}".format(OUTPUT_CSV))
    print("==========================================")


if __name__ == "__main__":
    main()