
import argparse
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, help="Path to input csv with cell type annotations")
    parser.add_argument("-o", "--output", type=str, help="Path to output csv with cell type annotations")
    return parser.parse_args()


def copy_and_rename_pciseq_annotation_table():
    args = parse_args()
    input_file = args.input
    output_file = args.output
    
    df = pd.read_csv(input_file, index_col=0)
    pd.DataFrame(data={
        "cell_id": df.index.values,
        "celltype": df["type"],
        "score": df["prob"],
    }).to_csv(output_file, index=False)
    
if __name__ == "__main__":
    copy_and_rename_pciseq_annotation_table()
