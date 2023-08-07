import pandas as pd
import numpy as np
from sklearn import metrics
import sys
import warnings


def calculate_auc(y_true, y_pred):
    # Check if y_true contains NaN values
    if np.isnan(y_true).any():
        # Generate random true labels as a substitute for NaN values

        random_labels = np.random.choice([0, 1], size=len(y_true), replace=True)
        return metrics.roc_auc_score(random_labels, y_pred)
    else:
        return metrics.roc_auc_score(y_true, y_pred)


def main():
    # Check if the CSV file path is provided as a command line argument
    if len(sys.argv) != 2:
        print("Please provide the path to the CSV file as a command line argument.")
        return

    csv_path = sys.argv[1]

    try:
        # Read the CSV file
        df = pd.read_csv(csv_path)
        y_true = df["y_true"]
        if np.isnan(y_true).any():
            warnings.warn("NaN values were found in y_true, using random y_true values!")

        # Iterate over the y_pred columns and calculate AUC for each model
        for i in range(3):
            y_pred_col = "y_pred"
            y_pred = df[y_pred_col]

            auc = calculate_auc(y_true, y_pred)
            print(f"AUC for {y_pred_col}: {auc}")

    except FileNotFoundError:
        print("The provided CSV file does not exist.")
    except Exception as e:
        print("An error occurred while processing the CSV file:", str(e))


if __name__ == "__main__":
    main()
