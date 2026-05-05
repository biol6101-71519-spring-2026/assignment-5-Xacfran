import os
import re

genes = [os.path.basename(p).replace(".model", "") for p in snakemake.input.models]

rows = []
for gene, model_file, report_file in zip(
    genes, snakemake.input.models, snakemake.input.reports
):
    # Read the one-line best model
    with open(model_file) as fh:
        best_model = fh.read().strip()

    # Pull LogL / AIC / BIC from the ModelFinder table row for the best model.
    # The table line looks like:
    #   TPM3u+R2   -12762.423   25622.847 - 0.000145   25624.679 - 0.000181   25912.430 + 0.191
    # Columns (0-based after split): 0=Model 1=LogL 2=AIC 3=±  4=w-AIC
    #                                 5=AICc  6=±  7=w-AICc  8=BIC  9=± 10=w-BIC
    aic = bic = lnl = "N/A"
    if os.path.exists(report_file):
        with open(report_file) as fh:
            for line in fh:
                parts = line.split()
                if parts and parts[0] == best_model and len(parts) >= 9:
                    try:
                        lnl = parts[1]
                        aic = parts[2]
                        bic = parts[8]
                    except IndexError:
                        pass
                    break

    rows.append((gene, best_model, lnl, aic, bic))

header = f"{'Gene':<15} {'Best Model':<25} {'lnL':<15} {'AIC':<15} {'BIC':<15}\n"
separator = "-" * 85 + "\n"

with open(snakemake.output[0], "w") as out:
    out.write("Best-fit substitution models per gene, IQ-TREE's ModelFinder)\n")
    out.write("=" * 85 + "\n")
    out.write(header)
    out.write(separator)
    for gene, model, lnl, aic, bic in rows:
        out.write(f"{gene:<15} {model:<25} {lnl:<15} {aic:<15} {bic:<15}\n")
    out.write(separator)
    out.write(f"\nTotal genes analysed: {len(rows)}\n")
    unique_models = {r[1] for r in rows}
    out.write(f"Unique models selected: {', '.join(sorted(unique_models))}\n")
