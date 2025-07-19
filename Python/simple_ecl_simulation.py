import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Simuliamo 1.000 mutui
np.random.seed(5)
n_loans = 1000

# Simuliamo i parametri
loans = pd.DataFrame(
    {
        "Loan_ID": np.arange(1, n_loans + 1),
        "EAD": np.random.uniform(50_000, 500_000, n_loans),  # Exposure at Default
        "PD": np.random.uniform(0.005, 0.10, n_loans),  # Probability of Default
        "LGD": np.random.uniform(0.2, 0.6, n_loans),  # Loss Given Default
    }
)

# Expected credit loss
loans["ECL"] = loans["PD"] * loans["LGD"] * loans["EAD"]
total_ecl = loans["ECL"].sum()
average_ecl = loans["ECL"].mean()

print(f"📊 Total Expected Credit Loss (Portfolio): €{total_ecl:,.2f}")
print(f"📈 Average ECL per loan: €{average_ecl:,.2f}")

# Plot distribuzione ECL
plt.hist(loans["ECL"], bins=50, edgecolor="black")
plt.title("Distribution of Expected Credit Loss (ECL) per loan")
plt.xlabel("ECL (€)")
plt.ylabel("Frequency")
plt.tight_layout()
plt.show()
