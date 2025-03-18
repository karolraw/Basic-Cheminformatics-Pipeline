import subprocess
import sys
def install_if_missing(package):
  try:
    __import__(package)
  except ImportError:
    print(f"{package} not found. Installing...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", package])

required_packages = ["matplotlib", "seaborn", "numpy", "pandas"]

for package in required_packages:
  install_if_missing(package)

import joblib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

#Make sure it runs from the root directory
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../"))
os.chdir(project_root)

#Specify location of the model and load it
model_path = "models/random_forest_optimized.pkl"
data_path = "data/processed/egfr_features.csv"
print("Loading model...")
rf_model = joblib.load(model_path)
df = pd.read_csv(data_path)

#Drop the non-fingerprint columns for feature importance checking
X = df.iloc[:,11:]

#Get the feature importance and show the top 20
feature_importance = rf_model.feature_importances_
top_n = 20

top_indices = np.argsort(feature_importance)[::-1][:top_n]
#Plot the data for visualization
plt.figure(figsize=(10,6))
sns.barplot(x=feature_importance[top_indices], y=[f"Feature {i}" for i in top_indices])
plt.xlabel("Feature importance score")
plt.ylabel("Top features")
plt.title("Top 20 Important Features in Random Forest")
plt.show()


feature_importance_df = pd.DataFrame({
  'Feature':[f'Feature {i}' for i in range(len(rf_model.feature_importances_))],
  'Importance': rf_model.feature_importances_
})

# Sort by importance
feature_importance_df = feature_importance_df.sort_values(by='Importance', ascending=False)

# Save to CSV
feature_importance_df.to_csv('data/processed/feature_importances.csv', index=False)

print("Feature importances saved to 'data/processed/feature_importances.csv'")
