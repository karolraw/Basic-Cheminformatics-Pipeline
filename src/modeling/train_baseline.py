#Check if optuna is installed
import subprocess
import sys
try:
  import optuna
except:
  print("Optuna not found, installing...")
  subprocess.check_call([sys.executable, "-m", "pip","install","optuna"])
  import optuna

import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, accuracy_score
import pickle
import os

#change the absolute working directory to the project root folder
project_root = os.path.dirname(os.path.abspath(__file__))
root_folder = os.path.abspath(os.path.join(os.path.join(project_root, '..'), '..'))
os.chdir(root_folder)

df = pd.read_csv("data/processed/egfr_features.csv")
#First 11 columns are not fingerprints
X = np.array(df.iloc[:, 11:].values)
#map active and inactive to 1 and 0 respectively for the model
y = df["bioactivity_class"].map({"inactive": 0, "active": 1}).values

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

#Fit a random forest classifier to training data
model = RandomForestClassifier(n_estimators=100, random_state=42)
model.fit(X_train,y_train)

#Set up our random forest classifier model for hyperparameter optimization
def objective(trial):
  n_estimators = trial.suggest_int("n_estimators", 50, 500)
  max_depth = trial.suggest_int("max_depth", 3, 30)
  min_samples_split = trial.suggest_int("min_samples_split", 2, 10)
  min_samples_leaf = trial.suggest_int("min_samples_leaf", 1, 5)

  model = RandomForestClassifier(
    n_estimators = n_estimators,
    max_depth = max_depth,
    min_samples_split = min_samples_split,
    min_samples_leaf = min_samples_leaf,
    random_state=42,
    n_jobs=-1
  )

  model.fit(X_train, y_train)
  y_prob = model.predict_proba(X_test)[:, 1]

  return roc_auc_score(y_test, y_prob)

#use optuna to run the hyperparameter optimization
study = optuna.create_study(direction="maximize")
study.optimize(objective, n_trials=30)
best_params = study.best_params
print("Best Hyperparameters:", best_params)

final_model = RandomForestClassifier(**best_params, random_state=42, n_jobs=-1)
final_model.fit(X_train,y_train)

#make predictions
y_pred = final_model.predict(X_test)
y_prob = final_model.predict_proba(X_test)[:, 1]

accuracy = accuracy_score(y_test,y_pred)
roc_auc = roc_auc_score(y_test,y_prob)

print(f"Optimized Random Forest - Accuracy {accuracy:.4f}, ROC-AUC {roc_auc:.4f}")

#Save the model
with open("models/random_forest_optimized.pkl", "wb") as f:
  pickle.dump(final_model, f)
print("Model saved to models/random_forest_optimized.pkl")




