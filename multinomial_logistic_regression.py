# %%
import pandas as pd
from sklearn.linear_model import LogisticRegressionCV
from sklearn.model_selection import StratifiedKFold
from sklearn.feature_selection import VarianceThreshold
from sklearn.model_selection import cross_val_score
import numpy as np
from sklearn.preprocessing import LabelEncoder
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.pipeline import Pipeline
from sklearn.model_selection import GridSearchCV
from sklearn.linear_model import LogisticRegression

# %% [markdown]
# SETTING UP TRAINING DATA (SANDRA'S SAMPLES)

# %%
staph_gc_private = pd.read_csv("/Users/nirushanbhag/Downloads/sandra_ML/sandraMLfinal/private_gc_matrix.csv")
private_acc2iso = pd.read_csv("/Users/nirushanbhag/Downloads/sandra_ML/sandraMLfinal/private_acc2isolation.txt", sep="\t", header=None)

# %%
private_acc2iso.columns = ['genome_name', 'Isolation_Source']
staph_gc_pa_merged_private = staph_gc_private.merge(private_acc2iso, on="genome_name") #merging with isolation source
staph_gc_pa_merged_private = staph_gc_pa_merged_private[staph_gc_pa_merged_private["Isolation_Source"] != "Oral"] #excluding oral
staph_gc_pa_merged_private["Isolation_Source"] = staph_gc_pa_merged_private["Isolation_Source"].replace({"Skin": 0, "Urine": 1, "Nasal": 2}) #label encoding
training = staph_gc_pa_merged_private #this is our training data

training

# %% [markdown]
# SETTING UP TESTING DATA (PUBLICLY AVAILABLE DATA)

# %%
staph_gc_public = pd.read_csv("/Users/nirushanbhag/Downloads/sandra_ML/sandraMLfinal/public_gc_matrix.csv")
public_acc2iso = pd.read_csv("/Users/nirushanbhag/Downloads/sandra_ML/sandraMLfinal/public_acc2isolation.txt", sep="\t", header=None)

# %%
public_acc2iso.columns = ['genome_name', 'Isolation_Source']
staph_gc_pa_mergedpublic = staph_gc_public.merge(public_acc2iso, on="genome_name") #merging with isolation site
staph_gc_pa_mergedpublic["Isolation_Source"] = staph_gc_pa_mergedpublic["Isolation_Source"].replace({"Skin": 0, "Urine": 1, "Nasal": 2}) #label encoding
testing = staph_gc_pa_mergedpublic #this is our testing data

testing

# %% [markdown]
# FEATURE SELECTION ON TRAINING AND TESTING

# %%
X_train= training.iloc[:, 1:-1]  # Feature columns
y_train = training.iloc[:, -1]  # Target column

X_test = testing.iloc[:, 1:-1] # Feature columns
X_test = X_test.filter(items = X_train.columns) #Keep Only Columns Found in X_Train
y_test = testing.iloc[:, -1] # Target column 

# Label encoding
le = LabelEncoder()
y_train = le.fit_transform(y_train)
y_test = le.transform(y_test)

# %% [markdown]
# FITTING MULTINOMIAL LOGISTIC REGRESSION ON TRAINING DATA USING GRIDSEARCHCV

# %%
#Splitting into Folds
cv = StratifiedKFold(n_splits=5, random_state = 1, shuffle = True)

#Making a Pipeline (first applying variance threshold filter and then fitting logistic regression model)
log_pipeline = Pipeline([("var_sel", VarianceThreshold()), ("fit_model", LogisticRegression(multi_class= 'multinomial', penalty='l2', solver='saga', max_iter=10000, random_state=1))])

#Defining a parameter grid to search (for gridsearch)
param_grid = {
    'var_sel__threshold': [0.05, 0.10, 0.15, 0.20],
    'fit_model__C': np.logspace(-4, 4, 10)
}

#Passing in the pipeline to GridSearch along with the defined 5 folds from CV 
search = GridSearchCV(log_pipeline, param_grid, cv = cv, scoring = "neg_log_loss", n_jobs = -1)

search.fit(X_train, y_train)


# %% [markdown]
# RESULTS

# %%
print(f"The best parameters: {search.best_params_}")
print(f"The best score: {search.best_score_}")
print(f"The best estimator: {search.best_estimator_}")

# %%
y_train_preds = search.predict(X_train)
comparison = (y_train_preds == y_train)
print("The number of correct predictions on the training data is", np.sum(comparison == True), len(comparison))
print("Training accuracy is", np.sum(comparison == True)/len(y_train_preds))

# %%
y_test_preds = search.predict(X_test)
comparison = (y_test_preds == y_test)
print("The number of correct predictions is", np.sum(comparison == True), "out of", len(comparison))
print("Testing accuracy is", np.sum(comparison == True)/len(y_test_preds))

# %% [markdown]
# INSPECT COEFFICIENT MATRIX

# %%
final_model = search.best_estimator_
features_mask = final_model.named_steps['var_sel'].get_support()
features_filtered = X_train.columns[features_mask].to_list()

logreg = final_model.named_steps['fit_model']

# %%
import seaborn as sns
import matplotlib.pyplot as plt

coef_matrix = logreg.coef_

coef_matrix_df = pd.DataFrame(coef_matrix, columns= features_filtered)

filt_coef_matrix_df = coef_matrix_df.loc[:, (abs(coef_matrix_df) > 0.025).any()] #first take the absolute value of all the values, then find the columns where in any of the four rows, there is at least one coefficient that is greater than 0.0205

index_names = ["Skin", "Urine", "Nasal"]
filt_coef_matrix_df.index = index_names

sns.clustermap(filt_coef_matrix_df, cmap = "plasma", linewidths = 0.5, xticklabels = True)
plt.title("Clustermap of Top LASSO Logistic Regression Coefficients Across 3 Staph Epi Niches", loc = "left")

filt_coef_listy = filt_coef_matrix_df.columns.to_list()
len(filt_coef_listy)

print(filt_coef_matrix_df)

# %% [markdown]
# LINK GENE CLUSTERS WITH COG FUNCTION

# %%
staph_gc_summary = pd.read_csv("/Users/nirushanbhag/Downloads/sandra_ML/sandraMLfinal/staph_PG_gene_clusters_summary.txt", sep = "\t") #Load the gene cluster summary file
staph_gc_summary.head() #Display the first few rows of the gene cluster summary file

gc_in_model = staph_gc_summary[staph_gc_summary["gene_cluster_id"].isin(filt_coef_listy)] #Model with gene clusters that are in the filtered coefficient list

gc_in_model                                                                                                                        


