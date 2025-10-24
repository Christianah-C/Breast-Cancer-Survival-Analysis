The analysis used the Breast Cancer METABRIC dataset, a comprehensive clinical and genomic dataset that contains information on breast cancer patients, including demographic, pathological, and treatment-related variables. Key variables include patient age, tumor size, tumor stage, hormone receptor status (ER, PR, HER2), and treatment types such as chemotherapy, radiotherapy, and type of breast surgery. The dataset also includes time-to-event (survival time) and patient status (living or deceased), making it suitable for survival analysis.

Analysis Description

The study employed a range of survival analysis techniques to explore factors influencing breast cancer patient survival:

Kaplan–Meier Estimator: Estimated and compared survival probabilities across treatment groups.

Cox Proportional Hazards Model: Examined the effect of clinical and demographic variables on survival risk.

Aalen’s Additive Model: Investigated how covariate effects vary over time.

Random Forest Survival Model: Ranked variable importance and captured complex, nonlinear effects.

Exponential and Weibull Models: Provided parametric estimates of survival curves for comparison with non-parametric results.

Together, these analyses identified tumor size, tumor stage, and patient age as the most significant predictors of survival, while estrogen receptor positivity and breast-conserving surgery were associated with improved outcomes.
