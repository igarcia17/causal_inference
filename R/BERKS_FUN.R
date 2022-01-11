library(dagitty)
if(!suppressWarnings(require("rethinking", quietly = TRUE))) {
  drawdag <- plot
}


DAG.Berkson <- dagitty("dag {
var1 -> hosp_patient
var2 -> hosp_patient
e_z -> hosp_patient
}")

coordinates(DAG.Berkson) <- list(x = c(var1 = 1,  var2 = 3, 
                                       hosp_patient = 2, e_z = 1.25),
                                 y = c(var1 = 1, var2 = 1, 
                                       hosp_patient = 3, e_z = 3))
drawdag(DAG.Berkson)


set.seed(11)


N=100

diabetes <- rbinom(n=100, size=1, prob=0.7)
cholecystitis <- rbinom(n=100, size=1, prob=0.3) 
hosp_patients <- diabetes*0.6 + cholecystitis*0.4 + rnorm(N, mean = 0, sd = 0.2)



df_pval_estimates <- function(var1, var2, var3) {
  dataset_2var <- data.frame(var1,var2)
  
  
  dataset_3var <- data.frame(var1,var2,var3)
  
  
  est_v1_v2 <- summary(glm(var1~var2, data = dataset_2var, family = 'binomial')
  )$coefficients['var2', 'Estimate']
  pval_v1_v2 <- summary(glm(var1~var2, data = dataset_2var, family = 'binomial')
  )$coefficients['var2', 'Pr(>|z|)']
  
  
  est_v1_v2_v3 <- summary(glm(var1~var2+var3, data=dataset_3var, 
                             family = 'binomial'))$coefficients['var2', 'Estimate']
  pval_v1_v2_v3 <-summary(glm(var1~var2+var3, data=dataset_3var, 
                             family = 'binomial'))$coefficients['var2', 'Pr(>|z|)']
  
  df_final <- data.frame('values_of' = c('Estimates', 'P_values'), 'diab_chol' = 
                               c(est_v1_v2, pval_v1_v2), 'diab_chol_HP' = 
                               c(est_v1_v2_v3, pval_v1_v2_v3))
  
print(df_final)
  
{if (pval_v1_v2 > 0.05)
  cat('The variables', var1, 'and', var2, 'do not show a association with a p value of',  
      pval_v1_v2, 'and an estimate of',est_v1_v2, '\n')
  else
    cat('The variables var1 and var2 show a association with a p value of',  
        pval_v1_v2, 'and an estimate of', est_v1_v2, '\n')}

{if (pval_v1_v2_v3 > 0.05)
  cat('The variables var1 and var2, adjusting by hospital patients 
    (a collider) do not show a association with a p value of',  
      pval_v1_v2_v3, 'and an estimate of',est_v1_v2_v3, '\n')
  else
    cat('The variables', var1, 'and var2 show a association with a p value of',  
        pval_v1_v2_v3, 'and an estimate of', est_v1_v2_v3, '\n')}
}

df_pval_estimates(diabetes, cholecystitis, hosp_patients)


