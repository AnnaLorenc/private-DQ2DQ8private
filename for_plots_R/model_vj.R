library(workflowsets)
library(tidymodels)
library(discrim)
library(ggrepel)
library(tidyverse)


#I am using lda, as PCA seems to well delineate the groups
###########
##
## PREPARE THE DATA
##
##########


##I use resampled frequencies as input:
#I will need to sample by individual when splitting data however
resampling_res_loc <- "/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/results/freqs_subsampled/productive"
annotation_loc <- "/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/data/collated_info.csv"
plot_dir <- "/Users/ania/Documents/DQ2DQ8/pipeline/260304/vj_model"


vj_all_resamples <-list()
for(i in list.files(resampling_res_loc)%>%grep(pat="tsv.gz", val=T)){
  j <- gsub(i, pat=".tsv.gz", rep="")
  vj_all_resamples[[j]] <- read_tsv(file.path(resampling_res_loc, i))%>%
    filter(param%in%c("vFamilyName", "jGeneName"))%>%
    select(param, group, contains("subsample"))%>%
    rename_with(~gsub(pat=".*_subsample_", rep="s_", .x))
    
}

vj_all_resamples <- vj_all_resamples%>%bind_rows(,.id="sample")%>%
  pivot_longer(., cols=starts_with("s_"), names_to = "subs")%>%
  replace_na(., list(value=0))

##also adding annotation


anno <- read_csv(annotation_loc)%>%
  mutate(source=substr(newname,1,1))%>%
  rename(patient=shortname)%>%
  unique

#adding anonymous IDs
anonym_patient_ids <- 
  anno %>%
  select(patient, source)%>%
  unique()%>%
  group_by(source)%>%
  mutate(anonym_patient_id=1:n()%>%str_pad(.,side = "left",width=2,"0")%>%paste0(substr(source,1,1),.))


anno <- anno %>%
  inner_join(., anonym_patient_ids%>%ungroup%>%select(patient, anonym_patient_id), by=c("patient"))%>%
  mutate(anonym_sample_id=paste0(anonym_patient_id, cells))%>%
  mutate(stage=ifelse(source=="F","stage 3","stage 2"),
         cells_long=ifelse(cells=="N","naive","memory"),
         geno=gsub(pat=".*o", rep="", genotype_short))

vj_all_resamples <- vj_all_resamples %>%
  left_join(., anno%>%select(sample_short, geno, cells_long, stage),
            by=c("sample"="sample_short"))


##some cleaning: take only V families available across all individuals:
groups_to_remove <- vj_all_resamples%>%
  group_by(group)%>%
  summarise(all=sum(value>0))%>%
  arrange(all)%>%
  filter(all< 1/3* 25*61*2 )%>%
  pull(group)

###Now I want to create a model which guesses genotype, given 
# VFam
#Vfam, Jgene
# J gene
#cell type added



vjarc <- vj_all_resamples%>%
  filter(!group%in% groups_to_remove)%>%
  group_by(sample, param, subs,geno, cells_long,stage)%>%
  mutate(valuen=value/sum(value))%>%
  pivot_wider(names_from=group, values_from=value, id_cols = c(sample, subs,geno,cells_long,stage))%>%
  ungroup()  %>%
mutate(across(where(is.character), as.factor))

colnames(vjarc) <-gsub(colnames(vjarc), pat="-", rep="_")

### small version, just averaged values per sample

vjarc_m <- vjarc%>%group_by(sample, geno, cells_long, stage)%>%
  summarise(across(starts_with("TCRBV"), mean))%>%
  ungroup


#splitting by individuals
#take 20% of individuals, stratified by genotype, cell type and stage
set.seed(42)

selected <- anno%>%select(geno, cells_long, stage, sample_short)%>%
  group_by(geno, cells_long, stage)%>%
  slice_sample(prop=0.2)


vjarc_train <- vjarc%>%filter(!sample%in% selected$sample_short)
vjarc_test  <-  vjarc%>%filter(sample%in% selected$sample_short)


make_manual_split <- function(train, test) {
  
  # Add unique row identifiers
  train <- train %>%ungroup%>% mutate(.row_id = row_number())
  test  <- test  %>%ungroup%>% mutate(.row_id = row_number() + nrow(train))
  
  # Combine into one dataset
  full_data <- bind_rows(train, test)
  
  # Create rsplit object
  split <- make_splits(
    list(
      analysis = train$.row_id,
      assessment = test$.row_id
    ),
    data = full_data
  )
  
  return(split)
}

vjarc_initial_split <- make_manual_split(vjarc_train, vjarc_test)


models = list(
lda_model = discrim_linear() %>%
  set_engine("MASS") %>%
  set_mode("classification"),

glm_model = logistic_reg() %>%
  set_engine("glm") %>%
  set_mode("classification"),

bt_model = boost_tree() %>%
  set_engine("xgboost") %>%
  set_mode("classification"),

rf_model = rand_forest() %>%
  set_engine("ranger") %>%
  set_mode("classification")
)




set.seed(42)
selected <- anno%>%select(geno, cells_long, stage, sample_short)%>%
  group_by(geno, cells_long, stage)%>%
  slice_sample(prop=0.2)

vjarc_m_train <- vjarc_m%>%filter(!sample%in% selected$sample_short)%>%ungroup
vjarc_m_test  <-  vjarc_m%>%filter(sample%in% selected$sample_short)%>%ungroup


basic_VfamJGene_m_rec <- recipe(geno ~ ., data = vjarc_m_train) %>%
  update_role(sample, new_role = "ID")%>%
  step_rm(c( cells_long, stage))%>%
  step_zv(all_predictors()) %>%
  step_nzv(all_predictors()) %>%
  step_normalize(all_numeric_predictors())

basic_Vfam_m_rec <- recipe(geno ~ ., data = vjarc_m_train) %>%
  update_role(sample, new_role = "ID")%>%
  step_rm(c(starts_with("TCRBJ"), cells_long, stage))%>%
  step_zv(all_predictors()) %>%
  step_nzv(all_predictors()) %>%
  step_normalize(all_numeric_predictors())

pca_VfamJGene_m_rec <- basic_VfamJGene_m_rec %>%
 step_pca (all_numeric(),threshold = 0.75)

pca_Vfam_m_rec <- basic_Vfam_m_rec %>%
  step_pca (all_numeric(),threshold = 0.75)

cellsstagepca_VfamJGene_m_rec <- recipe(geno ~ ., data = vjarc_m_train) %>%
  update_role(sample, new_role = "ID")%>%
  step_zv(all_predictors()) %>%
  step_nzv(all_predictors()) %>%
  step_normalize(all_numeric_predictors())%>%
  step_pca (all_numeric(),threshold = 0.75)%>%
  step_dummy(all_nominal_predictors())

cellsstage_VfamJGene_m_rec <- recipe(geno ~ ., data = vjarc_train%>%select(-subs)) %>%
update_role(sample, new_role = "ID")%>%
step_zv(all_predictors()) %>%
step_nzv(all_predictors()) %>%
step_normalize(all_numeric_predictors())%>%
step_dummy(all_nominal_predictors())




lda_models_m <- workflow_set(preproc=list(basic_VfamJGene=basic_VfamJGene_m_rec,
                                       pca_VfamJGene =pca_VfamJGene_m_rec ,
                                       basic_Vfam=basic_Vfam_m_rec,
                                       pca_Vfam= pca_Vfam_m_rec,
                                       cellsstagepca_VfamJGene =cellsstagepca_VfamJGene_m_rec ),
                           list(lda=models[[1]]), cross = FALSE)

folds_m <- vfold_cv(vjarc_m_train, v = 10,strata = "geno")
#this is to test on the same resampling each time
control <- control_resamples(save_pred = TRUE)


lda_models_m_res <- 
  lda_models_m%>%
  workflow_map("fit_resamples",seed=42,verbose=TRUE,
                control=control,
               resamples=folds_m,
               metrics=metric_set(accuracy, roc_auc, precision, f_meas, mcc))


lda_models_m_res%>%
  autoplot( select_best = TRUE)+
     geom_text_repel(aes(label = wflow_id), nudge_x = 1/8, nudge_y = 1/100) +theme_classic() +
       theme(legend.position = "none")+
  ggtitle("Model choice by 10fold xvalid, 98 samples",subtitle = "Features - which and how")

ggsave(file.path(plot_dir , "geno_from_vj_modelchoice_122.png"), width = 7.5, height = 6)

collect_metrics(lda_models_m_res, summarize = FALSE)%>% 
  mutate(wflow_id = reorder(wflow_id, .estimate)) %>% 
  ggplot(aes(x = id, y = .estimate, group = wflow_id, color = wflow_id)) + 
  geom_line(alpha = .5, linewidth = 1.25) + theme_classic()+
  theme(axis.text.x = element_text(angle=45, vjust=0.5)) +
  facet_wrap(~.metric, scales="free_y")+
  ggtitle("Model choice by 10fold xvalid, 98 samples",subtitle = "Features - which and how")

ggsave(file.path(plot_dir , "geno_from_vj_modelchoice_122B.png"), width = 7.5, height = 6)




###### Final fit basic

vjarc_m_initial_split <- make_manual_split(vjarc_m_train, vjarc_m_test)

m_workflow <- workflow()%>%
  add_recipe(basic_VfamJGene_m_rec)%>%
  add_model(models[[1]])

test_m_fit <- last_fit(m_workflow, split = vjarc_m_initial_split,  metrics=metric_set(accuracy, roc_auc, precision, f_meas, mcc))

collect_metrics( test_m_fit)

colors =c(DQ2="#FF4D4D",
          DQ2DQ8="#B455E0",
          DQ8="#2E86DE")

collect_predictions(test_m_fit ) %>%
  roc_curve(truth = geno, .pred_DQ2 ,.pred_DQ2DQ8, .pred_DQ8) %>%
  rename(genotype=.level)%>%
  ggplot(aes(x = 1 - specificity, y = sensitivity, col=genotype)) +
  geom_path(lwd=1) +
  geom_abline(lty = 3) +
  coord_equal() +
  theme_classic() + scale_color_manual(values=colors)+
  ggtitle("Predicting genotype from V and J genes frequencies",subtitle = "Modelled on the 98 samples, test on 24/122 samples")
ggsave(file.path(plot_dir , "geno_from_vj_model_122.png"), width = 7.5, height = 6)
##### seems that just V gene is enough/no improvement from other params.





# just repeated

preproc=list(basic_VfamJGene=basic_VfamJGene_m_rec,
             pca_VfamJGene =pca_VfamJGene_m_rec ,
             basic_Vfam=basic_Vfam_m_rec,
             pca_Vfam= pca_Vfam_m_rec,
             cellsstagepca_VfamJGene =cellsstagepca_VfamJGene_m_rec,
             cellsstage_VfamJGene_rec =cellsstage_VfamJGene_m_rec )

lda_models_full <- workflow_set(preproc=preproc,
                           list(lda=models[[1]]), cross = FALSE)

folds_full <- vfold_cv(vjarc_train%>%select(-subs), v = 10, strata = "geno")
#this is to test on the same resampling each time
control <- control_resamples(save_pred = TRUE)


lda_models_res_full <- 
  lda_models_full%>%
  workflow_map("fit_resamples",seed=42,verbose=TRUE,
               control=control,
               resamples=folds_full,
               metrics=metric_set(accuracy, roc_auc, precision, f_meas, mcc))



lda_models_res_full$wflow_id <- gsub(lda_models_res_full$wflow_id, pat="_lda$|_rec", rep="")

lda_models_res_full%>%autoplot( select_best = TRUE)+
  geom_text_repel(aes(label = wflow_id), nudge_x = 1/8, nudge_y = 1/100) +theme_classic() +
  theme(legend.position = "none") +
ggtitle("Model choice by 10fold xvalid, 2450 samples",subtitle = "Features - which and how")
ggsave(file.path(plot_dir , "geno_from_vj_modelchoice_3050A.png"), width = 7.5, height = 6)

collect_metrics(lda_models_res_full, summarize = FALSE)%>% 
  mutate(wflow_id = reorder(wflow_id, .estimate)) %>% 
  ggplot(aes(x = id, y = .estimate, group = wflow_id, color = wflow_id, lty=wflow_id)) + 
  geom_line(alpha = .5, linewidth = 1.25) + theme_classic()+
  theme(axis.text.x = element_text(angle=45, vjust=0.5)) +facet_wrap(~.metric, scales="free_y")+
  ggtitle("Model choice by 10fold xvalid, 2450 samples",subtitle = "Features - which and how")
ggsave(file.path(plot_dir , "geno_from_vj_modelchoice_3050B.png"), width = 7.5, height = 6)


###### basic_Vfam is the best one 

vjarc_initial_split <- make_manual_split(vjarc_train%>%select(-subs), vjarc_test%>%select(-subs))


cellsstage_VfamJGene_lda_wflow <- workflow()%>%
  add_recipe(cellsstage_VfamJGene_m_rec)%>%
  add_model(models[[1]])

test_fit <- last_fit(cellsstage_VfamJGene_lda_wflow , split = vjarc_initial_split,  metrics=metric_set(accuracy, roc_auc, precision, f_meas, mcc))

collect_metrics( test_fit)

colors =c(DQ2="#FF4D4D",
          DQ2DQ8="#B455E0",
          DQ8="#2E86DE")

collect_predictions(test_fit ) %>%
  roc_curve(truth = geno, .pred_DQ2 ,.pred_DQ2DQ8, .pred_DQ8) %>%
  rename(genotype=.level)%>%
  ggplot(aes(x = 1 - specificity, y = sensitivity, col=genotype)) +
  geom_path(lwd=1) +
  geom_abline(lty = 3) +
  coord_equal() +
  theme_classic() + scale_color_manual(values=colors)+
  ggtitle("Predicting genotype from V and J genes frequencies",subtitle = "Modelled on the 98*25 subsamples, test on 24*25/122*25 samples")
ggsave(file.path(plot_dir , "geno_from_vj_model.png"), width = 7.5, height = 6)

##Looking into the discriminants
test_fit$.workflow[[1]]$fit$fit$fit$scaling%>%as_tibble(rownames="predictor")%>%pivot_longer(cols = c(LD1, LD2), names_to = "LD")%>%
  group_by(LD)%>%
  mutate(new_LD=fct_reorder(predictor,value, .desc = TRUE))%>%
  ggplot(aes(y=value, x= new_LD))+geom_col()+coord_flip() +theme_classic() +facet_wrap(~LD, scales="free_y")+
  ggtitle("Linear discriminants from the model fitted to 3050 samples")
ggsave(file.path(plot_dir , "LD_from_vj_model.png"), width = 7.5, height = 6)


#first five predictors

lda_fit <- extract_fit_engine(test_fit$.workflow[[1]]$fit$fit)
w <- lda_fit$scaling[,1]   # LD1 weights
X <- as.matrix(vjarc_test%>%mutate( cells_long =as.numeric(factor(cells_long))-1, stage=as.numeric(factor(stage))-1)%>%select(-c(subs, sample, geno )))
contrib <- sweep(X, 2, w, `*`)
contrib_rel <- contrib / abs(rowSums(contrib))

top_five_predictors <- test_fit$.workflow[[1]]$fit$fit$fit$scaling%>%as_tibble(rownames="predictor")%>%pivot_longer(cols = c(LD1, LD2), names_to = "LD")%>%
  group_by(LD)%>%
  mutate(new_LD=fct_reorder(predictor,value, .desc = TRUE))%>%
  slice_head(n=5)%>%
  pull(predictor)%>%
  unique()

#top predictors by relative contribution
rel_top_five_predictors <-contrib_rel%>%as_tibble%>%
  add_column(geno=vjarc_test$geno, sample=vjarc_test$sample)%>%select(-c( cells_long, stage))%>%pivot_longer(cols=-c(geno,sample), names_to="predictor")%>%group_by(sample,geno, predictor)%>%
  summarise(value=mean(value))%>%group_by(geno, predictor)%>%summarise(val=mean(value))%>%group_by(geno)%>%slice_max(n=5, order_by=abs(val))%>%pull(predictor)%>%unique()


contrib_rel%>%as_tibble%>%
  add_column(geno=vjarc_test$geno, sample=vjarc_test$sample)%>%
  pivot_longer(cols=-c(geno,sample), names_to="predictor")%>%
  filter(predictor%in%top_five_predictors)%>%
  group_by(sample, geno, predictor)%>%
  summarise(value=mean(value))%>%
  ggplot(aes(x=geno, y=value, col=geno, fill=geno))+geom_boxplot(alpha=0.2,outliers = FALSE)+
  geom_jitter(height=0, width=0.2, size=2)+
  facet_wrap(~predictor, scales="free_y")+ theme_classic()+ scale_fill_manual(values=colors)+scale_color_manual(values=colors) +guides(fill="none", colour="none")+
  ggtitle("Contributions of top5 predictors, average for test sample subsamples")

ggsave(file.path(plot_dir , "mean_contrib_of_top5preds.png"), width = 7.5, height = 6)

vjarc_m%>%pivot_longer(cols=-c(geno,sample,cells_long, stage ), names_to="predictor")%>%
  filter(predictor%in%top_five_predictors)%>%
  ggplot(aes(x=geno, y=value, col=geno, fill=geno, shape=cells_long, group=geno))+geom_boxplot(alpha=0.2,outliers = FALSE)+
  geom_jitter(height=0, width=0.2, size=2)+
  facet_wrap(~predictor, scales="free_y")+ theme_classic()+ scale_fill_manual(values=colors)+scale_color_manual(values=colors) +guides(fill="none", colour="none")+scale_shape_manual(values=c(1,19))+
ggtitle("Averaged subsampled values top5 predictors,actual values")

ggsave(file.path(plot_dir , "actual_values_of_top5preds.png"), width = 7.5, height = 6)

vjarc_m%>%pivot_longer(cols=-c(geno,sample,cells_long, stage ), names_to="predictor")%>%
  filter(predictor%in%rel_top_five_predictors)%>%
  ggplot(aes(x=geno, y=value, col=geno, fill=geno, shape=cells_long, group=geno))+geom_boxplot(alpha=0.2,outliers = FALSE)+
  geom_jitter(height=0, width=0.2, size=2)+
  facet_wrap(~predictor, scales="free_y")+ theme_classic()+ scale_fill_manual(values=colors)+scale_color_manual(values=colors) +guides(fill="none", colour="none")+scale_shape_manual(values=c(1,19))+
  ggtitle("Averaged subsampled values top5 RELATIVE predictors,actual values")

ggsave(file.path(plot_dir , "actual_values_of_top5RELATIVEpreds.png"), width = 7.5, height = 6)
