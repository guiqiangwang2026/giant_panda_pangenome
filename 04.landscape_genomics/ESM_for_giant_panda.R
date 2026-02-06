# 2024-4-22
# Encoding UTF-8
# R 4.3.1
# packages verison 
# install.packages("biomod2", dependencies = TRUE)
# Important packages
library(biomod2)    # version 4.2-4 
library(flexsdm)    # version 1.3.4
library(ecospat)    # version 4.0.0
# Other packages
library(dplyr)      # version 1.1.4
library(data.table) # version 1.14.8
library(terra)      # version 1.7-55
library(tidyterra)
library(ggtext)
library(caret)
library(ENMeval)
library(sf)
# SHIFTING FROM RASTER TO TERRA
# https://oceanhealthindex.org/news/raster_to_terra/
# Sys.setenv(PROJ_LIB = "/public/home/haobinhao/.conda/envs/jupyter/share/proj") # 指定proj的位置，要不然总是报错

# 设置环境变量路径
envir_path <- "./CHELSA/WGS84_output"

# 获取所有.tif文件（不区分大小写）
tif_files <- list.files(path = envir_path, 
                       pattern = "\\.tif$", 
                       full.names = TRUE, 
                       ignore.case = TRUE)

# 检查是否找到文件
if (length(tif_files) == 0) {
  stop("在指定路径下未找到任何TIFF文件")
} else {
  message(paste("找到", length(tif_files), "个TIFF文件"))
}

# 方法1：逐个读取并重命名（推荐）
envir_stack <- rast()
for (file in tif_files) {
  # 获取不带扩展名的文件名
  layer_name <- tools::file_path_sans_ext(basename(file))
  # 读取单个文件并设置名称
  r <- rast(file)
  names(r) <- layer_name
  # 添加到堆栈
  envir_stack <- c(envir_stack, r)
}
# 读取所有栅格文件
# envir_stack <- rast(tif_files)

# 查看加载的变量信息
print("已加载的环境变量：")
print(names(envir_stack))
print(paste("共加载", nlyr(envir_stack), "个图层"))

# 绘制所有环境变量（分多个图显示）
par(mfrow = c(5, 5))  # 设置5行5列的绘图布局
plot(envir_stack)      # 绘制所有图层
par(mfrow = c(1, 1))   # 恢复默认绘图布局

# 将环境变量赋值给myExpl_or
myExpl_or <- envir_stack

# 检查最终变量
print("最终环境变量数据集结构：")
print(myExpl_or)


# Load species occurrences 
DataSpecies <- read.csv("./input/panda_point.csv")
head(DataSpecies)

# Get corresponding XY coordinates

myRespName  <- "Ailuropoda_melanoleuca"
myRespXY <- DataSpecies[DataSpecies$scientificName == myRespName, c('decimalLongitude', 'decimalLatitude')]

head(myRespXY)

# ca <- terra::vect("./panda_range/panda_range_contry.shp")

# How are the points distributed across our study area?
plot(myExpl_or[[1]], main = "Buffered minimum convex polygon")
# plot(ca, add = TRUE)
points(myRespXY, pch = 19, cex = 0.5, col = "red")

# new predictor layers 
#myExpl_new <- terra::crop(myExpl_or, ca)
#myExpl_new <-  terra::mask(myExpl_new, ca) # 正方形研究范围可能好看点

myExpl_new <- myExpl_or  # 直接使用，不再 crop/mask

plot(myExpl_new[[1]])
# plot(ca, add = TRUE)
points(myRespXY, pch = 19, cex = 0.5, col = "red")

## # 填补缺失值，要不然总是报错Error in mm %*% object$betas : non-conformable arguments
## # 检查每个图层的格子数量
# cell_counts <- sapply(names(myExpl), function(name) ncell(myExpl[[name]]))
# print(cell_counts)
## 检查了每个图层的格子数量是一致的，但是缺失值格子数量不一致，所以将NA值进行填补
# # 填补缺失值，线性插值
# myExpl_filled <- terra::approximate(myExpl_new, method = "linear")
# 
# # 归一化处理（将值缩放到0到1之间）
# myExpl_new <- (myExpl_filled - min(myExpl_filled)) / (max(myExpl_filled) - min(myExpl_filled))
# cat("缺失值填补，归一化完成")
# summary(myExpl_new)

# 填补缺失值，线性插值
myExpl_filled <- terra::approximate(myExpl_new, method = "linear")

# 方法1：分图层单独归一化（推荐）----------------------------------
# 定义归一化函数
norm_layer <- function(r) {
  vals <- values(r)
  min_val <- min(vals, na.rm = TRUE)
  max_val <- max(vals, na.rm = TRUE)
  return((r - min_val)/(max_val - min_val))
}

# 应用归一化
myExpl_norm <- lapply(myExpl_filled, norm_layer) %>% rast()

# 保持变量名一致
names(myExpl_norm) <- names(myExpl_filled)

# 替换原来的归一化结果
myExpl_new <- myExpl_norm
# ----------------------------------------------------------------

cat("缺失值填补，归一化完成")
summary(myExpl_new)

# Geographical filtering occurrences
filt_geo <- flexsdm::occfilt_geo(
  data = myRespXY,
  x = "decimalLongitude",
  y = "decimalLatitude",
  env_layer = myExpl_new,
  method = c('cellsize', factor = '1'),
  prj = crs(myExpl_or)
)

par(mfrow = c(1, 1))
myExpl_new[[1]] %>% plot(main = "Original occurrence data")
terra::points(myRespXY %>% 
                dplyr::select(decimalLongitude, decimalLatitude), col = "blue", cex = 1)
terra::points(filt_geo %>% 
                dplyr::select(decimalLongitude, decimalLatitude), col = "red",cex = 0.5, pch = 17)

myRespXY = filt_geo

# pearson correlations
terra::pairs(myExpl_new)

# 这里hbh，进行修改，注释掉过滤信息，加载所有环境变量。不行，还是要共线性计算过滤
# dealing with multicollinearity
vif_var <- flexsdm::correct_colinvar(myExpl_new, method = c("vif", th = "10"))
vif_var$env_layer
vif_var$removed_variables
vif_var$vif_table
myExpl <- vif_var$env_layer # also can be calculated using: my_Expl <- myExpl_new[[names(vif_var$env_layer)]]



# # 不使用VIF过滤，直接使用所有环境变量
# myExpl <- myExpl_new
# print("使用所有环境变量，不进行VIF过滤")
# print(names(myExpl))


# Get corresponding presence/absence data
myResp.PA  <- rep(1,dim(myRespXY)[1])

# if(length(myResp.PA) >= 10){
  
################################################################################
# Format Data with pseudo-absences : random method
myBiomodData.r <- BIOMOD_FormatingData(resp.var = myResp.PA,
                                       expl.var = myExpl,
                                       resp.xy = myRespXY,
                                       resp.name = myRespName,
                                       PA.nb.rep = 1,
                                       PA.nb.absences = 10000,
                                       PA.strategy = 'random')

myBiomodData.r
plot(myBiomodData.r)
sum(is.na(myBiomodData.r@data.env.var))

# Print default modeling options
bm_DefaultModelingOptions()

# Create default modeling options
# myBiomodOptions <- BIOMOD_ModelingOptions(MAXENT = list(path_to_maxent.jar =file.path("/Share/user/fanhzh/panda_sdm_new/maxent.jar")))
myBiomodOptions <- BIOMOD_ModelingOptions(
  MAXENT = list(path_to_maxent.jar = "./maxent.jar")
)


myBiomodOptions

# 定义一个函数来处理每个模型的参数
process_model <- function(model_params, model_name) {
  df <- data.frame()
  for (param_name in names(model_params)) {
    param_value <- model_params[[param_name]]
    if (is.list(param_value)) {
      for (sub_param_name in names(param_value)) {
        sub_param_value <- param_value[[sub_param_name]]
        # 检查子参数值的类型和长度
        if (is.function(sub_param_value)) {
          sub_param_value <- "<function>"
        } else if (is.null(sub_param_value) || length(sub_param_value) == 0) {
          sub_param_value <- NA_character_
        } else {
          sub_param_value <- as.character(sub_param_value)
        }
        new_row <- data.frame(model = model_name, 
                              param = paste0(param_name, "_", sub_param_name), 
                              value = sub_param_value, 
                              stringsAsFactors = FALSE)
        df <- bind_rows(df, new_row)
      }
    } else {
      # 检查参数值的类型和长度
      if (is.function(param_value)) {
        param_value <- "<function>"
      } else if (is.null(param_value) || length(param_value) == 0) {
        param_value <- NA_character_
      } else {
        param_value <- as.character(param_value)
      }
      new_row <- data.frame(model = model_name, 
                            param = param_name, 
                            value = param_value, 
                            stringsAsFactors = FALSE)
      df <- bind_rows(df, new_row)
    }
  }
  return(df)
}

# 使用lapply函数处理myBiomodOptions对象中的每个槽
model_dfs <- lapply(slotNames(myBiomodOptions), function(model_name) {
  model_params <- slot(myBiomodOptions, model_name)
  process_model(model_params, model_name)
})

################################################################################
## Model parameters can also be automatically tuned for your specific
## dataset with an optimization algorithm. The tuning can however
## be quite long. Duration for tuning all models sequentially
## with default optimization settings :
## on 8 x 2.5 GHz processor: approx. 15 min tuning all models

library(doParallel)
cl <- makeCluster(10)
doParallel::registerDoParallel(cl)

time.seq <- system.time(
  bm.tuning <- BIOMOD_Tuning(bm.format = myBiomodData.r, 
                             ME.env = myExpl, 
                             ME.n.bg = ncell(myExpl),
                             models = c("ANN", "MAXNET")
                             )) 

stopCluster(cl)

# plot(bm.tuning$tune.CTA.rpart)
# plot(bm.tuning$tune.CTA.rpart2)
# plot(bm.tuning$tune.RF)
# plot(bm.tuning$tune.ANN)
# plot(bm.tuning$tune.MARS)
# plot(bm.tuning$tune.FDA)
# plot(bm.tuning$tune.GBM)
# plot(bm.tuning$tune.GAM)

# Get tuned modeling options
# save all the tunning results
if(!dir.exists(paste0("./ESM_output/",myRespName))){dir.create(paste0("./ESM_output/",myRespName))}
saveRDS(bm.tuning$models.options,paste0("./ESM_output/",myRespName,"/",myRespName,"_model_tunning.rds"))
parameter_setting <- readRDS(paste0("./ESM_output/",myRespName,"/",myRespName,"_model_tunning.rds"))
myBiomodOptions <- parameter_setting

# 使用lapply函数处理myBiomodOptions对象中的每个槽
model_dfs_tuning <- lapply(slotNames(myBiomodOptions), function(model_name) {
  model_params <- slot(myBiomodOptions, model_name)
  process_model(model_params, model_name)
})


# 保存算法调优后的参数
for (i in 1: length(model_dfs)) {
  names(model_dfs[[i]]) = c("Model","Perset_parameter","Perset_value")
  names(model_dfs_tuning[[i]]) = c("Model","Tuned_parameter","Tuned_value")
  dat <-
    tryCatch({
      # 尝试绑定数据框
      cbind(model_dfs[[i]], model_dfs_tuning[[i]])
    }, error = function(e) {
      # 如果错误消息包含特定的文本，则仅保留model_dfs_tuning[[i]]的结果
      if (grepl("参数值意味着不同的行数", e$message, fixed = TRUE)) {
        return(model_dfs_tuning[[i]])
      } else {
        return(dat) 
      }
    })
  print(dat)
  fwrite(dat, paste0("./ESM_output/",myRespName,"/",dat$Model[1],"_para_tuning_res.csv"))
}

# I got this error because of an issue with the parallel computing back-end for foreach/dopar. 
# You may have some parallel computing going on in the background that is not getting cleaned up fully between runs. 
# The easiest way I have found to fix it is to call this function:
unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
unregister_dopar()
# 
gc()
################################################################################
# Several cross-validation methods are available and can be selected with the BIOMOD_Modeling function, which calls the bm_CrossValidation function to do so. 
# More examples are presented on the Secundary functions webpage.
# Model single models
# GLM : Generalized Linear Model (glm)
# GAM : Generalized Additive Model (gam, gam or bam)
# GBM : Generalized Boosting Model, or usually called Boosted Regression Trees (gbm)
# CTA : Classification Tree Analysis (rpart)
# ANN : Artificial Neural Network (nnet)
# SRE : Surface Range Envelop or usually called BIOCLIM
# FDA : Flexible Discriminant Analysis (fda)
# MARS : Multiple Adaptive Regression Splines (earth)
# RF : Random Forest (randomForest)
# MAXENT : Maximum Entropy (https://biodiversityinformatics.amnh.org/open_source/maxent/)
# MAXNET : Maximum Entropy (maxnet)
# XGBOOST : eXtreme Gradient Boosting Training (xgboost)
my.ESM <- ecospat.ESM.Modeling(data = myBiomodData.r,                  
                               NbRunEval =5,
                               DataSplit = 80, 
                               models.options = myBiomodOptions,             
                               modeling.id = 'AllModels',                
                               models = c("GLM","GBM", "ANN", "MAXNET"),
                               Prevalence = 0.5,                          
                               weighting.score=c("AUC"),
                               parallel=FALSE)
my.ESM 

## Evaluation and average of simple bivariate models to ESMs
## Bivariate small models were finally ensembled via Somer D: D = 2*(AUC-0.5)
## This is because sometimes, maxent cannot get Boyce values
my.ESM_EF <- ecospat.ESM.EnsembleModeling(my.ESM, weighting.score=c("SomersD"),threshold=0)

### Ensamble of Small Models: Projects Simple Bivariate Models Into Current period
my.ESM_proj_current<-ecospat.ESM.Projection(ESM.modeling.output=my.ESM,
                                            new.env=myExpl,
                                            name.env = "current",
                                            parallel = FALSE)

### Projection of calibrated ESMs into new space or time
my.ESM_EFproj_current <- ecospat.ESM.EnsembleProjection(ESM.prediction.output=my.ESM_proj_current,
                                                        ESM.EnsembleModeling.output=my.ESM_EF)

####save the model
mypath<-paste0("./ESM_output/",myRespName)
if(!file.exists(mypath)) dir.create(mypath)
save(my.ESM,my.ESM_EF,my.ESM_proj_current,my.ESM_EFproj_current,file=paste0(mypath,"/",myRespName,"_ESM.Rdata"))

## get the model performance of ESMs
my.ESM_EF$ESM.evaluations
write.csv(my.ESM_EF$ESM.evaluations,paste0("./ESM_output/",myRespName,"/","ESM_evaluation.csv"))

## get the weights of the single bivariate models used to build the ESMs
my.ESM_EF$weights
write.csv(my.ESM_EF$weights,paste0("./ESM_output/",myRespName,"/","ESM_weight.csv"))

## get the variable contributions of ESMs
ecospat.ESM.VarContrib(my.ESM,my.ESM_EF)
write.csv(ecospat.ESM.VarContrib(my.ESM,my.ESM_EF),paste0("./ESM_output/",myRespName,"/","ESM_varContrib.csv"))

## Evaluation of the ensemble models based on the pooling procedure 
my.ESM_evaluations <- ecospat.ESM.EnsembleEvaluation(my.ESM,
                               my.ESM_EF,
                               metrics = c("SomersD","AUC","MaxTSS","MaxKappa","Boyce"), 
                               EachSmallModels = FALSE)
my.ESM_evaluations
write.csv(my.ESM_evaluations$ESM.evaluations,paste0("./ESM_output/",myRespName,"/","ESM_Ensemble_evaluation_pooling.csv"))

## thresholds to produce binary maps
# Modified function for calculating the correct matrix
ecospat.ESM.threshold_mod1 <- function(ESM.EnsembleModeling.output, PEplot=FALSE){
  
  Full.models <- grep('Full',colnames(ESM.EnsembleModeling.output$ESM.fit),value = TRUE)
  
  EVAL <- NULL
  for(i in Full.models){
    DATA <- cbind(1:nrow(ESM.EnsembleModeling.output$ESM.fit),
                  resp.var = ESM.EnsembleModeling.output$ESM.fit$resp.var, 
                  ESM.EnsembleModeling.output$ESM.fit[,i]/1000)
    
    
    EVAL1 <- PresenceAbsence::presence.absence.accuracy(DATA[, ], threshold = as.vector(PresenceAbsence::optimal.thresholds(DATA[, ], opt.methods = "MaxSens+Spec")[-1], mode = "numeric"))
    if(length(EVAL1) > 1){
    TSS.th = EVAL1$threshold
    EVAL1 <- EVAL1[c(1, 4:7, 9:12)]
    EVAL1$SomersD <- EVAL1$AUC * 2 - 1
    boyce <- ecospat.boyce(DATA[,3], DATA[DATA[, 2] == 1,3],PEplot=PEplot)
    EVAL1$Boyce <- boyce$cor
    EVAL1$TSS <- EVAL1$sensitivity + EVAL1$specificity - 1
    EVAL1$TSS.th <- TSS.th
    EVAL1$Boyce.th.max <- EVAL1$Boyce.th.min <- EVAL1$MPA0.90 <-EVAL1$MPA0.95 <-EVAL1$MPA1.0 <- NA
    
    EVAL1$MPA1.0 <- ecospat.mpa(DATA[DATA[, 2] == 1,],perc=1)
    EVAL1$MPA0.95 <- ecospat.mpa(DATA[DATA[, 2] == 1,],perc=.95)
    EVAL1$MPA0.90 <- ecospat.mpa(DATA[DATA[, 2] == 1,],perc=.9)
    
    pos.F <- which(boyce$F.ratio>1)
    neg.F <- which(boyce$F.ratio<=1)
    if(max(neg.F) < min(pos.F)){
      EVAL1$Boyce.th.max <- EVAL1$Boyce.th.min <- mean(boyce$HS[c(max(neg.F),min(pos.F))])
      
    }else{
      EVAL1$Boyce.th.max <- mean(boyce$HS[c(max(neg.F),max(neg.F)+1)])
      EVAL1$Boyce.th.min <- mean(boyce$HS[c(min(pos.F),min(pos.F)-1)])
    }
    EVAL1$model <- i
    EVAL <- rbind(EVAL,EVAL1)
    } else {
      EVAL1 <- data.frame("model"=NA,   
                          "sensitivity"= NA,
                          "specificity"= NA,
                          "Kappa" = NA,
                          "AUC" = NA,
                          "sensitivity.sd" =NA,
                          "specificity.sd" = NA,
                          "Kappa.sd" =NA,
                          "AUC.sd" = NA,
                          "SomersD" = NA,
                          "Boyce" = NA,
                          "TSS" = NA,
                          "TSS.th" = NA,
                          "MPA1.0" = NA,
                          "MPA0.95" = NA,
                          "MPA0.90" = NA,
                          "Boyce.th.min" = NA,  
                          "Boyce.th.max" = NA)
      EVAL1$model <- i 
      EVAL <- rbind(EVAL,EVAL1)
    }
  }
  return(EVAL)
}
# 
my.ESM_thresholds <- ecospat.ESM.threshold_mod1(my.ESM_EF)

## Response curve
pdf(file=paste0("./ESM_output/",myRespName,"/","ESM_response.PDF"))
output.plot <- ecospat.ESM.responsePlot(ESM.EnsembleModeling.output = my.ESM_EF,
                                      ESM.modeling.output = my.ESM)
saveRDS(output.plot,paste0("./ESM_output/",myRespName,"/","all_Raw_response_curve.rds"))
dev.off()

## thresholds for the ensemble small models 
# Modified function for calculating the correct matrix
ecospat.ESM.threshold_mod_TH <- function(DATA, PEplot=FALSE){
    EVAL <- NULL
    EVAL1 <- PresenceAbsence::presence.absence.accuracy(DATA[, ], threshold = as.vector(PresenceAbsence::optimal.thresholds(DATA[, ], opt.methods = "MaxSens+Spec")[-1], mode = "numeric"))
    TSS.th <- EVAL1$threshold
    EVAL1 <- EVAL1[c(1, 4:7, 9:12)]
    EVAL1$SomersD <- EVAL1$AUC * 2 - 1
    boyce <- ecospat.boyce(DATA[,3], DATA[DATA[, 2] == 1,3],PEplot=PEplot)
    EVAL1$Boyce <- boyce$cor
    EVAL1$TSS <- EVAL1$sensitivity + EVAL1$specificity - 1
    EVAL1$TSS.th <- TSS.th
    EVAL1$Boyce.th.max <- EVAL1$Boyce.th.min <- EVAL1$MPA0.90 <-EVAL1$MPA0.95 <-EVAL1$MPA1.0 <- NA
    
    EVAL1$MPA1.0 <- ecospat.mpa(DATA[DATA[, 2] == 1,],perc=1)
    EVAL1$MPA0.95 <- ecospat.mpa(DATA[DATA[, 2] == 1,],perc=.95)
    EVAL1$MPA0.90 <- ecospat.mpa(DATA[DATA[, 2] == 1,],perc=.9)
    
    
    
    pos.F <- which(boyce$F.ratio>1)
    neg.F <- which(boyce$F.ratio<=1)
    if(max(neg.F) < min(pos.F)){
      EVAL1$Boyce.th.max <- EVAL1$Boyce.th.min <- mean(boyce$HS[c(max(neg.F),min(pos.F))])
      
    }else{
      EVAL1$Boyce.th.max <- mean(boyce$HS[c(max(neg.F),max(neg.F)+1)])
      EVAL1$Boyce.th.min <- mean(boyce$HS[c(min(pos.F),min(pos.F)-1)])
    }
    EVAL1$model <- "ESM"
    EVAL <- rbind(EVAL,EVAL1)
    return(EVAL)
}
# 

DATA <- cbind(1:dim(myBiomodData.r@coord)[1],
              resp.var = my.ESM_EF$ESM.fit$resp.var, 
              terra::extract(my.ESM_EFproj_current$EF,myBiomodData.r@coord)[,2]
              /1000)

Final_TH <- rbind(my.ESM_thresholds, ecospat.ESM.threshold_mod_TH(DATA, PEplot=FALSE))
write.csv(Final_TH,paste0("./ESM_output/",myRespName,"/","ESM_TH.csv"))

# predictions for present
writeRaster(my.ESM_EFproj_current,paste0("./ESM_output/",myRespName,"/",myRespName,"_present_continuous.tiff"),overwrite=T)

par(mfrow = c(2, 2))
plot(my.ESM_EFproj_current)
plot(my.ESM_EFproj_current[[1:3]] > c(my.ESM_thresholds$TSS.th[1]*1000,
                                      my.ESM_thresholds$TSS.th[2]*1000,
                                      my.ESM_thresholds$TSS.th[3]*1000))
plot(my.ESM_EFproj_current$EF > Final_TH[which(Final_TH$model == "ESM"),"TSS.th"]*1000, main = "Ensemble results")
terra::points(filt_geo %>% 
                dplyr::select(decimalLongitude, decimalLatitude), col = "red",cex = 1, pch = 17)
plot(my.ESM_EFproj_current$EF)
terra::points(filt_geo %>% 
                dplyr::select(decimalLongitude, decimalLatitude), col = "red",cex = 1, pch = 17)
