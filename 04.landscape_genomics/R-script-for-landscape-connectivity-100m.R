################################################################################
# This script used for Landscape connectivity analysis.
# The first part shows how to use ResistanceGA optimizing resistance surfaces
# R 4.3.1
# https://www.nature.com/articles/s41558-023-01772-8

# ———— Bug workaround for R 4.3.1 ————
# 如果 options("digits") 或 options("scipen") 是 NA，就重置为默认值
od <- getOption("digits")
os <- getOption("scipen")
if (is.na(od)) options(digits = 7)    # 默认是 7 :contentReference[oaicite:2]{index=2}
if (is.na(os)) options(scipen = 0)    # 默认是 0 :contentReference[oaicite:3]{index=3}

############
#first part#
############
#install.packages("devtools")
#devtools::install_github("wpeterman/ResistanceGA") # v4.2-10
#install.packages("JuliaCall")
#install.packages("progress")
#library(progress)
library(sp)                                                 # v2.1-1
library(raster)                                             # v3.6-26
library(ResistanceGA)                                       # v4.2-10
library(foreach)                                            # v1.5.2
library(parallel)                                           # v4.3.1
library(iterators)                                          # v1.0.14
library(doParallel)                                         # v1.0.17
library(gdistance)                                          # v1.6.4
library(JuliaCall)                                          # v0.17.6

################################################################################
# all_layer <- raster("./MultipleSurfaces_mantzrum_current_result_comb/all.combrep_1/elevation.habitiescurrent.land.river.roading/elevation.habitiescurrent.land.river.roading.asc")
# land_layer <- raster("./MultipleSurfaces_mantzrum_current_result_comb/all.combrep_1/single.surface/land.asc")
# river_layer <- raster("./MultipleSurfaces_mantzrum_current_result_comb/all.combrep_1/single.surface/river.asc")
# elevation_layer <- raster("./MultipleSurfaces_mantzrum_current_result_comb/all.combrep_1/single.surface/elevation.asc")
# roading_layer <- raster("./MultipleSurfaces_mantzrum_current_result_comb/all.combrep_1/single.surface/roading.asc")
# habitiescurrent_layer <- raster("./MultipleSurfaces_mantzrum_current_result_comb/all.combrep_1/single.surface/habitiescurrent.asc")
# 
# stack_layer <- stack(all_layer,land_layer,river_layer,elevation_layer,roading_layer,habitiescurrent_layer)
# plot(stack_layer)
# 
# model1 <- lm(values(all_layer)~values(land_layer)+values(river_layer)+values(elevation_layer)+values(roading_layer)+values(habitiescurrent_layer))
# summary(model1)


# Create a subdirectory for the second example
dir.create(file.path("giant_panda_300m"))

# Directory to write .asc files and results
write.dir <- paste0(getwd(),"/", "giant_panda_300m")
setwd(write.dir)

################################################################################
#Environmental raster stack
# 导入环境栅格图层
habitatcurrent <- raster("../Resistance_layer/habitat_resistance_100m.tif")
grazint <- raster("../Resistance_layer/Grazing_Intensity_resistance_100m.tif")
elevation <- raster("../Resistance_layer/DEM_resistance_100m.tif")
land <- raster("../Resistance_layer/land_cover_resistance_100m.tif")
slope <- raster("../Resistance_layer/slope_resistance_100m.tif")
all_layer <- stack(grazint,elevation,land,slope)
names(all_layer) <- c("grazing","ele","land","slope")
plot(all_layer)
all_layer_crop <- resample(all_layer, habitatcurrent, method = "bilinear")
all_layers <- stack(all_layer_crop,habitatcurrent)
all_layers_agg <- all_layers  %>% aggregate(fact = 3, fun = mean)
plot(all_layers_agg)

# habitatcurrent <- raster("../Input/monticola_resample.tif") %>% aggregate(fact = 5, fun = mean)
# roading <- raster("../Input/road_resample.tif")%>% aggregate(fact = 5, fun = mean)
# river <- raster("../Input/river_resample.tif")%>% aggregate(fact = 5, fun = mean)
# elevation <- raster("../Input/elev.tif")%>% aggregate(fact = 5,fun = mean)
# land <- raster("../Input/land_resample.tif")%>% aggregate(fact = 5, fun = mean)
# all_layer <- stack(habitatcurrent,roading,river,elevation,land)
# plot(all_layer)



# 检查分辨率
#res_list <- list(res(habitat), res(elevation), res(roading), res(river), res(land))
#print(res_list)

# 检查投影
#crs_list <- list(crs(habitat), crs(elevation), crs(roading_density), crs(water_distance), crs(land_use))
#print(crs_list)

# 检查边界
#ext_list <- list(extent(habitat), extent(elevation), extent(roading_density), extent(water_distance), extent(land_use))
#print(ext_list)

# 转换为 ASC 格式
writeRaster(all_layers_agg[[1]], paste0(write.dir,"/grazing.asc"), format = "ascii", overwrite = TRUE)
writeRaster(all_layers_agg[[2]], paste0(write.dir,"/ele.asc"), format = "ascii", overwrite = TRUE)
writeRaster(all_layers_agg[[3]], paste0(write.dir,"/land.asc"), format = "ascii", overwrite = TRUE)
writeRaster(all_layers_agg[[4]], paste0(write.dir,"/slope.asc"), format = "ascii", overwrite = TRUE)
writeRaster(all_layers_agg[[5]], paste0(write.dir,"/habitat.asc"), format = "ascii", overwrite = TRUE)
################################################################################

# 导入坐标数据
# 读取样本数据
xy_MN <- read.csv("../populations_locations.csv", header = TRUE)

# 检查经纬度列的名称是否正确
# 假设经纬度列名为 lon 和 lat
if (!all(c("lon", "lat") %in% names(xy_MN))) {
  stop("Error: The columns 'lon' and 'lat' must be present in the CSV file.")
}

# 创建 SpatialPoints 对象
coordinates(xy_MN) <- ~ lon + lat

# 设置投影
proj4string(xy_MN) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

# 将其转换为 SpatialPoints 对象
xy_MN_points <- as(xy_MN, "SpatialPoints")


# 导入遗传距离矩阵
gen_mat <- read.csv("../Genetic_distance.csv", header = F)

#准备外接交互
Sys.time()
JULIA_HOME <- "F:/Program Files/Julia-1.11.5/bin"
# JULIA_HOME <- "C:/Users/admin/AppData/Local/Programs/Julia-1.11.5/bin"
JuliaCall::julia_setup(JULIA_HOME)


# 1 = "Inverse-Reverse Monomolecular"
# 2 = "Inverse-Reverse Ricker"
# 3 = "Monomolecular"
# 4 = "Ricker"
# 5 = "Reverse Monomolecular"
# 6 = "Reverse Ricker"
# 7 = "Inverse Monomolecular"
# 8 = "Inverse Ricker"
# 9 = "Distance"

# 准备优化输入
# 20 线程 4.625344 mins
# $layer.names
# [1] "ele"     "grazing" "habitat" "land"    "slope"  
t1 <- Sys.time()
GA.inputs <- GA.prep(ASCII.dir = write.dir,  # 指定ASC文件的路径
                     method = "LL",
                     select.trans = list(c(1,3), c(1,3), c(5,7),c(1,3),c(1,3)),
                     max.cat = 1000,
                     max.cont = 1000,
                     run = 10,
                     seed = 555,
                     parallel = 5)

# 2) Quick test setup to find unreachable pairs
jl_test <- jl.prep(
  n.Pops          = length(xy_MN_points),
  response        = gen_mat[lower.tri(gen_mat)],
  CS_Point.File   = xy_MN_points,
  JULIA_HOME      = JULIA_HOME,
  run_test        = FALSE
)

flat.parm <- rep(1, sum(GA.inputs$parm.type$n.parm))

test.RAST <- Combine_Surfaces(
  PARM      = flat.parm,
  jl.inputs = jl_test,
  GA.inputs = GA.inputs,
  rescale   = TRUE
)

test.cs <- Run_CS.jl(jl_test, test.RAST, full.mat = TRUE)


# 3) Build full‐length response & mask
all_response <- gen_mat[lower.tri(gen_mat)]
keep_mask    <- as.integer(lower(test.cs) != -1)

# 4) Correct jl.prep() call using mask
jl.inputs <- jl.prep(
  n.Pops           = length(xy_MN_points),
  response         = all_response,       # full vector of length 78
  CS_Point.File    = xy_MN_points,
  pairs_to_include = keep_mask,          # same length 78
  JULIA_HOME       = JULIA_HOME,
  run_test         = FALSE
)

# 5) Run the optimizations
comb.optim <- MS_optim(jl.inputs = jl.inputs, GA.inputs = GA.inputs)
# comb.optim <- SS_optim(jl.inputs = jl.inputs, GA.inputs = GA.inputs)


t2 <- Sys.time()
t_used <- t2 - t1
print(t_used)

results <- raster("./Results/ele.grazing.habitat.land.slope.asc")
plot(results)

# ==== 新增部分：提取每个阻力图层的转换参数 ==== #

# 查看综合优化结果汇总
summary_results <- comb.optim@Summary
print(summary_results)

# 提取每个变量的最优转换信息
all_surfaces <- jl.ms.optim@MS_results

# 输出每个阻力面的转换函数类型及其参数
for (i in seq_along(all_surfaces)) {
  cat("\n======== 表层", i, "========\n")
  cat("原始图层名: ", names(all_surfaces)[i], "\n")
  cat("变换类型: ", all_surfaces[[i]]$Transformation, "\n")
  cat("参数1: ", all_surfaces[[i]]$shape, "\n")
  cat("参数2: ", all_surfaces[[i]]$max, "\n")  # 有些变换有max参数
}

# ==== 综合阻力图输出 ==== #

# 综合阻力面栅格（加权求和结果）
final_resistance_surface <- jl.ms.optim@Resistance

# 保存最终综合阻力图
writeRaster(final_resistance_surface,
            filename = file.path(write.dir, "composite_resistance_surface.tif"),
            format = "GTiff", overwrite = TRUE)
################################################################################

