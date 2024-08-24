#setwd('D:/R/zymogram') # set current working directory
rm(list=ls()) # clear the environment
cat("\014") # clear the console

# Loading packages
if (!require(pacman)) {install.packages('pacman')}
pacman::p_load(mixtools)

# 导入灰度值数据
Gray <- read.csv('GrayValue.csv', header = FALSE)
Gray <- as.data.frame(lapply(Gray[-c(1, 2), ], as.integer))


# 自定义函数，进行混合高斯模型拟合，并计算阈值、热点比例
process_data <- function(data, col1, col2) {
  if (all(is.na(data[[col1]])) || all(is.na(data[[col2]]))) {
    result_df <- data.frame(
      Proportion_Percent = NA,
      Threshold = NA,
      Comp1_Mean = NA,
      Comp1_SD = NA,
      Comp1_Weight = NA,
      Comp2_Mean = NA,
      Comp2_SD = NA,
      Comp2_Weight = NA
    )
  } else {
    vec <- rep(data[[col1]], times = data[[col2]])
    set.seed(123)
    mod <- normalmixEM(vec)
    
    index_highest_weight <- which.max(mod$lambda)
    mu_background <- mod$mu[index_highest_weight]
    sigma_background <- mod$sigma[index_highest_weight]
    threshold <- mu_background + 2 * sigma_background
    values_above_threshold <- sum(data[[col2]][which(data[[col1]] > threshold)])
    proportion_above_threshold <- values_above_threshold / sum(data[[col2]]) * 100
    
    comp1 <- c(mu = mod$mu[index_highest_weight], sigma = mod$sigma[index_highest_weight], lambda = mod$lambda[index_highest_weight])
    comp2 <- c(mu = mod$mu[-index_highest_weight], sigma = mod$sigma[-index_highest_weight], lambda = mod$lambda[-index_highest_weight])
    
    result_df <- data.frame(
      Proportion_Percent = proportion_above_threshold,
      Threshold = threshold,
      Comp1_Mean = comp1["mu"],
      Comp1_SD = comp1["sigma"],
      Comp1_Weight = comp1["lambda"],
      Comp2_Mean = comp2["mu"],
      Comp2_SD = comp2["sigma"],
      Comp2_Weight = comp2["lambda"]
    )
  }
  return(result_df)
}

combined_df <- data.frame(stringsAsFactors = F)

# 批量计算1:n所有样品的混合高斯分布、阈值和热点比例
for (i in 1:4) {
  col1 <- paste0("V", i*2 - 1)
  col2 <- paste0("V", i*2)
  result_df <- process_data(Gray, col1, col2)
  combined_df <- rbind(combined_df, result_df)
}

rownames(combined_df) <- NULL
write.csv(combined_df, file = "Hotspots.csv", row.names = T, na = "") # 导出结果
