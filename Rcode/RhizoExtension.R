#setwd('D:/R/zymogram') # set current working directory
rm(list=ls()) # clear the environment
cat("\014") # clear the console

# 载入函数包
if (!require(pacman)) {install.packages('pacman')}
pacman::p_load(minpack.lm, ggplot2)

# 导入灰度值数据
Gray <- read.csv('PixelGray.csv', header = FALSE)
Gray <- as.data.frame(lapply(Gray[-c(1, 2), ], as.numeric))


# 自定义函数，进行曲线拟合，输出关键参数x0和带拟合线的散点图
process_data <- function(data, col1, col2) {
  plot_list <- list()
  if (all(is.na(data[[col1]])) || all(is.na(data[[col2]]))) {
    x0_result <- data.frame(x0 = NA)
    plot_list[[paste0(i)]] <- list()
  } else {
    x <- unlist(na.omit(Gray[col1]))
    y <- unlist(na.omit(Gray[col2]))
    origin_data <- data.frame(x, y)
    
    # 定义四参数logistic函数用于拟合根际灰度值/活性散点
    sigmoid <- function(x, a, b, x0, y0) {y0 + a / (1 + exp(-(x - x0) / b))}
    fit <- nlsLM(y ~ sigmoid(x, a, b, x0, y0), data = origin_data,
                 start = c(a = 0, b = -1, x0 = 0, y0 = 0), control = nls.lm.control(maxiter = 1000, maxfev = 10000))
    rsquared <- 1 - sum(residuals(fit)^2) / sum((y - mean(y))^2)
    full_model <- nlsLM(y ~ sigmoid(x, a, b, x0, y0), start = coef(fit), control = nls.lm.control(maxiter = 1000, maxfev = 10000))
    p_value_model <- 1 - pchisq((AIC(fit) - AIC(full_model)), df = 1)
    
    fit_a <- coef(fit)["a"]
    fit_b <- coef(fit)["b"]
    fit_x0 <- coef(fit)["x0"]
    fit_y0 <- coef(fit)["y0"]
    if (fit_a < 0) {fit_y0 <- fit_y0 + fit_a}
    if (fit_a < 0) {fit_a <- abs(fit_a)}
    if (fit_b > 0) {fit_b <- -abs(fit_b)}
    if (fit_x0 < 0) {fit <- nlsLM(y ~ sigmoid(x, a, b, x0, y0), data = origin_data,
                                  start = c(a = 0, b = -1, x0 = 0, y0 = 0),
                                  control = nls.lm.control(maxiter = 1000, maxfev = 10000),
                                  lower = c(a = 0, b = -Inf, x0 = 0, y0 = 0),
                                  upper = c(a = 1.1 * (max(y) - min(y)), b = 0, x0 = Inf, y0 = Inf))
    rsquared <- 1 - sum(residuals(fit)^2) / sum((y - mean(y))^2)
    full_model <- nlsLM(y ~ sigmoid(x, a, b, x0, y0), start = coef(fit), control = nls.lm.control(maxiter = 1000, maxfev = 10000),
                        lower = c(a = 0, b = -Inf, x0 = 0, y0 = 0),
                        upper = c(a = 1.1 * (max(y) - min(y)), b = 0, x0 = Inf, y0 = Inf))
    p_value_model <- 1 - pchisq((AIC(fit) - AIC(full_model)), df = 1)
    
    fit_a <- coef(fit)["a"]
    fit_b <- coef(fit)["b"]
    fit_x0 <- coef(fit)["x0"]
    fit_y0 <- coef(fit)["y0"]}
    
    x0_result <- data.frame(x0 = fit_x0) # x0表示根际延伸量（根际范围）
    
    curve_data <- data.frame(x = seq(min(x) - 1 , max(x) + 1, length.out = 100),
                             y = sigmoid(seq(min(x) - 1 , max(x) + 1, length.out = 100), fit_a, fit_b, fit_x0, fit_y0))
    x0_point <- data.frame(x = fit_x0, y = predict(fit, data.frame(x = fit_x0)))
    
    plot_list[[paste0(i)]] <- ggplot(origin_data, aes(x, y)) +
      geom_point() +
      geom_line(data = curve_data, aes(x, y), color = "black") +
      scale_x_continuous(limits = c(min(x) - 1, max(x) + 1), expand = c(0, 0)) +
      scale_y_continuous(limits = c(fit_y0 - 0.1 * fit_a, fit_y0 + fit_a + 0.1 * fit_a), expand = c(0, 0)) +
      geom_text(aes(label = paste("y =", round(fit_y0, 2), "+", round(fit_a, 2), "/(1 + exp(-(x -", round(fit_x0, 2), ")/", round(fit_b, 2), "))"), 
                    x = max(x) + 0.9, y = fit_y0 + fit_a, hjust = 1, vjust = -0.5)) +
      geom_text(aes(label = paste("R² =", round(rsquared, 4), ", P =", round(p_value_model, 4)),
                    x = max(x) + 0.9, y = fit_y0 + fit_a, hjust = 1, vjust = 1)) +
      geom_hline(yintercept = fit_y0, linetype = "dashed", color = "orange") +
      geom_hline(yintercept = fit_y0 + fit_a, linetype = "dashed", color = "blue") +
      geom_point(data = x0_point, aes(x, y), color = "red", size = 3) +
      geom_segment(aes(x = fit_x0, xend = fit_x0, y = fit_y0 - 0.1 * fit_a, yend = predict(fit, x0_point), color = "red")) +
      geom_segment(aes(x = min(x) - 1, xend = fit_x0, y = predict(fit, x0_point), yend = predict(fit, x0_point), color = "red")) +
      labs(title = paste0("Scatter Plot with Fitted Curve", "-p", i), x = "Distance_(pixels)", y = "Gray_Value") +
      theme(legend.position = "none") +
      theme(aspect.ratio = 1) +
      geom_text(aes(label = paste("(", round(fit_x0, 2), ",", round(predict(fit, x0_point), 2), ")"),
                    x = fit_x0, y = predict(fit, x0_point), hjust = 0, vjust = -1), color = "red") +
      geom_text(aes(label = paste("y0 + a =", round(fit_y0 + fit_a, 2)), x = min(x) - 0.9, y = fit_y0 + fit_a, hjust = 0, vjust = -0.5), color = "blue") +
      geom_text(aes(label = paste("y0 =", round(fit_y0, 2)), x = min(x) - 0.9, y = fit_y0, hjust = 0, vjust = -0.5), color = "orange")
  }
  list_x0_plot <- list(x0 = x0_result, plot = plot_list)
  return(list_x0_plot)
}

combined_results <- data.frame(stringsAsFactors = F)
combined_plots <- list()

# 批量运算1:n所有lines的根际范围
for (i in 1:12) {
  col1 <- paste0("V", i*2 - 1)
  col2 <- paste0("V", i*2)
  list_x0_plot <- process_data(Gray, col1, col2)
  combined_results <- rbind(combined_results, list_x0_plot$x0)
  combined_plots <- c(combined_plots, list_x0_plot$plot)
}

rownames(combined_results) <- NULL
write.csv(combined_results, file = "RhizoExtension_Lines.csv", row.names = T, na = "") # 导出所有lines的根际范围

# combined_plots[1]


# 将同一样品的lines的根际范围取平均，以获得样品的根际范围
# 自定义函数，利用四分位距法判定异常值，进而将非异常值取平均
detect_outliers_and_calc_mean <- function(x) {
  q1 <- quantile(x, 0.25, na.rm = T)
  q3 <- quantile(x, 0.75, na.rm = T)
  iqr <- q3 - q1
  lower_bound <- q1 - 1.5 * iqr
  upper_bound <- q3 + 1.5 * iqr
  non_outliers <- x[x >= lower_bound & x <= upper_bound]
  mean_non_outliers <- mean(non_outliers, na.rm = T)
  return(mean_non_outliers)
}

n <- 4 # 规定每n个为一组，批量计算均值
grouped_results <- tapply(combined_results$x0, rep(1:(nrow(combined_results) %/% n), each = n), detect_outliers_and_calc_mean)
grouped_results <- data.frame(x0 = grouped_results)
write.csv(grouped_results, file = "RhizoExtension_Samples.csv", row.names = T, na = "") # 导出所有样品的根际范围
