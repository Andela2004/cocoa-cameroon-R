# =====================================================================
# Statistical analysis of cocoa production in Cameroon (2000-2024)
# Author: Ezechiel Marrion Andela Amougou
# Run from the project root:  Rscript scripts/analysis.R
# Input : data/cocoa_cameroon_2000_2024.csv
#         columns: year, cocoa_production_t, harvested_area_ha,
#                  temperature_c, precipitation_mm
# Output: figures/*.png and outputs/model_results.txt
# =====================================================================

library(nlme)       # gls(), corARMA()
library(lmtest)     # bptest()
library(sandwich)   # NeweyWest()

dir.create("figures", showWarnings = FALSE)
dir.create("outputs", showWarnings = FALSE)

# ---- 1. Data ---------------------------------------------------------
d <- read.csv("data/cocoa_cameroon_2000_2024.csv")
stopifnot(all(c("year", "cocoa_production_t", "harvested_area_ha",
                "temperature_c", "precipitation_mm") %in% names(d)))
d <- d[order(d$year), ]
cat("Observations:", nrow(d), " Missing values:", sum(is.na(d)), "\n")

d$yield_t_ha <- d$cocoa_production_t / d$harvested_area_ha

# ---- 2. Exploratory figures -----------------------------------------
png("figures/fig1_production.png", 900, 550, res = 120)
plot(d$year, d$cocoa_production_t, type = "b", pch = 16,
     xlab = "Year", ylab = "Cocoa production (t)",
     main = "Annual cocoa production, Cameroon")
dev.off()

png("figures/fig2_temperature.png", 900, 550, res = 120)
plot(d$year, d$temperature_c, type = "b", pch = 16,
     xlab = "Year", ylab = "Mean temperature (deg C)",
     main = "Annual mean temperature, Cameroon")
abline(lm(temperature_c ~ year, d), lty = 2)
dev.off()

png("figures/fig3_precipitation.png", 900, 550, res = 120)
plot(d$year, d$precipitation_mm, type = "b", pch = 16,
     xlab = "Year", ylab = "Precipitation (mm)",
     main = "Annual precipitation, Cameroon")
abline(lm(precipitation_mm ~ year, d), lty = 2)
dev.off()

# ---- 3. Yield and temperature trend ---------------------------------
cat("\nYield (t/ha): mean =", round(mean(d$yield_t_ha), 3),
    " min =", round(min(d$yield_t_ha), 3),
    " max =", round(max(d$yield_t_ha), 3), "\n")

m_temp <- lm(temperature_c ~ year, data = d)
nw <- NeweyWest(m_temp, lag = 2, prewhite = FALSE)
b  <- coef(m_temp)["year"]; se <- sqrt(diag(nw))["year"]
tcrit <- qt(0.975, df.residual(m_temp))
cat("Temperature trend (Newey-West SE, lag 2): ",
    round(b, 4), " C/year; 95% CI [",
    round(b - tcrit * se, 4), ",", round(b + tcrit * se, 4), "]; p =",
    signif(2 * pt(-abs(b / se), df.residual(m_temp)), 3), "\n")

# ---- 4. Common-trend check ------------------------------------------
raw <- cor.test(d$cocoa_production_t, d$temperature_c)
dt  <- cor.test(resid(lm(cocoa_production_t ~ year, d)),
                resid(lm(temperature_c ~ year, d)))
cat("\nRaw correlation production-temperature: r =", round(raw$estimate, 3),
    " p =", signif(raw$p.value, 3), "\n")
cat("After removing linear time trends:      r =", round(dt$estimate, 3),
    " p =", signif(dt$p.value, 3), "\n")

# ---- 5. First-differenced data --------------------------------------
n <- nrow(d)
dd <- data.frame(
  year   = d$year[-1],
  dlogP  = diff(log(d$cocoa_production_t)),
  dTemp  = diff(d$temperature_c),
  dlogPr = diff(log(d$precipitation_mm)),
  dlogA  = diff(log(d$harvested_area_ha))
)

f <- dlogP ~ dTemp + dlogPr + dlogA

# ---- 6. GLS with AR(1) and AR(2) errors (maximum likelihood) --------
ctrl <- glsControl(maxIter = 200, msMaxIter = 200)
m_ar1 <- gls(f, data = dd, correlation = corARMA(p = 1, q = 0, form = ~ year),
             method = "ML", control = ctrl)
m_ar2 <- gls(f, data = dd, correlation = corARMA(p = 2, q = 0, form = ~ year),
             method = "ML", control = ctrl)

cat("\n=== AIC comparison (ML) ===\n")
print(AIC(m_ar1, m_ar2))
cat("Delta AIC (AR2 - AR1):", round(AIC(m_ar2) - AIC(m_ar1), 3), "\n")

# Retained model: the one with the lower AIC (AR(1) kept if the gap is < 2)
final <- if (AIC(m_ar2) < AIC(m_ar1) - 2) m_ar2 else m_ar1
p_ar  <- if (identical(final, m_ar2)) 2 else 1
cat("Retained model: AR(", p_ar, ")\n", sep = "")

# ---- 7. Coefficients with 95% t-based confidence intervals ----------
tt  <- summary(final)$tTable
dfr <- final$dims$N - final$dims$p
tc  <- qt(0.975, dfr)
coef_tab <- data.frame(
  term     = rownames(tt),
  estimate = tt[, "Value"],
  ci_low   = tt[, "Value"] - tc * tt[, "Std.Error"],
  ci_high  = tt[, "Value"] + tc * tt[, "Std.Error"],
  p_value  = tt[, "p-value"],
  row.names = NULL
)
cat("\n=== Coefficients (df =", dfr, ") ===\n")
print(format(coef_tab, digits = 4), row.names = FALSE)
cat("\nCorrelation parameters:\n"); print(coef(final$modelStruct$corStruct,
                                               unconstrained = FALSE))

# ---- 8. Residual diagnostics of the FINAL model ---------------------
res <- residuals(final, type = "normalized")
fit <- fitted(final)

for (L in c(5, 10)) {
  lb <- Box.test(res, lag = L, type = "Ljung-Box", fitdf = p_ar)
  cat("Ljung-Box lag", L, ": X-squared =", round(lb$statistic, 3),
      " df =", lb$parameter, " p =", round(lb$p.value, 3), "\n")
}
bp <- bptest(res ~ fit, data = data.frame(res = res, fit = fit))
cat("Breusch-Pagan (vs fitted): p =", round(bp$p.value, 3), "\n")
sw <- shapiro.test(res)
cat("Shapiro-Wilk on residuals: p =", round(sw$p.value, 3), "\n")

png("figures/fig4_resid_vs_fitted.png", 900, 550, res = 120)
plot(fit, res, pch = 16, xlab = "Fitted values",
     ylab = "Normalized residuals", main = "Residuals vs fitted (final model)")
abline(h = 0, lty = 2); dev.off()

png("figures/fig5_qq.png", 900, 550, res = 120)
qqnorm(res, pch = 16, main = "Normal Q-Q plot (final model)")
qqline(res); dev.off()

png("figures/fig6_acf.png", 900, 550, res = 120)
acf(res, main = "ACF of residuals (final model)"); dev.off()

png("figures/fig7_resid_time.png", 900, 550, res = 120)
plot(dd$year, res, type = "b", pch = 16, xlab = "Year",
     ylab = "Normalized residuals", main = "Residuals over time (final model)")
abline(h = 0, lty = 2); dev.off()

# ---- 9. Save text summary -------------------------------------------
sink("outputs/model_results.txt")
cat("AIC:\n"); print(AIC(m_ar1, m_ar2))
cat("\nFinal model summary:\n"); print(summary(final))
cat("\nCoefficients with t-based 95% CI:\n"); print(coef_tab, digits = 4)
sink()
cat("\nDone. Figures in figures/, summary in outputs/model_results.txt\n")
