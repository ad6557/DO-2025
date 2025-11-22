################## when stage 2 enroll 20/arm/month, lambda = 0.06 ###############
N2_list = c(350,370,390,410,430,450,470,490)
dur_list = vector("numeric",8)
IA_list = vector("numeric",8)
FA_list = vector("numeric",8)
er_list = vector("numeric",8)

load("T1S1_N2_350_HR_0.7_0.7.RData")
er_list[1] = max(data_s$dif)
dur_list[1] = duration1
IA_list[1] = mean(output_suc$`IA time`)
FA_list[1] = mean(output_suc$`FA time`)
load("T1S1_N2_370_HR_0.7_0.7.RData")
er_list[2] = max(data_s$dif)
dur_list[2] = duration1
IA_list[2] = mean(output_suc$`IA time`)
FA_list[2] = mean(output_suc$`FA time`)
load("T1S1_N2_390_HR_0.7_0.7.RData")
er_list[3] = max(data_s$dif)
dur_list[3] = duration1
IA_list[3] = mean(output_suc$`IA time`)
FA_list[3] = mean(output_suc$`FA time`)
load("T1S1_N2_410_HR_0.7_0.7.RData")
er_list[4] = max(data_s$dif)
dur_list[4] = duration1
IA_list[4] = mean(output_suc$`IA time`)
FA_list[4] = mean(output_suc$`FA time`)
load("T1S1_N2_430_HR_0.7_0.7.RData")
er_list[5] = max(data_s$dif)
dur_list[5] = duration1
IA_list[5] = mean(output_suc$`IA time`)
FA_list[5] = mean(output_suc$`FA time`)
load("T1S1_N2_450_HR_0.7_0.7.RData")
er_list[6] = max(data_s$dif)
dur_list[6] = duration1
IA_list[6] = mean(output_suc$`IA time`)
FA_list[6] = mean(output_suc$`FA time`)
load("T1S1_N2_470_HR_0.7_0.7.RData")
er_list[7] = max(data_s$dif)
dur_list[7] = duration1
IA_list[7] = mean(output_suc$`IA time`)
FA_list[7] = mean(output_suc$`FA time`)
load("T1S1_N2_490_HR_0.7_0.7.RData")
er_list[8] = max(data_s$dif)
dur_list[8] = duration1
IA_list[8] = mean(output_suc$`IA time`)
FA_list[8] = mean(output_suc$`FA time`)


# Define colors and line types
cols <- c("firebrick", "forestgreen", "royalblue", "goldenrod")
ltys <- c(1, 2, 3, 4)
lwds <- c(2.5, 2.5, 2.5, 2.5)
pch_vals <- c(16, 17, 15, 18)

# Create the plot
plot(N2_list, dur_list, type = "l", col = cols[1], lty = ltys[1], lwd = lwds[1],
     xlim = c(350, 490), ylim = c(25, 45),
     xlab = "Number of People Enrolled per Arm in Stage 2", ylab = "Time (months)", 
     main = "Trial Time of Operational Phase 2/3 Design",
     cex.lab = 1.2, cex.main = 1.3)

# Add grid
grid()
# duration to success of inf ph2/3
abline(h = 33.39478, col = "gray30", lwd = 1.5)
# Annotate the line
text(x = 345, y = 33, labels = "Duration to Success of Inferential Phase 2/3 Design (350/arm in Stage 2)", col = "gray30", cex = 0.8, pos = 4)

# Add other lines
lines(N2_list, FA_list, col = cols[2], lty = ltys[2], lwd = lwds[2])
lines(N2_list, IA_list, col = cols[3], lty = ltys[3], lwd = lwds[3])
lines(N2_list, er_list, col = cols[4], lty = ltys[4], lwd = lwds[4])

points(N2_list, dur_list, col = cols[1], pch = pch_vals[1], cex = 0.8)
points(N2_list, FA_list, col = cols[2], pch = pch_vals[2], cex = 0.8)
points(N2_list, IA_list, col = cols[3], pch = pch_vals[3], cex = 0.8)
points(N2_list, er_list, col = cols[4], pch = pch_vals[4], cex = 0.8)

# Add legend
legend("topleft",
       legend = c("Duration to Success", "Duration to FA", "Duration to IA", "Enrollment End Time"),
       col = cols, lty = ltys, lwd = lwds, pch = pch_vals,
       cex = 0.85,  # Keep font size decent
       inset = 0.005,         # Reduce padding from plot edge
       x.intersp = 0.3,       # Reduce horizontal spacing between symbols and text
       y.intersp = 0.3,       # Reduce vertical spacing between lines
       seg.len = 1.2,         # Shorter line segments
       bg = "white", box.lty = 0)

################## when stage 2 enroll 40/arm/month, lambda = 0.06 ###############
N2_list = c(350,370,390,410,430,450,470,490)
dur_list = vector("numeric",8)
IA_list = vector("numeric",8)
FA_list = vector("numeric",8)
er_list = vector("numeric",8)

load("T1S1_N2_350_HR_0.7_0.7_40.RData")
er_list[1] = max(data_s$dif)
dur_list[1] = duration1
IA_list[1] = mean(output_suc$`IA time`)
FA_list[1] = mean(output_suc$`FA time`)
load("T1S1_N2_370_HR_0.7_0.7_40.RData")
er_list[2] = max(data_s$dif)
dur_list[2] = duration1
IA_list[2] = mean(output_suc$`IA time`)
FA_list[2] = mean(output_suc$`FA time`)
load("T1S1_N2_390_HR_0.7_0.7_40.RData")
er_list[3] = max(data_s$dif)
dur_list[3] = duration1
IA_list[3] = mean(output_suc$`IA time`)
FA_list[3] = mean(output_suc$`FA time`)
load("T1S1_N2_410_HR_0.7_0.7_40.RData")
er_list[4] = max(data_s$dif)
dur_list[4] = duration1
IA_list[4] = mean(output_suc$`IA time`)
FA_list[4] = mean(output_suc$`FA time`)
load("T1S1_N2_430_HR_0.7_0.7_40.RData")
er_list[5] = max(data_s$dif)
dur_list[5] = duration1
IA_list[5] = mean(output_suc$`IA time`)
FA_list[5] = mean(output_suc$`FA time`)
load("T1S1_N2_450_HR_0.7_0.7_40.RData")
er_list[6] = max(data_s$dif)
dur_list[6] = duration1
IA_list[6] = mean(output_suc$`IA time`)
FA_list[6] = mean(output_suc$`FA time`)
load("T1S1_N2_470_HR_0.7_0.7_40.RData")
er_list[7] = max(data_s$dif)
dur_list[7] = duration1
IA_list[7] = mean(output_suc$`IA time`)
FA_list[7] = mean(output_suc$`FA time`)
load("T1S1_N2_490_HR_0.7_0.7_40.RData")
er_list[8] = max(data_s$dif)
dur_list[8] = duration1
IA_list[8] = mean(output_suc$`IA time`)
FA_list[8] = mean(output_suc$`FA time`)


# Define colors and line types
cols <- c("firebrick", "forestgreen", "royalblue", "goldenrod")
ltys <- c(1, 2, 3, 4)
lwds <- c(2.5, 2.5, 2.5, 2.5)
pch_vals <- c(16, 17, 15, 18)

# Create the plot
plot(N2_list, dur_list, type = "l", col = cols[1], lty = ltys[1], lwd = lwds[1],
     xlim = c(350, 490), ylim = c(20, 40),
     xlab = "Number of People Enrolled per Arm in Stage 2", ylab = "Time (months)", 
     main = "Trial Time of Operational Phase 2/3 Design",
     cex.lab = 1.2, cex.main = 1.3)

# Add grid
grid()
# duration to success of inf ph2/3
abline(h = 29.13423, col = "gray30", lwd = 1.5)
# Annotate the line
text(x = 345, y = 28.6, labels = "Duration to Success of Inferential Phase 2/3 Design (350/arm in Stage 2)", col = "gray30", cex = 0.8, pos = 4)

# Add other lines
lines(N2_list, FA_list, col = cols[2], lty = ltys[2], lwd = lwds[2])
lines(N2_list, IA_list, col = cols[3], lty = ltys[3], lwd = lwds[3])
lines(N2_list, er_list, col = cols[4], lty = ltys[4], lwd = lwds[4])

points(N2_list, dur_list, col = cols[1], pch = pch_vals[1], cex = 0.8)
points(N2_list, FA_list, col = cols[2], pch = pch_vals[2], cex = 0.8)
points(N2_list, IA_list, col = cols[3], pch = pch_vals[3], cex = 0.8)
points(N2_list, er_list, col = cols[4], pch = pch_vals[4], cex = 0.8)

# Add legend
legend("topleft",
       legend = c("Duration to Success", "Duration to FA", "Duration to IA", "Enrollment End Time"),
       col = cols, lty = ltys, lwd = lwds, pch = pch_vals,
       cex = 0.85,  # Keep font size decent
       inset = 0.02,         # Reduce padding from plot edge
       x.intersp = 0.3,       # Reduce horizontal spacing between symbols and text
       y.intersp = 0.8,       # Reduce vertical spacing between lines
       seg.len = 2.5,         # Shorter line segments
       bg = "white", box.lty = 0)