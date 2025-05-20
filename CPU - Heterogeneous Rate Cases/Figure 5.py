import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# ------------------ 全局设置 ------------------
plt.rcParams["font.family"] = "Times New Roman"  # 设置字体
plt.rcParams["font.size"] = 14

# ------------------ 数据加载 ------------------
df = pd.read_csv("heterResults.csv")

# ------------------ 数据预处理 ------------------
# 转换长表
df_long = pd.melt(
    df,
    id_vars=["Location", "ρ"],
    value_vars=[f"N={n} ({metric})" for n in range(11, 16) for metric in ["Time Saving", "MPRE"]],
    var_name="Metric",
    value_name="Value"
)

# 提取N值和指标类型
df_long["N"] = df_long["Metric"].apply(lambda x: int(x.split("=")[1].split(" ")[0]))
df_long["Metric"] = df_long["Metric"].apply(lambda x: x.split("(")[1].strip(")"))

# 提取子数据集
st_paul_data = df_long[df_long["Location"] == "St. Paul"].copy()
greenville_data = df_long[df_long["Location"] == "Greenville"].copy()

# ------------------ 创建图形 ------------------
fig, axes = plt.subplots(1, 2, figsize=(16, 6))
plt.subplots_adjust(wspace=0.3)  # 控制子图间距

# ------------------ 定义通用参数 ------------------
box_width = 0.3
position_offset = 0.2
outlier_colors = {"Greenville": "green", "St. Paul": "blue"}
ylim_config = {
    "Time Saving": (40, 102),
    "MPRE": (0, 0.5)
}

# ------------------ 绘制Time Saving对比图 ------------------
# 准备数据
time_saving_st_paul = [st_paul_data[(st_paul_data["Metric"] == "Time Saving") & (st_paul_data["N"] == n)]["Value"].values for n in range(11, 16)]
time_saving_greenville = [greenville_data[(greenville_data["Metric"] == "Time Saving") & (greenville_data["N"] == n)]["Value"].values for n in range(11, 16)]

# 绘制St. Paul箱线图
bp1 = axes[0].boxplot(
    time_saving_st_paul,
    positions=np.arange(11, 16) - position_offset,
    widths=box_width,
    patch_artist=True,
    boxprops=dict(facecolor="skyblue"),
    medianprops=dict(color="black", linewidth=1.5),
    flierprops=dict(marker="o", markerfacecolor="none", markeredgecolor=outlier_colors["St. Paul"]),
#     showmeans=True,
    meanprops=dict(marker="D", markeredgecolor="black", markerfacecolor="black")
)

# 绘制Greenville箱线图
bp2 = axes[0].boxplot(
    time_saving_greenville,
    positions=np.arange(11, 16) + position_offset,
    widths=box_width,
    patch_artist=True,
    boxprops=dict(facecolor="lightgreen"),
    medianprops=dict(color="black", linewidth=1.5),
    flierprops=dict(marker="o", markerfacecolor="none", markeredgecolor=outlier_colors["Greenville"]),
#     showmeans=True,
    meanprops=dict(marker="D", markeredgecolor="black", markerfacecolor="black")
)

# 设置子图0参数
axes[0].set_title("(a) Percentage Time Saving", y=-0.2, fontsize=16)  # 标题在正下方
axes[0].set_xlabel("Number of Units (N)", fontsize=14)
axes[0].set_ylabel("Time Saving (%)", fontsize=14)
axes[0].set_xticks(range(11, 16))
axes[0].set_xticklabels(range(11, 16))
axes[0].set_ylim(ylim_config["Time Saving"])
for ax in axes:
    ax.tick_params(axis='both', which='major', labelsize=14)
axes[0].legend([bp1["boxes"][0], bp2["boxes"][0]], ["St. Paul", "Greenville"], fontsize=14)

# ------------------ 绘制MPRE对比图 ------------------
# 准备数据
mpre_st_paul = [st_paul_data[(st_paul_data["Metric"] == "MPRE") & (st_paul_data["N"] == n)]["Value"].values for n in range(11, 16)]
mpre_greenville = [greenville_data[(greenville_data["Metric"] == "MPRE") & (greenville_data["N"] == n)]["Value"].values for n in range(11, 16)]

# 绘制St. Paul箱线图
bp3 = axes[1].boxplot(
    mpre_st_paul,
    positions=np.arange(11, 16) - position_offset,
    widths=box_width,
    patch_artist=True,
    boxprops=dict(facecolor="skyblue"),
    medianprops=dict(color="black", linewidth=1.5),
    flierprops=dict(marker="o", markerfacecolor="none", markeredgecolor=outlier_colors["St. Paul"]),
#     showmeans=True,
    meanprops=dict(marker="D", markeredgecolor="black", markerfacecolor="black")
)

# 绘制Greenville箱线图
bp4 = axes[1].boxplot(
    mpre_greenville,
    positions=np.arange(11, 16) + position_offset,
    widths=box_width,
    patch_artist=True,
    boxprops=dict(facecolor="lightgreen"),
    medianprops=dict(color="black", linewidth=1.5),
    flierprops=dict(marker="o", markerfacecolor="none", markeredgecolor=outlier_colors["Greenville"]),
#     showmeans=True,
    meanprops=dict(marker="D", markeredgecolor="black", markerfacecolor="black")
)

# 设置子图1参数
axes[1].set_title("(b) MPRE for Steady State Probability Distribution", y=-0.2, fontsize=16)  # 标题在正下方
axes[1].set_xlabel("Number of Units (N)", fontsize=14)
axes[1].set_ylabel("MPRE (%)", fontsize=14)
axes[1].set_xticks(range(11, 16))
axes[1].set_xticklabels(range(11, 16))
axes[1].set_ylim(ylim_config["MPRE"])
for ax in axes:
    ax.tick_params(axis='both', which='major', labelsize=14)
axes[1].legend([bp3["boxes"][0], bp4["boxes"][0]], ["St. Paul", "Greenville"], fontsize=14)

# ------------------ 保存图片 ------------------
plt.savefig("TimeSaving&MPRE-colored.pdf", dpi=600, bbox_inches="tight")
plt.show()
