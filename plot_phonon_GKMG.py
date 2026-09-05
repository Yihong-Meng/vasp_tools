# phonopy 4.4.0 band.yaml 绘图脚本 —— PRB/APL 期刊风格
# 用法: MPLBACKEND=Agg python3 plot_band.py  ->  band.png (300 dpi)

import yaml
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ---------- 可调参数 ----------
FIG = (3.5, 3.2)        # 单栏尺寸(英寸)
COLOR = 'k'             # 支线颜色: k=黑(期刊标准), 想要红改 'tab:red'
LW = 0.9                # 支线粗细
DPI = 300               # 分辨率
SHOW_CM1 = True         # True=右侧 cm-1 副轴, False=去掉
# -----------------------------

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['STIXGeneral', 'Times New Roman', 'DejaVu Serif']
plt.rcParams['mathtext.fontset'] = 'stix'

d = yaml.safe_load(open('band.yaml'))
ph = d['phonon']
seg = list(d['segment_nqpoint'])

dist = np.array([q['distance'] for q in ph])
freqs = np.array([[b['frequency'] for b in q['band']] for q in ph])
print('nqpoint:', dist.size, 'branches:', freqs.shape[1], 'segments:', seg)

boundary_idx = np.concatenate([[0], np.cumsum(seg) - 1])
tick_pos = dist[boundary_idx]
tick_labels = [chr(915), 'K', 'M', chr(915)]   # G K M G

fig, ax = plt.subplots(figsize=FIG)
for i in range(freqs.shape[1]):
    ax.plot(dist, freqs[:, i], color=COLOR, lw=LW)

ax.axhline(0.0, color='0.6', lw=0.6)                    # 零线(淡灰)
for x in tick_pos[1:-1]:
    ax.axvline(x, color='0.6', lw=0.6, ls='--')          # 高对称点虚线

ax.tick_params(direction='in', top=False, bottom=False, left=True, length=4, width=1, labelsize=11)
ax.set_xticks(tick_pos)
ax.set_xticklabels(tick_labels, fontsize=14)
ax.set_xlim(dist.min(), dist.max())
ax.set_xlabel('Wave vector', fontsize=12)
ax.set_ylabel('Frequency (THz)', fontsize=12)

if SHOW_CM1:
    ax2 = ax.twinx()
    ax2.set_ylim(np.array(ax.get_ylim()) * 33.356)
    ax2.tick_params(direction='in', top=False, labelsize=10)
    ax2.set_ylabel('cm$^{-1}$', fontsize=11)

plt.tight_layout()
plt.savefig('band.png', dpi=DPI)
print('saved band.png  ticks:', np.round(tick_pos, 4))
