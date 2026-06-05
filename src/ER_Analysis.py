import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from data_parser import load_parsed
from ref_poly import ref_poly

zip_path = "../dat/HY202103"
base_res_dir = "../res/png"
target_wafers = ['D07', 'D08', 'D23', 'D24']

print("🚀 평탄화 로직 기반 ER 추출 및 날짜별/Center vs Edge 분석을 시작합니다...")

# 1. WaferMap과 BoxPlot 최상위 경로 설정
wafer_map_dir = os.path.join(base_res_dir, "WaferMap")
box_plot_dir = os.path.join(base_res_dir, "BoxPlot")
os.makedirs(wafer_map_dir, exist_ok=True)
os.makedirs(box_plot_dir, exist_ok=True)

er_data_list = []
count = 0

for d in load_parsed(zip_path, target_wafers):
    date_str = d.get('date', 'Unknown_Date')

    m = (d['ref_data']['wl'] >= d['wl_min']) & (d['ref_data']['wl'] <= d['wl_max'])
    v_ref_wl, v_ref_il = d['ref_data']['wl'][m], d['ref_data']['il'][m]
    if len(v_ref_wl) < 4: continue

    poly_func = ref_poly(v_ref_wl, v_ref_il, smooth=False)
    max_er = 0.0

    for b in d['bias_data_list']:
        m_b = (b['wl'] >= d['wl_min']) & (b['wl'] <= d['wl_max'])
        v_wl, v_il = b['wl'][m_b], b['il'][m_b]
        if len(v_wl) == 0: continue

        flat_il = v_il - poly_func(v_wl)
        peaks, _ = find_peaks(flat_il, prominence=3.0, distance=30)

        if len(peaks) >= 2:
            trend = np.poly1d(np.polyfit(v_wl[peaks], flat_il[peaks], 1))
            final_il = flat_il - trend(v_wl)
        else:
            final_il = flat_il

        er = np.percentile(final_il, 99) - np.percentile(final_il, 1)
        max_er = max(max_er, er)

    er_data_list.append({
        'Wafer': d['wafer_id'], 'Band': d['band'], 'Date': date_str, 'Column': d['die_c'], 'Row': d['die_r'],
        'Radius': np.sqrt(d['die_c'] ** 2 + d['die_r'] ** 2), 'ER': max_er
    })

    count += 1
    if count % 100 == 0: print(f"현재 {count}개 Die 파싱 및 평탄화 ER 추출 완료...")

df = pd.DataFrame(er_data_list)

# ==========================================================
# 2. 자정을 넘긴 측정 데이터(0603, 0604) 병합 처리
# ==========================================================
def merge_midnight_dates(date_val):
    date_str = str(date_val)
    if '0603' in date_str or '0604' in date_str:
        return '20190603'
    return date_str

df['Date'] = df['Date'].apply(merge_midnight_dates)
print("🕒 0603 및 0604 날짜 데이터를 '0603_0604_Combined'로 통합했습니다.")
# ==========================================================

filtered_df = pd.DataFrame()

for (wafer, band, date), group in df.groupby(['Wafer', 'Band', 'Date']):
    if len(group) > 5:
        m_val, s_val = group['ER'].mean(), group['ER'].std()
        filtered_df = pd.concat(
            [filtered_df, group[(group['ER'] >= m_val - 3 * s_val) & (group['ER'] <= m_val + 3 * s_val)]])
    else:
        filtered_df = pd.concat([filtered_df, group])

max_rad = filtered_df['Radius'].max()
edge_thresh = max_rad * 0.75
filtered_df['Region'] = np.where(filtered_df['Radius'] > edge_thresh, 'Edge', 'Center')

# ==========================================================
# [통일된 디자인] Wafer Map 그리기 (격자 및 좌표 추가 수정본)
# ==========================================================
print("▶ 웨이퍼 및 날짜별 Wafer Map 생성 중...")

# 격자 범위를 결정하기 위해 물리적 반경(`max_rad`)을 기준으로 축 Limits 설정 (대칭)
map_limit = np.ceil(max_rad) + 0.5

for b in filtered_df['Band'].unique():
    band_limits = {'min': np.floor(filtered_df[filtered_df['Band'] == b]['ER'].min()),
                   'max': np.ceil(filtered_df[filtered_df['Band'] == b]['ER'].max())}

    for (w, date), group in filtered_df[filtered_df['Band'] == b].groupby(['Wafer', 'Date']):
        if group.empty: continue

        # figsize를 정방형으로 유지
        plt.figure(figsize=(10, 10))
        ax = plt.gca()  # 현재 Axes 객체 가져오기

        # 1. 웨이퍼 외곽선 및 Edge 영역 구분선 (zorder=1: 맨 아래)
        theta = np.linspace(0, 2 * np.pi, 100)
        ax.plot((max_rad + 0.5) * np.cos(theta), (max_rad + 0.5) * np.sin(theta), color='#555555', lw=2, zorder=1)
        ax.plot(edge_thresh * np.cos(theta), edge_thresh * np.sin(theta), color='#FF8888', ls='--', lw=2, alpha=0.7,
                zorder=1)

        # 2. 바둑판 격자 및 축 설정
        ax.set_aspect('equal')  # 정방형 비율 유지

        # 축 범위 설정 (0,0 기준 대칭)
        ax.set_xlim(-map_limit, map_limit)
        ax.set_ylim(-map_limit, map_limit)

        # Ticks 설정: 정수 좌표에 눈금을 매김
        ticks = np.arange(-np.floor(map_limit), np.ceil(map_limit) + 1, 1)
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)

        # 격자 추가 (zorder=2: 외곽선 위, 데이터 아래)
        ax.grid(True, which='major', axis='both', color='#DDDDDD', linestyle='--', linewidth=1, zorder=2)

        # 축 라벨 스타일 및 표시 설정
        ax.tick_params(axis='both', which='major', labelsize=11, labelcolor='#333333')
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontweight('bold')

        # 축 이름 설정
        ax.set_xlabel("Column (X Coordinate)", fontsize=14, fontweight='bold', labelpad=10)
        ax.set_ylabel("Row (Y Coordinate)", fontsize=14, fontweight='bold', labelpad=10)

        # 3. 데이터 포인트 플로팅 (zorder=5: 격자 위)
        # ER 특성에 맞게 cmap='coolwarm_r' 유지
        sc = ax.scatter(group['Column'], group['Row'], c=group['ER'], cmap='coolwarm_r',
                        vmin=band_limits['min'], vmax=band_limits['max'], s=500, edgecolor='black', alpha=0.9,
                        zorder=5)

        # 4. 데이터 값 텍스트 표시 (zorder=6: 데이터 포인트 위)
        for _, row in group.iterrows():
            ax.text(row['Column'], row['Row'], f"{row['ER']:.1f}", ha='center', va='center', fontsize=9,
                    color='black', fontweight='bold', zorder=6)

        # 컬러바 설정
        cb = plt.colorbar(sc, shrink=0.8, pad=0.03)
        cb.set_label('Extinction Ratio [dB]', fontsize=13, fontweight='bold')
        cb.ax.tick_params(labelsize=11)
        for l in cb.ax.yaxis.get_ticklabels(): l.set_weight("bold")

        # 제목 설정
        plt.title(f"Wafer Map: {w} / {b} (Flattened ER)\nDate: {date}", fontsize=17, fontweight='bold', pad=20)

        # 경로: WaferMap / 웨이퍼 / 날짜
        w_dir = os.path.join(wafer_map_dir, w, date)
        os.makedirs(w_dir, exist_ok=True)
        # bbox_inches='tight'로 축 라벨이 그림 밖으로 잘 나가지 않게 보호
        plt.savefig(os.path.join(w_dir, f"Map_{w}_{b}_{date}_ER.png"), bbox_inches='tight', dpi=100)
        plt.close()
# ==========================================================
# 3. [통일된 디자인] Box Plot 그리기 (개별 웨이퍼 단위)
# ==========================================================
print("▶ 웨이퍼 및 날짜별 Box Plot 생성 중...")
for (w, b, date), wbd_df in filtered_df.groupby(['Wafer', 'Band', 'Date']):
    if wbd_df.empty: continue

    plt.figure(figsize=(8, 8))  # 개별 웨이퍼 플롯 사이즈
    pos, data, labels, colors = [], [], [], []

    c = wbd_df[wbd_df['Region'] == 'Center']['ER'].values
    e = wbd_df[wbd_df['Region'] == 'Edge']['ER'].values

    if len(c) > 0:
        pos.append(1)
        data.append(c)
        labels.append(f"Center\nn={len(c)}")
        colors.append('#3498db')
    if len(e) > 0:
        pos.append(2)
        data.append(e)
        labels.append(f"Edge\nn={len(e)}")
        colors.append('#e74c3c')

    if not data: continue

    # Y축 범위를 미리 계산하여 배경색 영역 지정에 사용
    all_data = np.concatenate(data)
    y_min = min(all_data.min(), 15.0) - 2  # 데이터 최소값 또는 15dB 중 작은 값 기준 안전마진
    y_max = max(all_data.max(), 25.0) + 2  # 데이터 최대값 또는 25dB 중 큰 값 기준 안전마진
    plt.ylim(y_min, y_max)

    # ------------------------------------------------------
    # 🌟 성능별 배경색(수평 영역) 지정
    # ------------------------------------------------------
    # 타겟(20.0) 이상: 연한 초록색 (Good)
    plt.axhspan(20.0, y_max, facecolor='#e8f8f5', alpha=0.6, zorder=0, label='Good Region')

    # 타겟(20.0) 미만: 연한 붉은/코랄 계열 (Poor) - 눈에 띄면서 부정적인 의미 전달
    plt.axhspan(y_min, 20.0, facecolor='#fadbd8', alpha=0.6, zorder=0, label='Poor Region')
    # ------------------------------------------------------

    box = plt.boxplot(data, positions=pos, patch_artist=True, widths=0.5,
                      flierprops=dict(marker='d', markerfacecolor='black', markersize=6, alpha=0.6),
                      zorder=2)
    for p, c_hex in zip(box['boxes'], colors): p.set_facecolor(c_hex); p.set_alpha(0.7)

    # Jitter
    for p, d_arr in zip(pos, data):
        plt.scatter(np.random.normal(p, 0.05, len(d_arr)), d_arr, color='black', alpha=0.5, s=20, zorder=3)

    # 평균 및 타겟 라인
    avg_er = wbd_df['ER'].mean()
    plt.axhline(avg_er, color='blue', ls='--', lw=2.5, label=f'Avg: {avg_er:.2f} dB', zorder=4)
    plt.axhline(20.0, color='red', ls='-', lw=2.5, label='Target: 20.00 dB', zorder=4)

    plt.title(f"ER Analysis : {w} ({b})\nDate: {date}", fontsize=18, fontweight='bold', pad=15)
    plt.xticks(pos, labels, fontsize=14, fontweight='bold')
    plt.yticks(fontsize=14, fontweight='bold')
    plt.ylabel('Extinction Ratio [dB]', fontsize=16, fontweight='bold')

    plt.legend(loc='upper right', prop={'size': 11, 'weight': 'bold'})
    plt.grid(True, axis='y', ls=':', alpha=0.4, zorder=1)
    plt.xlim(0.5, max(pos) + 0.5)
    plt.tight_layout()

    # 경로: BoxPlot / 웨이퍼 / 날짜
    box_dir = os.path.join(box_plot_dir, w, date)
    os.makedirs(box_dir, exist_ok=True)
    plt.savefig(os.path.join(box_dir, f"Box_{w}_{b}_{date}_ER_Flattened.png"), bbox_inches='tight')
    plt.close()