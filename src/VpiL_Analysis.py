import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, savgol_filter
from concurrent.futures import ProcessPoolExecutor
from data_parser import load_parsed
from ref_poly import q_sub, ref_poly


def _extract_vpil_data(args):
    d, L_length = args
    wafer, band, c, r = d['wafer_id'], d['band'], d['die_c'], d['die_r']
    radius = np.sqrt(c ** 2 + r ** 2)
    date_str = d.get('date', 'Unknown_Date')

    m = (d['ref_data']['wl'] >= d['wl_min']) & (d['ref_data']['wl'] <= d['wl_max'])
    v_ref_wl, v_ref_il = d['ref_data']['wl'][m], d['ref_data']['il'][m]

    if len(v_ref_wl) < 31:
        return None

    poly_func = ref_poly(v_ref_wl, v_ref_il, smooth=True)

    # 💡 [수정됨] 부동소수점 오차 방지를 위해 == 0.0 대신 절대값 비교 사용
    z_data = next((b for b in d['bias_data_list'] if b['bias'] is not None and abs(b['bias']) < 1e-3), None)
    if not z_data:
        return None

    m0 = (z_data['wl'] >= d['wl_min']) & (z_data['wl'] <= d['wl_max'])
    w0, i0 = z_data['wl'][m0], z_data['il'][m0]
    flat0 = savgol_filter(i0, 31, 3) - poly_func(w0)
    v0, _ = find_peaks(-flat0, prominence=0.3, distance=20)

    if len(v0) < 2:
        return None

    fsr = np.mean(np.diff(w0[v0]))
    cwl = w0[v0[np.argmin(np.abs(w0[v0] - d['target_wl']))]]

    volts, p_pi = [], []
    sh = fsr / 2.5
    for b in d['bias_data_list']:
        if b['bias'] is None: continue
        w, i = b['wl'], b['il']
        mb = (w >= d['wl_min']) & (w <= d['wl_max'])
        wb, ib = w[mb], i[mb]
        mloc = (wb >= cwl - sh) & (wb <= cwl + sh)
        if np.sum(mloc) < 5: continue
        flat = savgol_filter(ib[mloc], 11, 3) - poly_func(wb[mloc])
        vwl = q_sub(wb[mloc], flat)
        volts.append(b['bias'])
        p_pi.append(2.0 * (vwl - cwl) / fsr)

    if len(volts) < 5:
        return None

    volts, p_pi = np.array(volts), np.array(p_pi)

    fit = np.poly1d(np.polyfit(volts, p_pi, 2))
    slope_at_0 = max(np.abs(np.polyder(fit)(0.0)), 1e-5)
    vpil_0V = L_length / slope_at_0

    # 💡 데이터가 안 나온다면 여기가 원인일 수 있습니다! (범위 확인 필요)
    if not (0.1 <= vpil_0V <= 10.0):
        return None

    return {
        'Wafer': wafer, 'Band': band, 'Date': date_str, 'Column': c, 'Row': r,
        'Radius': radius, 'VpiL_0V': vpil_0V
    }


def main():
    zip_path = "../dat/HY202103"
    base_res_dir = "../res/png"
    target_wafers = ['D07', 'D08', 'D23', 'D24']
    L_length = 0.05
    TARGET_VPIL = {'LMZO': 1.4, 'LMZC': 2.0}

    wafer_map_dir = os.path.join(base_res_dir, "WaferMap")
    box_plot_dir = os.path.join(base_res_dir, "BoxPlot")
    os.makedirs(wafer_map_dir, exist_ok=True)
    os.makedirs(box_plot_dir, exist_ok=True)

    print("🚀 분석용 데이터 추출을 시작합니다...")

    parsed_data_list = load_parsed(zip_path, target_wafers)
    tasks = [(d, L_length) for d in parsed_data_list]
    summary_rows = []

    with ProcessPoolExecutor(max_workers=None) as ex:
        futures = [ex.submit(_extract_vpil_data, t) for t in tasks]
        for f in futures:
            res = f.result()
            if res is not None:
                summary_rows.append(res)

    if not summary_rows:
        print("❌ 유효한 데이터를 찾지 못했습니다. _extract_vpil_data의 필터 조건을 확인하세요.")
        return

    df = pd.DataFrame(summary_rows)

    def merge_midnight_dates(date_val):
        date_str = str(date_val)
        if '0603' in date_str or '0604' in date_str:
            return '20190603'
        return date_str

    df['Date'] = df['Date'].apply(merge_midnight_dates)
    print("🕒 날짜 병합 처리 완료.")

    filtered_df = pd.DataFrame()
    df_grouped = df.groupby(['Wafer', 'Band', 'Date', 'Column', 'Row', 'Radius'], as_index=False)['VpiL_0V'].mean()

    for (wafer, band, date), group in df_grouped.groupby(['Wafer', 'Band', 'Date']):
        if len(group) > 5:
            m, s = group['VpiL_0V'].mean(), group['VpiL_0V'].std()
            valid_group = group[(group['VpiL_0V'] >= m - 3 * s) & (group['VpiL_0V'] <= m + 3 * s)]
            filtered_df = pd.concat([filtered_df, valid_group])
        else:
            filtered_df = pd.concat([filtered_df, group])

    max_r = filtered_df['Radius'].max()
    edge_limit = max_r * 0.75
    filtered_df['Region'] = np.where(filtered_df['Radius'] > edge_limit, 'Edge', 'Center')

    print(f"✅ 데이터 추출 및 필터링 완료! (총 {len(filtered_df)}개 다이 분석)")

    # ==========================================================
    # Wafer Map 그리기 (격자 및 좌표 추가 수정본)
    # ==========================================================
    print("▶ 웨이퍼 및 날짜별 Wafer Map 생성 중...")
    band_limits = {}
    for b in filtered_df['Band'].unique():
        b_data = filtered_df[filtered_df['Band'] == b]['VpiL_0V']
        band_limits[b] = {'min': b_data.min() - 0.05, 'max': b_data.max() + 0.05}

    # 격자 범위를 결정하기 위해 물리적 반경(`max_r`)을 기준으로 축 Limits 설정 (대칭)
    # 0.5 여유를 두어 원 가장자리 데이터가 격자에 걸치지 않게 함
    map_limit = np.ceil(max_r) + 0.5

    for (wafer, band, date), group in filtered_df.groupby(['Wafer', 'Band', 'Date']):
        # figsize를 정방형으로 유지
        plt.figure(figsize=(10, 10))
        ax = plt.gca()  # 현재 Axes 객체 가져오기

        # 1. 웨이퍼 외곽선 및 Edge 영역 구분선 (zorder=1: 맨 아래)
        theta = np.linspace(0, 2 * np.pi, 100)
        ax.plot((max_r + 0.5) * np.cos(theta), (max_r + 0.5) * np.sin(theta), color='#555555', lw=2, zorder=1)  # 외곽선
        ax.plot(edge_limit * np.cos(theta), edge_limit * np.sin(theta), color='#FF8888', ls='--', lw=2, alpha=0.7,
                zorder=1)  # Edge 구분

        # 2. 바둑판 격자 및 축 설정
        ax.set_aspect('equal')  # 정방형 비율 유지

        # 축 범위 설정 (0,0 기준 대칭)
        ax.set_xlim(-map_limit, map_limit)
        ax.set_ylim(-map_limit, map_limit)

        # Ticks 설정: 정수 좌표에 눈금을 매김 (zorder=0으로 설정하기 위해 grid 사용 권장)
        ticks = np.arange(-np.floor(map_limit), np.ceil(map_limit) + 1, 1)
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)

        # 격자 추가: 주요 눈금(Major Ticks) 위에 회색 점선 격자 (zorder=2: 외곽선 위, 데이터 아래)
        ax.grid(True, which='major', axis='both', color='#DDDDDD', linestyle='--', linewidth=1, zorder=2)

        # 축 라벨 스타일 및 표시 설정 (plt.axis('off')를 제거했으므로 보임)
        ax.tick_params(axis='both', which='major', labelsize=11, labelcolor='#333333')
        # Ticks 굵게 (Optional)
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontweight('bold')

        # 축 이름 설정
        ax.set_xlabel("Column (X Coordinate)", fontsize=14, fontweight='bold', labelpad=10)
        ax.set_ylabel("Row (Y Coordinate)", fontsize=14, fontweight='bold', labelpad=10)

        # 컬러 바 범위
        v_min, v_max = band_limits[band]['min'], band_limits[band]['max']

        # 3. 데이터 포인트 플로팅 (zorder=5: 격자 위)
        # 데이터가 많아 격자가 좁아 보일 경우 s(size)를 약간 줄여도 좋습니다 (예: 400).
        scatter = ax.scatter(group['Column'], group['Row'], c=group['VpiL_0V'], cmap='coolwarm',
                             vmin=v_min, vmax=v_max, s=500, edgecolor='black', alpha=0.9, zorder=5)

        # 4. 데이터 값 텍스트 표시 (zorder=6: 데이터 포인트 위)
        for _, row in group.iterrows():
            ax.text(row['Column'], row['Row'], f"{row['VpiL_0V']:.2f}",
                    ha='center', va='center', fontsize=9, weight='bold', color='black', zorder=6)

        # 컬러바 설정
        cb = plt.colorbar(scatter, shrink=0.8, pad=0.03)
        cb.set_label('Vpi*L @ 0V [V*cm]', fontsize=13, fontweight='bold')
        cb.ax.tick_params(labelsize=11)
        for label in cb.ax.yaxis.get_ticklabels(): label.set_weight("bold")

        # 제목 설정
        plt.title(f"Wafer Map: {wafer} / {band} (VpiL)\nDate: {date}", fontsize=17, fontweight='bold', pad=20)

        # [수정] plt.axis('off') 제거 (축을 보이게 함)
        # plt.tight_layout() # 격자와 라벨이 잘리지 않게 조정

        # 저장
        w_dir = os.path.join(wafer_map_dir, wafer, date)
        os.makedirs(w_dir, exist_ok=True)
        # bbox_inches='tight'는 축 라벨이 그림 밖으로 잘 나가지 않게 해줍니다.
        plt.savefig(os.path.join(w_dir, f"Map_{wafer}_{band}_{date}_VpiL_0V.png"), bbox_inches='tight', dpi=100)
        plt.close()
    # ==========================================================
    # Box Plot 그리기 (VpiL 맞춤형 디자인)
    # ==========================================================
    print("▶ 웨이퍼 및 날짜별 Box Plot 생성 중...")

    for (wafer, band, date), wbd_df in filtered_df.groupby(['Wafer', 'Band', 'Date']):
        current_target = TARGET_VPIL.get(band, 1.5)
        avg_vpil = wbd_df['VpiL_0V'].mean()

        plt.figure(figsize=(8, 8))
        pos, data, labels, colors = [], [], [], []

        c_data = wbd_df[wbd_df['Region'] == 'Center']['VpiL_0V'].values
        e_data = wbd_df[wbd_df['Region'] == 'Edge']['VpiL_0V'].values

        if len(c_data) > 0:
            pos.append(1)
            data.append(c_data)
            labels.append(f"Center\nn={len(c_data)}")
            colors.append('#3498db')
        if len(e_data) > 0:
            pos.append(2)
            data.append(e_data)
            labels.append(f"Edge\nn={len(e_data)}")
            colors.append('#e74c3c')

        if not data: continue

        all_data = np.concatenate(data)
        y_min = min(all_data.min(), current_target) - 0.2
        y_max = max(all_data.max(), current_target) + 0.2
        plt.ylim(y_min, y_max)

        # [VpiL 맞춤] 작을수록 좋으므로 Target 아래가 초록색(Good)
        plt.axhspan(y_min, current_target, facecolor='#e8f8f5', alpha=0.6, zorder=0, label='Good Region')
        plt.axhspan(current_target, y_max, facecolor='#fadbd8', alpha=0.6, zorder=0, label='Poor Region')

        box = plt.boxplot(data, positions=pos, patch_artist=True, widths=0.5,
                          flierprops=dict(marker='d', markerfacecolor='black', markersize=6, alpha=0.6),
                          zorder=2)
        for p, c in zip(box['boxes'], colors): p.set_facecolor(c); p.set_alpha(0.7)

        for p, d_arr in zip(pos, data):
            plt.scatter(np.random.normal(p, 0.05, len(d_arr)), d_arr, color='black', alpha=0.5, s=20, zorder=3)

        plt.axhline(avg_vpil, color='blue', ls='--', lw=2.5, label=f'Avg: {avg_vpil:.2f} V*cm', zorder=4)
        plt.axhline(current_target, color='red', ls='-', lw=2.5, label=f'Target: {current_target:.2f} V*cm', zorder=4)

        plt.title(f"VpiL Analysis : {wafer} ({band})\nDate: {date}", fontsize=18, fontweight='bold', pad=15)
        plt.xticks(pos, labels, fontsize=14, fontweight='bold')
        plt.yticks(fontsize=14, fontweight='bold')
        plt.ylabel("Vpi*L @ 0V [V*cm]", fontsize=16, fontweight='bold')

        plt.legend(loc='upper right', prop={'size': 11, 'weight': 'bold'})
        plt.grid(True, axis='y', ls=':', alpha=0.4, zorder=1)
        plt.xlim(0.5, max(pos) + 0.5)
        plt.tight_layout()

        box_dir = os.path.join(box_plot_dir, wafer, date)
        os.makedirs(box_dir, exist_ok=True)
        plt.savefig(os.path.join(box_dir, f"Box_{wafer}_{band}_{date}_VpiL_0V.png"), bbox_inches='tight')
        plt.close()

    print("✅ 모든 그래프가 웨이퍼별/날짜별 폴더 구조로 저장되었습니다!")


if __name__ == "__main__":
    main()