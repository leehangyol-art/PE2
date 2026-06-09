# MZM 측정 데이터 분석용 공용 모듈.
# 무거운 계산(파싱/피팅/지표 추출)을 여기에 모아두고,
# 주피터 노트북은 이 함수들만 import 해서 화면에 결과를 바로 띄운다.
# 경로와 웨이퍼 목록을 인자로 받기 때문에 데이터가 달라져도 코드 수정이 필요 없다.

import os
import re
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, savgol_filter
from data_parser import parse_wafer_data   # 기존 XML 파서 재사용
from ref_poly import q_sub, ref_poly       # 기존 공용 헬퍼 재사용

VPIL_MIN, VPIL_MAX = 0.1, 10.0   # 유효 VpiL 범위 (V*cm)
DEFAULT_L_UM = 500               # 길이 못 읽었을 때 기본값 (um)


def _show_fig(fig):
    # 그림을 주피터 화면에 띄우고 닫는다. matplotlib 백엔드 설정에 의존하지 않아
    # VS Code / 주피터 어디서든 안정적으로 표시된다.
    try:
        from IPython.display import display
        display(fig)
    except Exception:
        pass
    plt.close(fig)


# ============================================================
# 1. 데이터 로딩 / 목록 (드롭다운용)
# ============================================================
def load(folder, wafers):
    # 폴더 경로 + 웨이퍼 목록 -> 다이별 데이터 리스트
    return list(parse_wafer_data(folder, wafers))


def build_length_map(folder, wafers):
    # 다이별 arm 길이(um)를 XML 모듈레이터 이름에서 읽어 표로 만든다.
    # 이름 예: 'MZMCTE_LULAB_450_500' -> 뒤 500 = 길이(um)
    out = {}
    for root_dir, _, files in os.walk(folder):
        for fname in files:
            if not fname.lower().endswith('.xml'):
                continue
            fp = os.path.join(root_dir, fname)
            if not any(w in fp for w in wafers):
                continue
            try:
                root = ET.parse(fp).getroot()
            except Exception:
                continue
            info = root.find('.//TestSiteInfo')
            if info is None:
                continue
            ts = info.get('TestSite', '').upper()
            band = 'LMZO' if 'LMZO' in ts else ('LMZC' if 'LMZC' in ts else None)
            if band is None:
                continue
            die_c = int(info.get('DieRow', 0)) if info.get('DieRow') else 0
            die_r = int(info.get('DieColumn', 0)) if info.get('DieColumn') else 0
            wafer = next((w for w in wafers if w in fp), None)
            if wafer is None:
                continue
            for mod in root.iter('Modulator'):
                nm = mod.attrib.get('Name', '')
                if 'ALIGN' in nm.upper():
                    continue
                m = re.search(r'_LULAB_\d+_(\d+)', nm)
                if m:
                    out[(wafer, die_c, die_r, band)] = int(m.group(1))
                break
    return out


def list_wafers(data):
    return sorted({d['wafer_id'] for d in data})


def list_dates(data, wafer):
    return sorted({d['date'] for d in data if d['wafer_id'] == wafer})


def list_dies(data, wafer, date):
    # 그 웨이퍼/날짜의 (밴드, 열, 행) 목록
    items = {(d['band'], d['die_c'], d['die_r'])
             for d in data if d['wafer_id'] == wafer and d['date'] == date}
    return sorted(items)


def get_die(data, wafer, date, band, c, r):
    # 선택한 다이 한 개의 데이터 찾기
    for d in data:
        if (d['wafer_id'], d['date'], d['band'], d['die_c'], d['die_r']) == (wafer, date, band, c, r):
            return d
    return None


# ============================================================
# 2. 지표 추출 (계산)
# ============================================================
def extract_vpil(d, L_cm):
    # 다이 1개의 0V 기준 V_pi*L 계산. 반환: (vpil_0V, volts, vpil_curve) 또는 None
    m = (d['ref_data']['wl'] >= d['wl_min']) & (d['ref_data']['wl'] <= d['wl_max'])
    ref_wl, ref_il = d['ref_data']['wl'][m], d['ref_data']['il'][m]
    if len(ref_wl) < 31:
        return None
    poly = ref_poly(ref_wl, ref_il, smooth=True)            # 평탄화 기준선

    z = next((b for b in d['bias_data_list']
              if b['bias'] is not None and abs(b['bias']) < 1e-3), None)  # 0V 측정
    if not z:
        return None
    m0 = (z['wl'] >= d['wl_min']) & (z['wl'] <= d['wl_max'])
    w0, i0 = z['wl'][m0], z['il'][m0]
    flat0 = savgol_filter(i0, 31, 3) - poly(w0)             # 0V 평탄화
    v0, _ = find_peaks(-flat0, prominence=0.3, distance=20) # 깊은 골(null)
    if len(v0) < 2:
        return None
    fsr = np.mean(np.diff(w0[v0]))                          # 골 간격 평균 = FSR
    cwl = w0[v0[np.argmin(np.abs(w0[v0] - d['target_wl']))]]  # 기준 null

    volts, p_pi = [], []
    sh = fsr / 2.5
    for b in d['bias_data_list']:
        if b['bias'] is None:
            continue
        mb = (b['wl'] >= d['wl_min']) & (b['wl'] <= d['wl_max'])
        wb, ib = b['wl'][mb], b['il'][mb]
        mloc = (wb >= cwl - sh) & (wb <= cwl + sh)          # null 주변만
        if np.sum(mloc) < 5:
            continue
        flat = savgol_filter(ib[mloc], 11, 3) - poly(wb[mloc])
        volts.append(b['bias'])
        p_pi.append(2.0 * (q_sub(wb[mloc], flat) - cwl) / fsr)  # 위상(pi 단위)

    if len(volts) < 5:
        return None
    volts, p_pi = np.array(volts), np.array(p_pi)
    deriv = np.polyder(np.poly1d(np.polyfit(volts, p_pi, 2)))   # 위상-전압 기울기
    vpil_0V = L_cm / max(abs(deriv(0.0)), 1e-5)             # V_pi*L = 길이 / 기울기
    if not (VPIL_MIN <= vpil_0V <= VPIL_MAX):
        return None
    vpil_curve = L_cm / np.maximum(np.abs(deriv(volts)), 1e-5)
    return vpil_0V, volts, vpil_curve


def extract_phase(d):
    # 전압별 null 이동을 라디안 위상으로. 반환: (volts, phase_rad, fsr) 또는 None
    m = (d['ref_data']['wl'] >= d['wl_min']) & (d['ref_data']['wl'] <= d['wl_max'])
    ref_wl, ref_il = d['ref_data']['wl'][m], d['ref_data']['il'][m]
    if len(ref_wl) < 31:
        return None
    poly = ref_poly(ref_wl, ref_il, smooth=True)
    z = next((b for b in d['bias_data_list']
              if b['bias'] is not None and abs(b['bias']) < 1e-3), None)
    if not z:
        return None
    m0 = (z['wl'] >= d['wl_min']) & (z['wl'] <= d['wl_max'])
    w0, i0 = z['wl'][m0], z['il'][m0]
    flat0 = savgol_filter(i0, 31, 3) - poly(w0)
    v0, _ = find_peaks(-flat0, prominence=0.3, distance=20)
    if len(v0) < 2:
        return None
    fsr = np.mean(np.diff(w0[v0]))
    cwl = w0[v0[np.argmin(np.abs(w0[v0] - d['target_wl']))]]
    sh = fsr / 2.5
    pts = []
    for b in d['bias_data_list']:
        if b['bias'] is None:
            continue
        mb = (b['wl'] >= d['wl_min']) & (b['wl'] <= d['wl_max'])
        wb, ib = b['wl'][mb], b['il'][mb]
        mloc = (wb >= cwl - sh) & (wb <= cwl + sh)
        if np.sum(mloc) < 5:
            continue
        flat = savgol_filter(ib[mloc], 11, 3) - poly(wb[mloc])
        vwl = q_sub(wb[mloc], flat)
        # 위상(라디안): 한 FSR 이동 = 2*pi
        pts.append((b['bias'], 2.0 * np.pi * (vwl - cwl) / fsr))
    if len(pts) < 3:
        return None
    pts.sort()
    V = np.array([p[0] for p in pts])
    PH = np.array([p[1] for p in pts])
    return V, PH, fsr


def extract_il(d):
    # 삽입손실(IL) = 모든 전압 중 통과대역 최대 투과율 (가장 손실 적은 지점)
    mx = -np.inf
    for b in d['bias_data_list']:
        if b['bias'] is None:
            continue
        m = (b['wl'] >= d['wl_min']) & (b['wl'] <= d['wl_max'])
        if np.any(m):
            mx = max(mx, float(np.max(b['il'][m])))
    return mx if mx > -np.inf else None


def extract_er(d):
    # 소광비(ER) = 평탄화 후 (99퍼센타일 - 1퍼센타일) 의 전압별 최댓값
    m = (d['ref_data']['wl'] >= d['wl_min']) & (d['ref_data']['wl'] <= d['wl_max'])
    ref_wl, ref_il = d['ref_data']['wl'][m], d['ref_data']['il'][m]
    if len(ref_wl) < 4:
        return None
    poly = ref_poly(ref_wl, ref_il, smooth=False)
    best = 0.0
    for b in d['bias_data_list']:
        mb = (b['wl'] >= d['wl_min']) & (b['wl'] <= d['wl_max'])
        wl, il = b['wl'][mb], b['il'][mb]
        if len(wl) == 0:
            continue
        flat = il - poly(wl)
        pk, _ = find_peaks(flat, prominence=3.0, distance=30)
        if len(pk) >= 2:                                    # 봉우리 기울기 보정
            flat = flat - np.poly1d(np.polyfit(wl[pk], flat[pk], 1))(wl)
        best = max(best, float(np.percentile(flat, 99) - np.percentile(flat, 1)))
    return best if best > 0 else None


def metric_value(d, metric, length_map=None):
    # 지표 이름으로 스칼라 값 하나 반환 (웨이퍼맵/박스플롯용)
    if metric == 'vpil':
        L_cm = (length_map or {}).get(
            (d['wafer_id'], d['die_c'], d['die_r'], d['band']), DEFAULT_L_UM) * 1e-4
        res = extract_vpil(d, L_cm)
        return res[0] if res else None
    if metric == 'il':
        return extract_il(d)
    if metric == 'er':
        return extract_er(d)
    return None


METRIC_LABEL = {'vpil': 'Vpi*L (V*cm)', 'il': 'IL (dB)', 'er': 'ER (dB)'}


# ============================================================
# 3. 시각화 (그려서 화면에 보여줌 - 노트북 인라인용)
# ============================================================
def plot_die_vpil(data, length_map, wafer, date, band, c, r):
    # 선택한 다이 1개의 VpiL 곡선
    d = get_die(data, wafer, date, band, c, r)
    fig, ax = plt.subplots(figsize=(8, 5))
    res = extract_vpil(d, (length_map or {}).get((wafer, c, r, band), DEFAULT_L_UM) * 1e-4) if d else None
    if not res:
        ax.text(0.5, 0.5, '유효한 VpiL 데이터 없음', ha='center', va='center')
        _show_fig(fig)
        return
    vpil, volts, curve = res
    ax.plot(volts, curve, 's-', color='gray', label='VpiL curve')
    ax.plot(0.0, vpil, 'r*', markersize=16, label=f'@0V = {vpil:.3f} V*cm')
    ax.axhline(vpil, color='red', ls=':', alpha=0.5)
    ax.set_xlabel('Voltage (V)'); ax.set_ylabel('Vpi*L (V*cm)')
    ax.set_title(f'{wafer} {band} ({c},{r})  |  VpiL@0V = {vpil:.3f} V*cm')
    ax.grid(alpha=0.3); ax.legend()
    fig.tight_layout()
    _show_fig(fig)


def plot_die_phase(data, wafer, date, band, c, r):
    # 선택한 다이의 위상변화(라디안) + pi 기준선
    d = get_die(data, wafer, date, band, c, r)
    fig, ax = plt.subplots(figsize=(8, 5))
    res = extract_phase(d) if d else None
    if not res:
        ax.text(0.5, 0.5, '유효한 위상 데이터 없음', ha='center', va='center')
        _show_fig(fig)
        return
    V, PH, fsr = res
    ax.plot(V, PH, 'o-', color='steelblue', label='Δφ (rad)')
    ax.axhline(np.pi, color='red', ls=':', alpha=0.6)      # pi = V_pi 도달점
    ax.text(V.max(), np.pi, '  π', color='red', va='center')
    ax.axhline(0, color='gray', lw=0.5)
    ax.set_xlabel('Voltage (V)'); ax.set_ylabel('Phase shift Δφ (rad)')
    ax.set_title(f'{wafer} {band} ({c},{r})  |  FSR = {fsr:.2f} nm')
    ax.grid(alpha=0.3); ax.legend()
    fig.tight_layout()
    _show_fig(fig)


def plot_die_plot(data, length_map, wafer, date, band, c, r):
    # [추가] 선택한 다이 1개의 Raw 스펙트럼 스펙트럼 플롯
    d = get_die(data, wafer, date, band, c, r)
    fig, ax = plt.subplots(figsize=(8, 5))
    if d is None:
        ax.text(0.5, 0.5, '유효한 데이터 없음', ha='center', va='center')
        _show_fig(fig)
        return

    # bias 데이터 플로팅
    for b in d['bias_data_list']:
        ax.plot(b['wl'], b['il'], label=b['label'], linewidth=2)

    # Reference 데이터 강조 플로팅
    ax.plot(d['ref_data']['wl'], d['ref_data']['il'], label=d['ref_data']['label'],
            linewidth=2, color='black', alpha=0.8, linestyle='--')

    ax.set_title(f"Wafer: {wafer} / Coord: ({c}, {r}) / Band: {band}", fontsize=12, fontweight='bold', pad=10)
    ax.set_xlabel("Wavelength (nm)", fontsize=11, fontweight='bold')
    ax.set_ylabel("Transmission (dB)", fontsize=11, fontweight='bold')
    ax.set_ylim(-65, 5)
    ax.legend(loc='best', prop={'size': 9, 'weight': 'bold'})
    ax.grid(True, linestyle='--', alpha=0.6, linewidth=1)

    for spine in ax.spines.values():
        spine.set_linewidth(1.5)

    fig.tight_layout()
    _show_fig(fig)


def plot_die_fitting(data, length_map, wafer, date, band, c, r):
    # [추가] 선택한 다이 1개의 리플 제거 및 피팅 스펙트럼 (Zoomed)
    d = get_die(data, wafer, date, band, c, r)
    fig, ax = plt.subplots(figsize=(8, 5))
    if d is None:
        ax.text(0.5, 0.5, '유효한 데이터 없음', ha='center', va='center')
        _show_fig(fig)
        return

    m = (d['ref_data']['wl'] >= d['wl_min']) & (d['ref_data']['wl'] <= d['wl_max'])
    v_ref_wl, v_ref_il = d['ref_data']['wl'][m], d['ref_data']['il'][m]

    if len(v_ref_wl) < 31:
        ax.text(0.5, 0.5, '데이터 부족 (N < 31)', ha='center', va='center')
        _show_fig(fig)
        return

    # 평탄화용 다항식 피팅 및 R^2 계산
    sm_ref = savgol_filter(v_ref_il, 31, 3)
    poly = np.poly1d(np.polyfit(v_ref_wl, sm_ref, 3))
    y_fitted = poly(v_ref_wl)

    ss_res = np.sum((sm_ref - y_fitted) ** 2)
    ss_tot = np.sum((sm_ref - np.mean(sm_ref)) ** 2)
    r_squared = 1.0 if ss_tot == 0 else 1 - (ss_res / ss_tot)

    # 줌 범위 지정을 위한 Peak 탐색 (-2V 기준 혹은 첫 번째 bias)
    z_min, z_max = d['wl_min'], d['wl_max']
    tgt_bias = next((b for b in d['bias_data_list'] if b['bias'] == -2.0),
                    d['bias_data_list'][0] if d['bias_data_list'] else None)

    if tgt_bias:
        m_t = (tgt_bias['wl'] >= d['wl_min']) & (tgt_bias['wl'] <= d['wl_max'])
        w_t, i_t = tgt_bias['wl'][m_t], tgt_bias['il'][m_t]

        if len(w_t) >= 31:
            flat_t = savgol_filter(i_t, 31, 3) - poly(w_t)
            peaks, _ = find_peaks(flat_t, prominence=3.0, distance=30)

            if len(peaks) >= 2:
                centers = (w_t[peaks[:-1]] + w_t[peaks[1:]]) / 2.0
                band_str = str(d.get('band', '')).upper()
                target_wl = d.get('target_wl', 1310.0) if 'O' in band_str else d.get('target_wl', 1550.0)
                idx = np.argmin(np.abs(centers - target_wl))

                if idx + 1 < len(peaks):
                    z_min, z_max = w_t[peaks[idx]] - 0.5, w_t[peaks[idx + 1]] + 0.5

    # 전압별 평탄화 및 봉우리 기울기 보정 후 플로팅
    for b in d['bias_data_list']:
        m_b = (b['wl'] >= d['wl_min']) & (b['wl'] <= d['wl_max'])
        v_wl, v_il = b['wl'][m_b], b['il'][m_b]
        if len(v_wl) < 31:
            continue

        flat_il = savgol_filter(v_il, 31, 3) - poly(v_wl)
        peaks, _ = find_peaks(flat_il, prominence=3.0, distance=30)

        if len(peaks) >= 2:
            flat_il -= np.poly1d(np.polyfit(v_wl[peaks], flat_il[peaks], 1))(v_wl)

        ax.plot(v_wl, flat_il, label=b['label'], alpha=0.8, linewidth=2)

    # Reference선 플로팅
    ax.plot(v_ref_wl, sm_ref - y_fitted, label=f'REF (R^2={r_squared:.4f})',
            color='black', lw=2, linestyle='--', zorder=10)

    ax.set_title(f"Wafer: {wafer} / Coord: ({c}, {r}) / Band: {band}\nSmoothed, Flattened & Zoomed", fontsize=12,
                 fontweight='bold', pad=10)
    ax.set_xlabel('Wavelength (nm)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Transmission (dB)', fontsize=11, fontweight='bold')
    ax.axhline(0, color='gray', ls='--', alpha=0.6, linewidth=1.5)
    ax.set_xlim(z_min, z_max)
    ax.legend(loc='best', prop={'size': 9, 'weight': 'bold'})
    ax.grid(True, linestyle='--', alpha=0.6, linewidth=1)

    for spine in ax.spines.values():
        spine.set_linewidth(1.5)

    fig.tight_layout()
    _show_fig(fig)

# metric별 스타일 (원본 분석 스크립트와 동일)
_MAP_CFG = {
    'vpil': dict(cmap='coolwarm',   cbar='Vpi*L @ 0V [V*cm]', maptag='VpiL',
                 boxtag='VpiL', unit='V*cm', target={'LMZO': 1.4, 'LMZC': 2.0}, good_above=False),
    'il':   dict(cmap='coolwarm_r', cbar='IL [dB]', maptag='IL',
                 boxtag='IL', unit='dB', target=-8.75, good_above=True),
    'er':   dict(cmap='coolwarm_r', cbar='Extinction Ratio [dB]', maptag='Flattened ER',
                 boxtag='ER', unit='dB', target=20.0, good_above=True),
}


def _merge_date(date):
    # 06/03 자정 넘겨 06/04 새벽까지 이어진 측정을 같은 날(0603)로 (원본과 동일)
    s = str(date)
    return '20190603' if ('0603' in s or '0604' in s) else s


def _metric_df(data, length_map, metric):
    # 원본 스크립트와 동일 전처리: 중복좌표 평균 -> 날짜병합 -> (그룹>5면) 3시그마 -> Center/Edge
    rows = []
    for d in data:
        v = metric_value(d, metric, length_map)
        if v is None:
            continue
        rows.append({'Wafer': d['wafer_id'], 'Band': d['band'], 'Date': _merge_date(d['date']),
                     'Column': d['die_c'], 'Row': d['die_r'],
                     'Radius': np.hypot(d['die_c'], d['die_r']), 'val': v})
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df = df.groupby(['Wafer', 'Band', 'Date', 'Column', 'Row', 'Radius'], as_index=False)['val'].mean()
    keep = []
    for _, g in df.groupby(['Wafer', 'Band', 'Date']):
        if len(g) > 5:
            mu, sd = g['val'].mean(), g['val'].std()
            g = g[(g['val'] >= mu - 3 * sd) & (g['val'] <= mu + 3 * sd)]
        keep.append(g)
    df = pd.concat(keep, ignore_index=True)
    df['Region'] = np.where(df['Radius'] > df['Radius'].max() * 0.75, 'Edge', 'Center')
    return df


def plot_wafer_map(data, length_map, wafer, date, band, metric='vpil', selected_c=None, selected_r=None):
    # res 폴더의 웨이퍼맵과 동일한 스타일로 그린다.
    cfg = _MAP_CFG[metric]
    df = _metric_df(data, length_map, metric)
    qdate = _merge_date(date)
    fig, ax = plt.subplots(figsize=(10, 10))
    grp = df[(df['Wafer'] == wafer) & (df['Band'] == band) & (df['Date'] == qdate)] if not df.empty else df
    if df.empty or grp.empty:
        ax.text(0.5, 0.5, '데이터 없음', ha='center', va='center')
        _show_fig(fig)
        return

    # 색 범위(밴드 전체 기준) + 지도 크기(전체 반경 기준)
    bvals = df[df['Band'] == band]['val']
    if metric == 'vpil':
        v_min, v_max = bvals.min() - 0.05, bvals.max() + 0.05
    else:
        v_min, v_max = np.floor(bvals.min()), np.ceil(bvals.max())
    max_r = df['Radius'].max()
    edge_limit = max_r * 0.75
    map_limit = np.ceil(max_r) + 0.5

    th = np.linspace(0, 2 * np.pi, 100)
    ax.plot((max_r + 0.5) * np.cos(th), (max_r + 0.5) * np.sin(th), color='#555555', lw=2, zorder=1)
    ax.plot(edge_limit * np.cos(th), edge_limit * np.sin(th), color='#FF8888', ls='--', lw=2, alpha=0.7, zorder=1)
    ax.set_aspect('equal')
    ax.set_xlim(-map_limit, map_limit);
    ax.set_ylim(-map_limit, map_limit)
    ticks = np.arange(-np.floor(map_limit), np.ceil(map_limit) + 1, 1)
    ax.set_xticks(ticks);
    ax.set_yticks(ticks)
    ax.grid(True, which='major', color='#DDDDDD', linestyle='--', linewidth=1, zorder=2)

    for lb in ax.get_xticklabels() + ax.get_yticklabels():
        lb.set_fontweight('bold')
    ax.set_xlabel('Column (X Coordinate)', fontsize=14, fontweight='bold', labelpad=10)
    ax.set_ylabel('Row (Y Coordinate)', fontsize=14, fontweight='bold', labelpad=10)

    # 전체 다이 마커 그리기
    sc = ax.scatter(grp['Column'], grp['Row'], c=grp['val'], cmap=cfg['cmap'],
                    vmin=v_min, vmax=v_max, s=500, edgecolor='black', alpha=0.9, zorder=5)

    for _, row in grp.iterrows():
        ax.text(row['Column'], row['Row'], f"{row['val']:.2f}", ha='center', va='center',
                fontsize=9, weight='bold', color='black', zorder=6)

    # 🌟 [여기가 핵심 추가 코드!] 선택한 다이가 있다면 그 위에 굵은 핫핑크색 링을 씌웁니다.
    if selected_c is not None and selected_r is not None:
        ax.scatter(selected_c, selected_r, s=750, facecolors='none',
                   edgecolors='#FF007F', linewidths=4, zorder=10, label='Selected Die')

    cb = plt.colorbar(sc, ax=ax, shrink=0.8, pad=0.03)
    cb.set_label(cfg['cbar'], fontsize=13, fontweight='bold')
    ax.set_title(f"Wafer Map: {wafer} / {band} ({cfg['maptag']})\nDate: {qdate}",
                 fontsize=17, fontweight='bold', pad=20)
    _show_fig(fig)


def plot_box(data, length_map, wafer, date, band, metric='vpil', selected_c=None, selected_r=None):
    # res 폴더의 Center vs Edge 박스플롯과 동일한 스타일.
    cfg = _MAP_CFG[metric]
    df = _metric_df(data, length_map, metric)
    qdate = _merge_date(date)

    # 범례가 우측으로 빠지므로 가로 비율을 살짝 늘려주면(9, 8) 균형이 맞습니다.
    fig, ax = plt.subplots(figsize=(9, 8))

    grp = df[(df['Wafer'] == wafer) & (df['Band'] == band) & (df['Date'] == qdate)] if not df.empty else df
    if df.empty or grp.empty:
        ax.text(0.5, 0.5, '데이터 없음', ha='center', va='center')
        _show_fig(fig)
        return
    tgt = cfg['target'][band] if isinstance(cfg['target'], dict) else cfg['target']
    c_data = grp[grp['Region'] == 'Center']['val'].values
    e_data = grp[grp['Region'] == 'Edge']['val'].values
    pos, box_data, labels, colors = [], [], [], []
    if len(c_data) > 0:
        pos.append(1);
        box_data.append(c_data);
        labels.append(f'Center\nn={len(c_data)}');
        colors.append('#3498db')
    if len(e_data) > 0:
        pos.append(2);
        box_data.append(e_data);
        labels.append(f'Edge\nn={len(e_data)}');
        colors.append('#e74c3c')
    if not box_data:
        ax.text(0.5, 0.5, '데이터 없음', ha='center', va='center');
        _show_fig(fig);
        return

    allv = np.concatenate(box_data)
    y_min = min(allv.min(), tgt) - 0.2
    y_max = max(allv.max(), tgt) + 0.2
    ax.set_ylim(y_min, y_max)

    # Good/Poor 영역 (VpiL은 작을수록 좋음, IL/ER은 클수록 좋음)
    if cfg['good_above']:
        ax.axhspan(tgt, y_max, facecolor='#e8f8f5', alpha=0.6, zorder=0, label='Good Region')
        ax.axhspan(y_min, tgt, facecolor='#fadbd8', alpha=0.6, zorder=0, label='Poor Region')
    else:
        ax.axhspan(y_min, tgt, facecolor='#e8f8f5', alpha=0.6, zorder=0, label='Good Region')
        ax.axhspan(tgt, y_max, facecolor='#fadbd8', alpha=0.6, zorder=0, label='Poor Region')

    box = ax.boxplot(box_data, positions=pos, patch_artist=True, widths=0.5,
                     flierprops=dict(marker='d', markerfacecolor='black', markersize=6, alpha=0.6), zorder=2)
    for p, cl in zip(box['boxes'], colors):
        p.set_facecolor(cl);
        p.set_alpha(0.7)
    for p, arr in zip(pos, box_data):
        ax.scatter(np.random.normal(p, 0.05, len(arr)), arr, color='black', alpha=0.5, s=20, zorder=3)

    avg = grp['val'].mean()
    ax.axhline(avg, color='blue', ls='--', lw=2.5, label=f'Avg: {avg:.2f} {cfg["unit"]}', zorder=4)
    ax.axhline(tgt, color='red', ls='-', lw=2.5, label=f'Target: {tgt:.2f} {cfg["unit"]}', zorder=4)

    # 🌟 [수정 포인트 1] 가로지르는 긴 선을 없애고 짧은 지역 선과 텍스트 위치를 최적화했습니다.
    if selected_c is not None and selected_r is not None:
        target_die = grp[(grp['Column'] == selected_c) & (grp['Row'] == selected_r)]
        if not target_die.empty:
            die_val = target_die['val'].values[0]
            die_region = target_die['Region'].values[0]
            x_pos = 1 if die_region == 'Center' else 2

            # 전체를 관통하지 않고 해당 박스 기둥 너비만큼만 짧은 가로 점선을 그려줍니다.
            ax.hlines(die_val, xmin=x_pos - 0.28, xmax=x_pos + 0.28, color='#FF007F', ls='--', lw=2, zorder=5)

            # 왕별 표시
            ax.scatter(x_pos, die_val, color='#FF007F', marker='*', s=350, edgecolor='black',
                       linewidths=1.5, zorder=10, label=f'Die ({selected_c}, {selected_r})')

            # 별 바로 오른쪽에 수치 텍스트 배치
            ax.text(x_pos + 0.15, die_val, f"{die_val:.2f}", color='#FF007F', fontweight='bold', va='center', zorder=11)

    ax.set_title(f"{cfg['boxtag']} Analysis : {wafer} ({band})\nDate: {qdate}", fontsize=18, fontweight='bold', pad=15)
    ax.set_xticks(pos);
    ax.set_xticklabels(labels, fontsize=14, fontweight='bold')
    ax.set_ylabel(cfg['cbar'], fontsize=16, fontweight='bold')

    # 🌟 [수정 포인트 2] 범례를 그래프 우측 바깥(bbox_to_anchor)으로 완전히 쫓아내어 겹침을 방지합니다.
    ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1), prop={'size': 11, 'weight': 'bold'})

    ax.grid(True, axis='y', ls=':', alpha=0.4, zorder=1)
    ax.set_xlim(0.5, max(pos) + 0.5)

    # 바깥으로 나간 범례가 잘리지 않도록 레이아웃 여백을 자동 조정합니다.
    fig.tight_layout()
    _show_fig(fig)

def interactive(data, length_map=None, metric='vpil'):
    # 드롭다운으로 다이를 골라 그 다이의 그래프를 화면에 띄운다 (노트북 전용).
    # metric: 'vpil' 이면 VpiL 곡선, 'phase' 면 위상(라디안) 곡선.
    import ipywidgets as widgets
    from ipywidgets import interact

    opts = sorted({(d['wafer_id'], d['date'], d['band'], d['die_c'], d['die_r']) for d in data})
    labels = [(f"{w} / {dt} / {b} ({c},{r})", (w, dt, b, c, r)) for (w, dt, b, c, r) in opts]

    def _show(sel):
        w, dt, b, c, r = sel
        if metric == 'phase':
            plot_die_phase(data, w, dt, b, c, r)
        else:
            plot_die_vpil(data, length_map, w, dt, b, c, r)

    interact(_show, sel=widgets.Dropdown(options=labels, description='다이:'))


def plot_die_summary(data, length_map, wafer, date, band, c, r):
    # 한 다이의 주요 분석을 한 그림(2x2)에 모아서 보여준다: raw / 평탄화 / VpiL / 위상
    d = get_die(data, wafer, date, band, c, r)
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    fig.suptitle(f'{wafer} {band} {date} ({c},{r})', fontsize=14, fontweight='bold')
    if d is None:
        axes[0, 0].text(0.5, 0.5, '데이터 없음', ha='center', va='center')
        _show_fig(fig)
        return

    m = (d['ref_data']['wl'] >= d['wl_min']) & (d['ref_data']['wl'] <= d['wl_max'])
    rwl, ril = d['ref_data']['wl'][m], d['ref_data']['il'][m]
    poly = ref_poly(rwl, ril, smooth=True) if len(rwl) >= 31 else None

    # (1) Raw 스펙트럼 (전압별)
    ax = axes[0, 0]
    for b in d['bias_data_list']:
        mb = (b['wl'] >= d['wl_min']) & (b['wl'] <= d['wl_max'])
        ax.plot(b['wl'][mb], b['il'][mb], lw=0.8, label=b.get('label', ''))
    ax.set_title('Raw spectra'); ax.set_xlabel('Wavelength (nm)'); ax.set_ylabel('IL (dB)')
    ax.grid(alpha=0.3); ax.legend(fontsize=7, ncol=2)

    # (2) 평탄화 스펙트럼
    ax = axes[0, 1]
    if poly is not None:
        for b in d['bias_data_list']:
            mb = (b['wl'] >= d['wl_min']) & (b['wl'] <= d['wl_max'])
            ax.plot(b['wl'][mb], b['il'][mb] - poly(b['wl'][mb]), lw=0.8)
    ax.set_title('Flattened spectra'); ax.set_xlabel('Wavelength (nm)'); ax.grid(alpha=0.3)

    # (3) VpiL 곡선
    ax = axes[1, 0]
    res = extract_vpil(d, (length_map or {}).get((wafer, c, r, band), DEFAULT_L_UM) * 1e-4)
    if res:
        vpil, volts, curve = res
        ax.plot(volts, curve, 's-', color='gray'); ax.plot(0.0, vpil, 'r*', markersize=14)
        ax.set_title(f'VpiL @0V = {vpil:.3f} V*cm')
    else:
        ax.set_title('VpiL: 유효값 없음')
    ax.set_xlabel('Voltage (V)'); ax.set_ylabel('Vpi*L (V*cm)'); ax.grid(alpha=0.3)

    # (4) 위상(라디안) 곡선
    ax = axes[1, 1]
    pr = extract_phase(d)
    if pr:
        V, PH, fsr = pr
        ax.plot(V, PH, 'o-', color='steelblue'); ax.axhline(np.pi, color='red', ls=':')
        ax.text(V.max(), np.pi, '  π', color='red', va='center')
        ax.set_title(f'Phase (rad)  FSR={fsr:.2f} nm')
    else:
        ax.set_title('Phase: 유효값 없음')
    ax.set_xlabel('Voltage (V)'); ax.set_ylabel('Δφ (rad)'); ax.grid(alpha=0.3)

    fig.tight_layout()
    _show_fig(fig)


def dashboard(data, length_map=None):
    # 버튼으로 웨이퍼/밴드/지표를 고르고 다이를 골라 그래프를 보는 UI (노트북 전용).
    import ipywidgets as widgets
    from IPython.display import display, clear_output

    wafer_b = widgets.ToggleButtons(options=list_wafers(data), description='웨이퍼')
    band_b = widgets.ToggleButtons(options=['LMZC', 'LMZO'], description='밴드')
    metric_b = widgets.ToggleButtons(options=[('VpiL', 'vpil'), ('IL', 'il'), ('ER', 'er')], description='지표')
    die_d = widgets.Dropdown(description='다이')
    out = widgets.Output()

    def _fill_dies():
        try:
            die_d.unobserve(_refresh, 'value')
        except ValueError:
            pass
        dies = sorted({(d['date'], d['die_c'], d['die_r']) for d in data
                       if d['wafer_id'] == wafer_b.value and d['band'] == band_b.value})
        die_d.options = [(f"{dt} ({c},{r})", (dt, c, r)) for (dt, c, r) in dies]
        die_d.observe(_refresh, 'value')

    def _refresh(*_):
        with out:
            clear_output(wait=True)
            if not die_d.value:
                print('이 웨이퍼/밴드 조합에는 다이가 없습니다.')
                return
            dt, c, r = die_d.value
            w, b, met = wafer_b.value, band_b.value, metric_b.value
            plot_die_summary(data, length_map, w, dt, b, c, r)
            plot_wafer_map(data, length_map, w, dt, b, met)
            plot_box(data, length_map, w, dt, b, met)

    def _on_top(*_):
        _fill_dies()
        _refresh()

    wafer_b.observe(_on_top, 'value')
    band_b.observe(_on_top, 'value')
    metric_b.observe(_refresh, 'value')

    _fill_dies()
    display(widgets.VBox([wafer_b, band_b, metric_b, die_d, out]))
    _refresh()
