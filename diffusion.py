import numpy as np
import pandas as pd
import numba
import matplotlib.pyplot as plt

def make_input_param(index_number):
    input_params = pd.read_csv("input_param.csv")
    list_element = input_params["Element"].values
    list_tempc = input_params["Temp(C)"].values
    list_temp_sd = input_params["TempSD(C)"].values
    list_pressmpa = input_params["Press(MPa)"].values
    list_oxidebuffer = input_params["Oxide buffer"].values
    list_dx = input_params["dx(um)"].values
    list_maxtime = input_params["MaxTime(yr)"].values
    list_orientation = input_params["Orientation(degree)"].values
    element = list_element[index_number]
    tempc = list_tempc[index_number]
    temp_sd = list_temp_sd[index_number]
    pressmpa = list_pressmpa[index_number]
    oxbuffer = list_oxidebuffer[index_number]
    dx = list_dx[index_number]
    maxtime = list_maxtime[index_number]
    # カンマで分割し，numpy array形式に変換
    orientation = list(map(float, list_orientation[index_number].split(',')))
    orientation = np.array(orientation)
    return element, tempc, temp_sd, pressmpa, oxbuffer, orientation, dx, maxtime

# 組成valueを0<=value<=1で規格化
def norm_value(value):
    normalized_value = (value - min(value)) / (max(value) - min(value))
    return normalized_value

# 規格化したものをもとに戻す
def restore_value(normalized_value, unit_value):
    restored_value = normalized_value * (max(unit_value) - min(unit_value)) \
                     + min(unit_value)
    return restored_value

# 数値計算
# 無次元化しているものとする
# 陽的差分法の精度を確保するためにr=1/6に固定
@numba.jit(nopython=True)
def calc_fdm(norm_initial_value, nt, r=(1/6)):
    nx = len(norm_initial_value)
    u = norm_initial_value
    u_next = np.zeros(nx)
    # 計算結果を格納する空の配列を事前に作成
    # 0行目に初期条件を入れるため，行数は+1している
    array_result = np.zeros((nt+1, nx))
    array_result[0] = u
    # 時間ステップを逐次計算
    for i in range(nt):
        # 空間ステップを逐次計算
        for j in range(nx):
            if j == 0:
                # リム側はDirichlet境界条件(組成一定)を仮定
                u_next[j] = u[j]
            elif j == nx-1:
                # コア側は対称境界条件を仮定
                u_next[j] = u[j] + r * (u[j] - 2 * u[j] + u[j-1])
            else:
                u_next[j] = u[j] + r * (u[j+1] - 2 * u[j] + u[j-1])
        u = u_next
        # nt行nx列行列のi行目(時間方向のインデックス)に計算結果を代入
        array_result[i+1] = u
    return array_result

# 計算結果のcsvファイル("result.csv")を読み込んでデータフレームを返す
def make_inputfile(filename):
    df = pd.read_csv(filename)
    dist = df["Distance(um)"].values
    values = df["Value"].values
    data = [dist, values]
    return data, dist, values

# 計算値の2点間の直線 y=m*x+n のパラメータ m, n を計算する
def calc_param(initial_dist, result_values, dx):
    params = np.zeros([len(initial_dist)-1, 2])
    for j in range(len(initial_dist)-1):
        m = (result_values[j+1] - result_values[j]) / dx
        n = result_values[j] - j * dx
        params[j] = [m,n]
    return params


def calc_bestfitindex(array_result, initial_dist, measured_dist, measured_value, dx):
    list_ssr = []
    num_result = len(array_result)
    for i in range(num_result):
        result_values = array_result[i]
        params = calc_param(initial_dist, result_values, dx)
        # 残差二乗和の変数を作成
        ssr = 0
        for j in range(len(initial_dist)-1):
            # 分析値を順番に入力
            for k in range(len(measured_dist)):
                # もし範囲内にあれば内挿して残差二乗を計算
                # 末尾に来た場合を想定している
                begin_dist, end_dist = initial_dist[j], initial_dist[j+1]
                cf_meas_dist = measured_dist[k]
                if (begin_dist <= cf_meas_dist < end_dist ) \
                or (cf_meas_dist == end_dist):
                    lerp_result_value = params[j,0] * cf_meas_dist + params[j,1]
                    ssr += (lerp_result_value - measured_value[k]) ** 2
                else:
                    pass
        list_ssr.append(ssr)
    # 測定値を最も説明する値を決定
    best_fit_index = list_ssr.index(min(list_ssr))
    return best_fit_index

# プロット
def make_image(measured_data, initial_data, fit_data, tempc, time_yrs):
    fig = plt.figure()
    plt.rcParams["font.size"] = 16
    ax = fig.add_subplot(111)
    ax.set_title(str(time_yrs)+" yr in " + str(tempc) + "$^\circ$C")
    ax.plot(*measured_data, "o", color="w", mec="k", label="Measured")
    ax.plot(*initial_data, "--", color="k", label="Initial")
    ax.plot(*fit_data, "-", color="r", label="Best fit")
    ax.set_xlabel("Distance (\u03bcm)", fontsize=18)
    ax.set_ylabel("Value", fontsize=18)
    fig.legend(loc=2, fancybox=False, framealpha=1, edgecolor="k", fontsize=14)
    fig.savefig('img.jpg', dpi=300, bbox_inches='tight')

# fO2計算クラス
class CalcOxBuffer():
    def __init__(self, tempk, press, oxbuffer):
        # 示強変数を格納
        self.tempk = tempk
        self.pressbar = press * 10 ** (-5)
        self.oxbuffer = oxbuffer
        # ref Frost (1991) Rev Mineral Geochem
        self.COEFS = {"FMQ": (-25096.3, 8.735, 0.110),
                      "NNO": (-24930, 9.36, 0.046),
                      "MH": (-25700.6, 14.558, 0.019)}
    # calc fO2 buffer
    def calc_fo2(self, ):
        coefs = self.COEFS[self.oxbuffer]
        log_fo2 = (coefs[0] / self.tempk) + coefs[1] \
                  + coefs[2] * (self.pressbar-1) / self.tempk
        # print("fO2: 10^", log_fo2)
        return 10 ** log_fo2

# 拡散係数クラス
class SetCoef():
    # 初期化メソッド
    def __init__(self, tempk, press, oxbuffer):
        self.tempk = tempk
        self.press = press
        self.fo2 = CalcOxBuffer(self.tempk, self.press, oxbuffer).calc_fo2()
        self.R_CONST = 8.31

    # 斜長石CaAl-NaSi (Liu and Yund, Am Mineral 1992)
    def calc_coef_plg_an(self, ):
        if self.tempk >= 985 + 273.15:
            # 1000-1050C
            coef = 4 * (10 ** (-16)) \
                   * np.exp(-103000 / (self.R_CONST * self.tempk))
        else:
            # 900-975C
            coef = 11 * (10 ** (-16)) \
                   * np.exp(-371000 / (self.R_CONST * self.tempk))
        return coef

    # 直方輝石Fe-Mg (Dohmen et al., Am Mineral 2016)
    def calc_coef_opx_femg(self, orientation=[0,90,90]):
        # [001]方向の拡散は[100]方向の3.5倍速いので係数を追加
        coef_c = 1.12 * (10 ** (-6)) * (self.fo2 ** 0.053) \
                 * np.exp(-308000 / (self.R_CONST*self.tempk))
        coef_a = (1 / 3.5) * coef_c
        coef_b = coef_a
        coefs = (coef_a, coef_b, coef_c)
        coef = calc_orientation(coefs, orientation)
        return coef
    
    # 単斜輝石Fe-Mg (Muller et al., Contrib Mineral Petrol)
    def calc_coef_cpx_femg(self, ):
        coef = 2.77 * (10 ** (-7)) \
               * np.exp(-320700 / (self.R_CONST*self.tempk))
        return coef

    # カンラン石Fe-Mg (Dohmen and Chakraborty, Phys Chem Miner 2007)
    def calc_coef_olv_femg(self, xfe, orientation=[0,90,90]):
        coef_c = (10 ** (-9.21)) * ((self.fo2 / (10 ** (-7))) ** (1 / 6))\
                 * (10 ** (3 * (xfe - 0.1))) * np.exp((-201000
                 + (self.press - (10 ** 5)) * (7 * (10 ** (-6))))
                 / (self.R_CONST * self.tempk))
        # [001]方向の拡散は[100]方向の6倍速いので係数を追加
        coef_a = (1 / 6) * coef_c
        coef_b = coef_a
        coefs = (coef_a, coef_b, coef_c)
        coef = calc_orientation(coefs, orientation)
        return coef

# 結晶軸の方位の影響を計算
def calc_orientation(coefs, orientation):
    coef_a, coef_b, coef_c = coefs[0], coefs[1], coefs[2]
    rads = np.deg2rad(orientation)
    alpha, beta, gamma = rads[0], rads[1], rads[2]
    coef = coef_a * (np.cos(alpha) ** 2) + coef_b * (np.cos(beta) ** 2) \
           + coef_c * (np.cos(gamma) ** 2)
    return coef

# 単位換算
def convert_units(tempc, pressmpa, oxbuffer):
    tempk = tempc + 273.15
    press = pressmpa * (10 ** 6)
    extensive_vars = (tempk, press, oxbuffer)
    return extensive_vars

# 計算結果から無次元化した時間をもとに戻す
# dt'/dx'^2 = 1/6, t = τ*t' = τ*nt*dt', x = xmax*x' (∵ 0<=x'max<=1)
# D*τ/xmax^2 = 1より，τ = xmax^2/D, またdt' = dx'^2/6
# したがって，t = (xmax^2/D)*nt*(dx'^2/6)
# 各単位は揃える(SI単位系として計算)
class ConvertTime:
    def __init__(self, initial_dist, coef):
        coef = coef
        xmax = max(initial_dist) * 10 ** (-6)
        norm_dx = (1 / (len(initial_dist) - 1))
        self.rr = ((xmax ** 2) / coef) * ((norm_dx ** 2) / 6)
    # 測定値から実時間を計算
    def restore_time(self, best_fit_index):
        time_s = best_fit_index * self.rr
        # 小数点以下を切り捨てていることに注意
        time_yrs = round(time_s / (60 * 60 * 24 * 365))
        return time_yrs
    # 入力値から時間ステップを計算
    def calc_nt(self, maxtime_yr):
        maxtime_sec = 60 * 60 * 24 * 365 * maxtime_yr
        # ここの切り捨て操作の関係上，入力したMaxTimeは計算のMaxTimeと厳密には一致しない．
        nt = int(maxtime_sec / self.rr)
        print("Time step is",nt)
        return nt

def main():
    index_number = 0
    element, md_tempc, temp_sd, pressmpa, oxbuffer, orientation, dx, maxtime \
        = make_input_param(index_number)
    print("maxtime is", maxtime)
    # 初期組成，分析値を入力
    # 初期組成は分析値よりコア側を同じか幅広く設定する
    initial_data, initial_dist, initial_value \
        = make_inputfile("initial_value.csv")
    measured_data, measured_dist, measured_value \
        = make_inputfile("measured_value.csv")
    # 数値計算のために初期組成を無次元化
    norm_initial_value = norm_value(initial_value)
    # 拡散係数を決定
    lower_tempc = md_tempc - temp_sd
    upper_tempc = md_tempc + temp_sd
    list_tempc = np.array([lower_tempc, md_tempc, upper_tempc])
    list_time = []
    for tempc in list_tempc:
        print("Temp C is", tempc)
        extensive_vars = convert_units(tempc, pressmpa, oxbuffer)
        if element == "pl-CaAlNaSi":
            coef = SetCoef(*extensive_vars).calc_coef_plg_an()
        elif element == "opx-FeMg":
            coef = SetCoef(*extensive_vars).calc_coef_opx_femg(orientation)
        elif element == "cpx-FeMg":
            coef = SetCoef(*extensive_vars).calc_coef_cpx_femg()
        elif element == "olv-FeMg":
            xfe = initial_value[0]
            coef = SetCoef(*extensive_vars).calc_coef_olv_femg(xfe, orientation)
        else:
            print("Not Found!")
            pass
        # 入力値から時間ステップを計算
        convert_time = ConvertTime(initial_dist, coef)
        nt = convert_time.calc_nt(maxtime)
        # 陽的差分法による数値計算を実行し，配列に格納する
        array_result = calc_fdm(norm_initial_value, nt)
        # 測定値を最も説明する値を決定
        best_fit_index = calc_bestfitindex(array_result, initial_dist, \
                            measured_dist, measured_value, dx)
        # 計算結果から無次元化した時間をもとに戻す
        time_yrs = convert_time.restore_time(best_fit_index)
        if nt == best_fit_index:
            print("Input longer time!")
        list_time.append(time_yrs)
        # 無次元化していた組成をもとに戻し，配列に格納
        if tempc == list_tempc[1]:
            fit_result_value \
            = restore_value(array_result[best_fit_index], initial_value)
    
    # 図を出力
    list_time = np.array(list_time)
    # 温度が高いほうが時間は短くなるため
    md_time = list_time[1]
    lower_sd = list_time[2] - md_time
    upper_sd = list_time[0] - md_time
    print("Diffusion time is", md_time, \
          "lower SD:", lower_sd, "upper SD", upper_sd, "yrs.")
    
    fit_data = (initial_dist, fit_result_value)
    list_plot_data = (measured_data, initial_data, fit_data)
    make_image(*list_plot_data, md_tempc, md_time)

if __name__ == '__main__':
    main()