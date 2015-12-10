#! coding:utf-8
"""
cepstrum_for_gwt
GWT配列からケプストラムを算出


Created by fifi  (2015/12/02 6:17)
"""
__version__ = '0.1'

import numpy as np


def cepstrum_for_gwt(gwtdata=None, lifter1=5, lifter2=None, type='det'):
    """
    GWTのデータ配列からケプストラムを算出
    dataは1次元配列のみ可能.
    ※ ケプストラムには周波数軸上での折り返し成分が必要.
    ※ GWTスペクトルには折り返し成分がないため、人為的に作成する。
    :param data: 振幅スペクトル 1D-ndarray
    :param lift: リフタリング次数
    :param lift2: リフタリング次数. type='band'のときのみ利用。
    :param type: 'low', 'high', 'band'
    :return: lifterd_data
    """

    # 左側データ
    gwt_left = np.atleast_2d(gwtdata)
    print gwt_left.shape

    # 右側データの生成：(人為的に折り返す)
    gwt_right = np.fliplr(gwt_left)
    print gwt_right.shape

    # numpy.hstack() で列方向に結合
    gwt_data = np.hstack((gwt_left, gwt_right))
    print gwt_data.shape

    del gwt_right
    del gwt_left
    del gwtdata

    # 周波数領域
    FFT_abs = np.abs(gwt_data)

    # ケプストラム
    _absSpec = FFT_abs
    ceps = np.real(np.fft.ifft(np.log(_absSpec)))

    # リフタリング
    # ローパス：'low'
    rown = ceps.shape[0]
    cols = ceps.shape[1]
    Ceps_low = ceps.copy()
    print Ceps_low
    # Ceps_low[:, lifter1:Ceps_low.shape[1]-lifter1] = 0
    # w = np.vstack((np.ones(rown, lifter1), np.zeros(rown, cols - 2*lifter1), np.ones(rown, lifter1)))
    w = np.hstack((np.ones((rown, lifter1)), np.zeros((rown, cols - 2 * lifter1)), np.ones((rown, lifter1))))
    print w.shape
    Ceps_low = Ceps_low * w
    print Ceps_low.shape
    print Ceps_low

    # 逆毛ぷす
    Ceps_sspec = np.abs(np.exp(np.fft.fft(Ceps_low)))

    # 反転成分除去
    Cepsed = Ceps_sspec[:, cols / 2:]
    return Cepsed

    pass


def main():
    import matplotlib.pyplot as plt
    Fs = 400
    f = 5
    t = np.linspace(0, 1, Fs)
    A = 100
    Y = 0
    for f in [5, 25, 50, 52, 55, 110, 120]:
        Y += A * np.sin(2 * np.pi * f * t)
    Y = np.array([Y])

    F = np.abs(np.fft.fft(Y))
    F = F[:F.size/2].copy()

    # ケプストラム
    lifted = cepstrum_for_gwt(gwtdata=F.copy(), lifter1=15, type='low')

    print Y.shape
    print lifted.shape
    # plt.plot(t, F[0])
    # plt.plot(t, lifted[0])
    plt.plot(t, 20*np.log10(F[0]))
    plt.plot(t, 20*np.log10(lifted[0]))
    plt.show()

    pass


if __name__ == "__main__":
    main()
