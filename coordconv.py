"""
座標変換用の関数

Copyright (c) 2016 Kenji Nakakuki
Released under the MIT license
"""

from collections import namedtuple
import numpy as np
from numpy import sin, cos, sqrt, arctan2, pi

D2R = pi / 180.0
R2D = 180.0 / pi

# 定数の設定
WGS84 = namedtuple('WGS84', ['re_a',  # [m] WGS84の長軸
                             'eccen1',  # First Eccentricity
                             'eccen1sqr',  # First Eccentricity squared
                             'one_f',  # 扁平率fの1/f（平滑度）
                             're_b',  # [m] WGS84の短軸
                             'e2',  # 第一離心率eの2乗
                             'ed2'  # 第二離心率e'の2乗
                             ])
wgs84 = WGS84(6378137.0, 8.1819190842622e-2, 6.69437999014e-3, 298.257223563,
              6356752.314245, 6.6943799901414e-3, 6.739496742276486e-3)


def blh2ecef(lat, lon, height):
    """緯度／経度／高度からECEF座標に変換

    緯度／経度／高度(Latitude/Longitude/Height):LLH座標系から
    地球中心地球固定座標ECEF(Earth Centered Earth Fixed)座標系に変換する

    引数 :
        lat : WGS84 緯度 [deg]
        lon : WGS84 経度 [deg]
        height : WGS84 ジオイド高 [m]
    """

    s_lat = sin(lat * D2R)
    c_lat = cos(lat * D2R)
    s_lon = sin(lon * D2R)
    c_lon = cos(lon * D2R)

    re_n = wgs84.re_a / sqrt(1 - (wgs84.eccen1 * s_lat) ** 2)
    ecef_x = (re_n + height) * c_lat * c_lon
    ecef_y = (re_n + height) * c_lat * s_lon
    ecef_z = (re_n * (1 - wgs84.eccen1sqr) + height) * s_lat

    return [ecef_x, ecef_y, ecef_z]


def ecef2blh(x, y, z):
    """ECEF座標から緯度/経度/高度に変換

    引数 :
        x,y,z: ECEF座標での位置[m]

    返り値 :
        phi: 緯度[deg]
        lam: 経度[deg]
        height: WGS84の平均海面高度[m]
    """

    p = sqrt(x ** 2 + y ** 2)  # 現在位置での地心からの距離[m]
    theta = arctan2(z * wgs84.re_a, p * wgs84.re_b)  # [rad]

    phi = R2D * arctan2(z + wgs84.ed2 * wgs84.re_b * sin(theta) ** 3,
                        p - wgs84.e2 * wgs84.re_a * cos(theta) ** 3)
    lam = R2D * arctan2(y, x)
    height = p / cos(D2R * phi) - wgs84.re_a / \
        sqrt(1.0 - wgs84.e2 * sin(D2R * phi) ** 2)

    return [phi, lam, height]


def launch2ecef(n, e, d, xr, yr, zr):
    """ 射点座標系からECEF座標系へ座標変換

    引数 :
        n,e,d : 射点中心座標系のNorth-East-Down座標[m]
        xr,yr,zr : ECEF-XYZ座標上の参照位置（射点）:[m]

    返り値 :
        x,y,z : ECEF座標系上の座標[m]
    """

    # 射点の緯度経度
    phi, lam, _ = ecef2blh(xr, yr, zr)
    phi *= D2R
    lam *= D2R
    s_phi = sin(phi)
    c_phi = cos(phi)
    s_lam = sin(lam)
    c_lam = cos(lam)

    x = -s_phi * c_lam * n - s_lam * e - c_phi * c_lam * d + xr
    y = -s_phi * s_lam * n + c_lam * e - c_phi * s_lam * d + yr
    z = c_phi * n - s_phi * d + zr

    return [x, y, z]


def dcm_x2n(phi, lam):
    """WGS84 ECEF-XYZからLocal tangent NED直交座標系への回転行列を計算

    引数
        phi : 緯度 [deg]
        lam : 経度 [deg]

    返り値
        dcm : WGS84 ECEF-XYZからLocal tangent NED直交座標系への回転行列

    単体試験（doctestによる単体試験の例）
    >>> dcm = dcm_x2n(38.54, 140.123)  # 3x3の行列を得る
    >>> round(dcm[0, 0], 13) == 0.4781509665478
    True
    >>> round(dcm[0, 1], 13) == -0.3994702417770
    True
    >>> round(dcm[0, 2], 13) == 0.7821733689688
    True
    >>> round(dcm[1, 0], 13) == -0.6411416200655
    True
    >>> round(dcm[1, 1], 13) == -0.7674225843822
    True
    >>> round(dcm[1, 2], 13) == 0.0
    True
    >>> round(dcm[2, 0], 13) == 0.6002575082490
    True
    >>> round(dcm[2, 1], 13) == -0.5014839009527
    True
    >>> round(dcm[2, 2], 13) == -0.6230608484538
    True
    """

    sphi, cphi = sin(phi * D2R), cos(phi * D2R)
    slam, clam = sin(lam * D2R), cos(lam * D2R)
    dcm = [[-sphi * clam, -sphi * slam, cphi],
           [-slam, clam, 0.0],
           [-cphi * clam, -cphi * slam, -sphi]]
    return np.array(dcm)


if __name__ == "__main__":
    import doctest

    doctest.testmod()
