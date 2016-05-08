"""
Quarternion Library

Copyright (c) 2016 Kenji Nakakuki
Released under the MIT license
"""

import numpy as np
from numpy import sin, cos
import unittest

D2R = np.pi / 180.0
R2D = 180.0 / np.pi


# モジュール関数の定義
def deltaquat(q, omega):
    """ クォータニオンの時間微分

    引数 :
        q : 微分を計算するクォータニオン (1x4)
        omega : Bodyの角速度 [rad/s]
    返り値 :
        dq/dt = - 1/2 * [0, omega] * q
    """
    return -0.5 * quatmultiply(np.r_[0, omega], q)


def quatnorm(q):
    """クォータニオンのノルム

    引数:
        q: クォータニオン (1x4)
    返り値:
        クォータニオンのノルム
    """
    return np.linalg.norm(q)


def quatnormalize(q):
    """クォータニオンの正規化

    引数:
        q: クォータニオン (1x4)
    返り値:
        正規化後クォータニオン
    """
    return q / (quatnorm(q))


def quatconj(q):
    """ 共役クォータニオン
    引数:
        q: クォータニオン (1x4)
    返り値:
        共役クォータニオン (1x4)
    """
    return [q[0], -q[1], -q[2], -q[3]]


def quatinv(p):
    """ 逆クォータニオン

    引数 :
        p: クォータニオン (1x4)
    返り値 :
        逆クォータニオン (1x4)
    """
    return quatconj(p) / (quatnorm(p) ** 2)


def quatmultiply(q, p):
    """ クォータニオン積

    引数 :
        q: クォータニオン (1x4)
        p: クォータニオン (1x4)
    返り値 :
        o: クォータニオン積 (1x4) 正規化してからリターンする
    """

    o1 = q[0] * p[0] - q[1] * p[1] - q[2] * p[2] - q[3] * p[3]
    o2 = q[1] * p[0] + q[0] * p[1] - q[3] * p[2] + q[2] * p[3]
    o3 = q[2] * p[0] + q[3] * p[1] + q[0] * p[2] - q[1] * p[3]
    o4 = q[3] * p[0] - q[2] * p[1] + q[1] * p[2] + q[0] * p[3]

    o = np.array([o1, o2, o3, o4])
    dq = 1 - np.sum(o * o)

    return o * (1 + 0.5 * dq)


def rot_coord_mat(q):
    """方向余弦行列をクォータニオンから計算する

    機体座標系のベクトルをNED local tangent frameに変換する際の回転行列に相当

    引数
        q: クォータニオン (1x4)
    返り値 :
        dcm: 方向余弦行列
    """
    dcm = [[1 - 2 * (q[2] ** 2 + q[3] ** 2), 2 * (q[1] * q[2] - q[0] * q[3]),
            2 * (q[1] * q[3] + q[0] * q[2])],
           [2 * (q[1] * q[2] + q[0] * q[3]), 1 - 2 * (q[1] ** 2 + q[3] ** 2),
            2 * (q[2] * q[3] - q[0] * q[1])],
           [2 * (q[1] * q[3] - q[0] * q[2]), 2 * (q[2] * q[3] + q[0] * q[1]),
            1 - 2 * (q[1] ** 2 + q[2] ** 2)]]
    return np.array(dcm)


def attitude(roll, pitch, yaw):
    """初期の方位角と仰角[deg]からクォータニオンと方向余弦行列を生成

    引数 :
        roll: ロール角 [deg]
        pitch: ピッチ角 [deg]
        yaw: ヨー角 [deg]
        備考) NED local tangent frameを yaw-->pitch-->rollの順に回して機体座標系を得る
         R(y) = [[cos_y, -sin_y, 0.], [sin_y, cos_y, 0.], [0., 0., 1.]]
         R(p) = [[cos_p, 0., sin_p], [0., 1., 0.], [-sin_p, 0., cos_p]]
         R(r) = [[1., 0., 0.], [0, cos_r, -sin_r], [0., sin_r, cos_r]]
         dcm_bn = R(y)*R(p)*R(r)
    返り値 :
        quat: クォータニオン[-] (4x1)
        dcm_bn : 方向余弦行列
    """

    cos_y, sin_y = cos(yaw * D2R), sin(yaw * D2R)
    cos_p, sin_p = cos(pitch * D2R), sin(pitch * D2R)
    cos_r, sin_r = cos(roll * D2R), sin(roll * D2R)
    rot_y = np.array([[cos_y, -sin_y, 0.], [sin_y, cos_y, 0.], [0., 0., 1.]])
    rot_p = np.array([[cos_p, 0., sin_p], [0., 1., 0.], [-sin_p, 0., cos_p]])
    rot_r = np.array([[1., 0., 0.], [0, cos_r, -sin_r], [0., sin_r, cos_r]])
    dcm_bn = rot_y @ rot_p @ rot_r

    c11, c12, c13 = dcm_bn[0, 0], dcm_bn[0, 1], dcm_bn[0, 2]
    c21, c22, c23 = dcm_bn[1, 0], dcm_bn[1, 1], dcm_bn[1, 2]
    c31, c32, c33 = dcm_bn[2, 0], dcm_bn[2, 1], dcm_bn[2, 2]
    tr = c11 + c22 + c33
    q1 = 0.5 * np.sqrt(1 + tr)
    k2 = np.sqrt(0.25 * (2 * c11 + 1 - tr))
    k3 = np.sqrt(0.25 * (2 * c22 + 1 - tr))
    k4 = np.sqrt(0.25 * (2 * c33 + 1 - tr))

    if k2 >= k3 and k2 >= k4:
        q2 = k2 * np.sign(c32 - c23)
        q3 = k3 * np.sign(q2 * (c21 + c12))
        q4 = k4 * np.sign(q2 * (c31 + c13))
    elif k3 >= k2 and k3 >= k4:
        q3 = k3 * np.sign(c13 - c31)
        q2 = k2 * np.sign(q3 * (c12 + c21))
        q4 = k4 * np.sign(q3 * (c32 + c23))
    else:
        q4 = k4 * np.sign(c21 - c12)
        q2 = k2 * np.sign(q4 * (c13 + c31))
        q3 = k3 * np.sign(q4 * (c23 + c32))

    quat = np.array([q1, q2, q3, q4])

    return quat, dcm_bn


class TestQuarternion(unittest.TestCase):
    """
    本モジュールのテストコード（for unittest）
    unittest.TestCaseのサブクラスとすることで自動的にunittestの
    テストケースであると認識される。
    """

    def setUp(self):
        """ 各テストメソッドの前にSet upとして実行するコード """
        self.q = [0.499524110790929, 0.865201139495554,
                  0.021809693682668, -0.037775497555895]
        self.p = [0.999942884777492, 0.00872648011789139,
                  0.0043632400589457, 0.0043632400589457]

    def tearDown(self):
        """ 各テストメソッドの最後に実行するコード """
        print('tearDown: テストメソッドを1つ実行しました')

    def test_quatconj(self):  # テストメソッドは test_ で始まる名前とする
        q_ans = quatconj(self.q)
        self.assertEqual(self.q[0], q_ans[0])
        for i in range(3):
            self.assertEqual(self.q[i + 1], -q_ans[i + 1])

    def test_quatmultiply(self):
        q_mul = quatmultiply(self.q, self.p)
        q_ans = [0.49201508245344, 0.869770795054512,
                 0.0198832642285152, -0.0320090379767411]
        for i in range(4):
            self.assertAlmostEqual(q_mul[i], q_ans[i], places=14)

    def test_deltaquat(self):
        omega = [0.017453292519943295, 0.008726646259971648,
                 0.008726646259971648]
        dq = deltaquat(self.q, omega)
        dq_ans = [0.0112192514401438, -0.0061478346436931,
                  -0.0094251502315234, 0.002107541287674]
        for i in range(4):
            self.assertAlmostEqual(dq[i], dq_ans[i], places=14)


# 本モジュールをスクリプトとして実行した場合に実行されるコード ＝ unittestによる単体試験
if __name__ == "__main__":
    unittest.main()
