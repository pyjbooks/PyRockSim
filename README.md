# PyRockSim

## 概要
Rocket Launch Simulater written in Python
本プログラムは、https://github.com/ina111/MatRockSim をPythonに移植したものです。
ただしロケットの打ち上げに関する機能のコア部分しか移植しておりません。
枠組みとしては6自由度のロケット打ち上げシミュレーターですが、現在は実質的に姿勢は固定になっているので
3自由度の質点モデルシミュレータとなっています。水平座標系での飛翔をシミュレーションしますが、プロット時には
きちんとした高度（楕円体高）に変換することもできます。

## 実行
Python 3.5の機能を使って書かれていますのでPython 3.5以降を使ってください。
とはいえ、行列の積の@演算子を使わなければPython 3.4などでも動きます。
必要なパッケージ類はimport文を見れば分かりますが、NumPy、SciPy、Matplotlib、Basemapが
必要になります。

IPythonコンソール（シェル）では以下のように実行してください。

> In [1]: %run rocket.py

もしくはスクリプトモードで実行するなら以下のようにしてください。

> C:\python\PyRockSim> python rocket.py


## パラメータ設定
- rocket.pyファイルの中の辞書変数 rocket_settings を変更してください。
- 抗力係数はRocketSimクラスのcd_rocketメソッドの中のテーブル（リスト）の値を変更してください。
- 方程式の積分精度はrocket.pyの中のatolおよびrtolを変更してください。

## 座標系について
MatRockSimでは局地水平座標系（Local tangent frame）をUEN（Up-East-North）座標で表していますが、
本プログラムではNED（North-East-Up）座標系にしてあります。
それ以外はほぼMatRockSimと同じですのでMatRockSimの"Matlab Rocket Flight Simulator.pdf"を参照して
ください。

## ソースコードの文字コード
UTF-8としています。したがって、ファイルの1行目にエンコーディング宣言は書かれていません。

## Future Works
- 特に今は考えていません。バグが見つかればFixします。

Copyright (C) 2016, Kenji Nakakuki
Released under MIT License
