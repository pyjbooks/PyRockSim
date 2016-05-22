# PyRockSim

## 概要
Rocket Launch Simulater written in Python
本プログラムは、https://github.com/ina111/MatRockSim をPythonに移植したものです。ただしロケットの打ち上げに関する機能のコア部分しか移植しておりません。枠組みとしては6自由度のロケット打ち上げシミュレーターですが、現在は実質的に姿勢は固定になっているので3自由度の質点モデルシミュレータとなっています。水平座標系での飛翔をシミュレーションしますが、プロット時にはきちんとした高度（楕円体高）に変換することもできます。

Python 3系を使う場合は本レポジトリの master ブランチを、Python 2系を使う場合は ver2 ブランチを使ってください。

## 必要なパッケージ
Python 2系の場合 : Version 2.7以上
Python 3系の場合 : Version 3.5以上
SciPy 0.17.0以上 (interp1dの外挿機能を利用するため)
NumPy
Matplotlit
basemap
（バージョン指定の無いものは適宜依存関係を満たすバージョン）

## 必要なパッケージのインストール
パッケージをひとつひとつインストールするのは手間なので、ディストリビューションパッケージの利用をお勧めします。以下に代表的なディストリビューションパッケージを示します。32bit/64bitの両プラットフォーム、およびWindows/OS X/Linuxの全てに対応しているAnacondaは特にお薦めです。

- Anaconda : https://www.continuum.io/
- Enthought Canopy : https://www.enthought.com/products/canopy/
- WinPython : http://winpython.github.io/
- Python(x,y) : http://python-xy.github.io/

なお、basemapというパッケージは自分で追加しない限りインストールされていないことが多いかと思います。 2016年5月の時点では、 Linux/OS X用のパッケージがAnacondaに用意されています。LinuxまたはOS Xを利用している方は、

 $ conda install basemap

とすればインストール可能です。

また、 Windowsの場合、Christoph Gohlkeが公開しているWindows用パッケージ（ http://www.lfd.uci.edu/~gohlke/pythonlibs/ ）を使えばほぼ問題なく使えると思われます。Christoph GohlkeはWheel形式のパッケージを配布していますので、Pythonのバージョンに対応するWheelをダウンロードして、例えば以下のようにインストールします。

 $ pip install basemap-1.0.8-cp35-none-win32.whl

## 実行
IPythonコンソール（シェル）では以下のように実行してください。

> In [1]: %run rocket.py

もしくはスクリプトモードで実行するなら以下のようにしてください。

> C:\python\PyRockSim> python rocket.py


## パラメータ設定
- 各種設定値はrocket.pyファイルの中の辞書変数 rocket_settings を変更してください。
- sp.integrate.odeintを使う積分では引数のatolおよびrtolを変更して精度と計算時間を調整できます。

## 座標系について
MatRockSimでは局地水平座標系（Local tangent frame）をUEN（Up-East-North）座標で表していますが、本プログラムではNED（North-East-Up）座標系にしてあります。それ以外はほぼMatRockSimと同じですのでMatRockSimの"Matlab Rocket Flight Simulator.pdf"を参照してください。

## ソースコードの文字コード
UTF-8としています。したがって、Python 3用のmasterブランチでは、ファイルの1行目にエンコーディング宣言は書かれていません。

## Future Works
- 特に今は考えていません。バグが見つかればFixします。

Copyright (C) 2016, Kenji Nakakuki
Released under MIT License
