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

なお、basemapというパッケージは自分で追加しない限りインストールされていないことが多いかと思います。 2017年2月の時点では、 Linux/OS X用のパッケージがAnaconda（Python 2.7/3.4/3.5/3.6の全てのバージョンでbasemap 1.0.7が使えます）に
用意されています。LinuxまたはOS Xを利用している方は、

    $ conda install basemap

とすればインストール可能です。

また、 Windowsの場合、Christoph Gohlkeが公開しているWindows用パッケージ( http://www.lfd.uci.edu/~gohlke/pythonlibs/ )を使えばほぼ問題なく使えると思われます。
Christoph GohlkeはWheel形式のパッケージを配布していますので、PythonのバージョンとOSの種類（32bit OS/64bit OS）に対応するWheelをダウンロードして、
例えば以下のようにインストールします(32bit版 WindowsでPython 3.5を使う場合)。

    $ pip install basemap-1.0.8-cp35-none-win32.whl

なお、2017年2月時点で最新のAnaconda 4.3.0(Python 3.6版)を使っている場合は、依存関係のあるパッケージも同時にインストールする必要（適切にC/C++コンパイラがインストールされて
おり、適切に環境設定されていれば、basemapのインストールをpipに指示するだけで、あとはpipが勝手に必要なパッケージをダウンロードして、コンパイル後にインストールしてくれます。
しかし、下記のように必要なコンパイル済みパッケージを全てダウンロードしてインストールする方がトラブルに遭遇する可能性は低くなります。）が
あるので、まずChristoph Gohlkeのパッケージから以下の3つのパッケージを取得します(64bit OSでPython 3.6を使っている場合の例)。

    pyproj-1.9.5.1-cp36-cp36m-win_amd64.whl
    pyshp-1.2.10-py2.py3-none-any.whl
    basemap-1.0.8-cp36-cp36m-win_amd64.whl

これらの3つのwheelパッケージを、pipを使ってインストールしていきます。
（以下、"Anaconda Prompt"(DOSコマンドプロンプト)のプロンプトを > で表します）

具体的には、3つのwheelパッケージを保存したフォルダにおいて、
以下の3つのコマンドを実行します。

    > pip install pyproj-1.9.5.1-cp36-cp36m-win_amd64.whl
    > pip install pyshp-1.2.10-py2.py3-none-any.whl
    > pip install basemap-1.0.8-cp36-cp36m-win_amd64.whl

以上でインストール終了です。32bit OSの場合や、Python 3.5などの場合の他、パッケージの更新によっても
wheelパッケージの名前が変わりますのでご注意ください。

ところで、Windowsユーザがbasemapをインストールする方法にはもう一つあります。
それは、conda（Anacondaのパッケージマネージャ）を使ってインストールする方法です。
先に述べたように、AnacondaはWindows用のbasemapパッケージを提供していませんが、conda-forge( https://conda-forge.github.io/ )
のconda用のパッケージレポジトリに、Windows用（32bit/64bit両方）の
パッケージが用意されているのです。

具体的には、以下の2つのインストールコマンドを実行します。「-c conda-forge」は、パッケージを「conda-forge」から取得してくる
ことを意味します。2つ目のコマンドでは、詳細な精度の高い地図データをインストールしています。

    > conda install -c conda-forge basemap=1.0.8.dev0
    > conda install -c conda-forge basemap-data-hires

これらを実行するだけで、インストール完了です。ただし、上記のコマンドで「1.0.8.dev0」とあるところは、
適宜最新のバージョンを指定するといいでしょう。最新のバージョンは、ここ( https://anaconda.org/conda-forge/basemap )で確認できます。

なお、basemapインストールの際に、conda自体のバージョンが古いものに変わってしまう場合があります。
他のパッケージの管理に影響すると考えられる場合は、次のようにしてcondaのバージョンを戻してしまいましょう。

    > conda update conda

これでcondaのバージョンを元に戻しても、basemapの利用上はなんら問題ありません。


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

## branches

- master : Python 3用のメインブランチ（Python 3.5以上）
- pdb_practice : 故意にバグを混入した、デバッグ練習用ブランチ（Python 3.5以上）
- ver2 : Python 2で動作するように修正したブランチ（Python 2.7以上）
- unittest : 主要な全関数の単体試験コードを、masterブランチに追加したブランチ（Python 3.5以上）

## ソースコードの文字コード
UTF-8としています。したがって、Python 3用のmasterブランチでは、ファイルの1行目にエンコーディング宣言は書かれていません。

## Future Works
- 特に今は考えていません。バグが見つかればFixします。

Copyright (C) 2016, Kenji Nakakuki
Released under MIT License
