# diffusion-calculator
斑晶鉱物中の元素拡散をシミュレーションするためのプログラムです．

## Requirement
以下の環境で開発・動作を確認しています．
- Windows 10 Home 20H2
- Python 3.8.5
- Numpy 1.19.5
- Pandas 1.2.0
- Maptlotlib 3.3.3

## Features
斑晶鉱物の元素拡散クロノメトリーはマグマプロセスの時間スケールを決定する目的で広く使われていますが，計算がやや面倒なことが多いのがネックです．
オリビンの計算を支援するプログラムとしてはDIPRA (Girona and Costa, Geochem Geophys Geosys 2013) がありますが，直方輝石や斜長石向けの使いやすいプログラムは公開されていませんでした．
diffusion-calculatorは (1) 多くの種類の鉱物に対応すること，(2) 複雑な初期条件についても計算可能であること，(3) Excel等のGUIでの入力（とターミナル上でのプログラムの実行）だけで操作が完結することを目指して作成しました．

計算手法の詳細は以下の通りです．
- アルゴリズム：単純な陽的差分法で実装，ただしdt/dx^2=1/6となるように無次元化を行っているのでまあまあ精度は高いです．
- 境界条件：リムはメルトと平衡（Dirichlet条件）とし，コアについて対称なゾーニングをしているものと仮定し計算を行います（対称境界条件）．
- フィッティング：計算結果を線形に内挿した上で，最小二乗法を適用して最も一致する時間を出力します．あらかじめ分析値の外れ値は除外することをおすすめします．
- 誤差の計算：温度の誤差＝拡散時間の誤差とみなして計算しています．間接的に温度がfO2に及ぼす影響も加味されますが，酸素バッファーそのものの誤差は計算されません．
フィッティングの精度については一貫した計算が難しいため無視しています．

## Usage
1. `diffusion.py`を`initial_value.csv`と`measured_value.csv`, `input_param.csv`と同じディレクトリに設置します．
1. `initial_value.csv`に初期プロファイルを，`measured_value.csv`に分析値を入力します．
"Distance(um)"の列で0 μmの点がリムに，最終行がコア部分の組成となるようにします．
1. `imput_param.csv`に計算する元素，温度(C)，温度の標準誤差(C)，圧力(MPa)，酸素バッファ(FQM, NNO, MHのどれか)，分析線と各結晶軸のなす角度(°)，計算時間の最大値(年)を入力します．
結晶軸のなす角度は，`a軸,b軸,c軸`の順で**半角カンマ，後ろに空白はつけずに入力してください**．デフォルトではa軸に平行に分析したものとして計算が行われます．
1. 計算する元素は現在時点で，斜長石CaAl-NaSi (`pl-CaAlNaSi`)，直方輝石Fe-Mg (`opx-FeMg`)，単斜輝石Fe-Mg (`cpx-FeMg`)，かんらん石Fe-Mg (`olv-FeMg`) が使用できます．
引用している拡散係数についてはReferenceも参照してください．
3. `diffusion.py`を実行します．計算中にターミナルには`maxtime`（最大計算時間），`Time step`（時間ステップの繰り返し計算回数，計算に要する時間の目安になります），`Diffusion time`（拡散時間）が表示されます．
4. 最大計算時間と拡散時間が一致した場合，`Input longer time!`と表示されます．最大計算時間が小さすぎる可能性があるのでもっと大きな値を入力してください．
5. 計算が終了すると同じディレクトリに計算結果とフィッティングした拡散プロファイルが描画された`img.jpg`が保存されます．

## Reference
- 拡散方程式の導出 (Crank, 1975)
- 酸素フガシティのパラメータ (Frost, Rev Mineral Geochem 1991)
- 拡散係数の計算
  - 斜長石CaAl-NaSi (Liu and Yund, Am Mineral 1992)
  - 直方輝石Fe-Mg (Dohmen et al., Am Mineral 2016)
  - 単斜輝石Fe-Mg (Muller et al., Contrib Mineral Petrol 2013)
  - カンラン石Fe-Mg (Dohmen and Chakraborty, Phys Chem Miner 2007)
- 拡散係数の異方性の計算 (Costa and Chakraborty, Earth Planet Sci Lett 2004)

## Author
* Motohiro Sato (佐藤初洋)
* E-mail: msatores "at" gmail.com

## License
Copyright (C) 2021 Motohiro Sato

diffusion-calculator is released under [MIT license](https://opensource.org/licenses/mit-license.php).
